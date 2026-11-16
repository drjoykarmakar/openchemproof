#!/usr/bin/env python3
"""Verify normalized rows against retained original source bytes, without a network."""

import hashlib
import json
import re
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

from rdkit import Chem

BASE = Path(__file__).resolve().parent
NS = {"c": "http://www.xml-cml.org/schema"}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    manifest = json.loads((BASE / "manifest.json").read_text())
    raw_jsonl = (BASE / "nmrshiftdb2_100.jsonl").read_bytes()
    assert digest(raw_jsonl) == manifest["normalized_sha256"]
    assert digest((BASE / "DATABASE-LICENSE.txt").read_bytes()) == manifest["license_sha256"]
    rows = [json.loads(line) for line in raw_jsonl.splitlines()]
    assert len(rows) == manifest["normalized_count"] == 100
    assert len({r["record_id"] for r in rows}) == 100
    counts = Counter()
    for row in rows:
        source = row["source"]
        xml_bytes = (BASE / "selected_sources" / (row["record_id"] + ".xml")).read_bytes()
        sdf_bytes = (BASE / "selected_sources" / (row["record_id"] + ".sdf")).read_bytes()
        assert digest(xml_bytes) == source["selected_xml_sha256"]
        assert digest(sdf_bytes) == row["structure"]["source_sdf_record_sha256"]
        assert Chem.MolFromSmiles(row["structure"]["canonical_smiles"]) is not None
        doc = ET.fromstring(xml_bytes)
        molecule = doc.find("c:molecule", NS)
        assert molecule.get("id") == row["record_id"]
        atom_ids = {a.get("id") for a in molecule.findall("c:atomArray/c:atom", NS)}
        mol_bytes = re.search(rb"<molecule .*?</molecule>\r?\n", xml_bytes, re.S)[0]
        assert digest(mol_bytes) == source["molecule_fragment_sha256"]
        spectra = {s.get("id"): s for s in doc.findall("c:spectrum", NS)}
        assert {s["nucleus"] for s in row["spectra"]} == {"1H", "13C"}
        for s in row["spectra"]:
            counts[s["nucleus"] + "_spectra"] += 1
            e = spectra[s["source_spectrum_id"]]
            assert e.get("moleculeRef") == row["record_id"]
            assert s["solvent"] == row["solvent"]
            assert s["measured_integrations"] == []
            assert s["bibliography_source"]
            assert "nmr:Program" not in s["metadata_source"]
            assert float(s["conditions_source"]["cml:field"]) > 0
            fragment = next(
                b
                for b in re.findall(rb"<spectrum .*?</spectrum>\r?\n", xml_bytes, re.S)
                if ('id="' + s["source_spectrum_id"] + '"').encode() in b
            )
            assert digest(fragment) == s["source_fragment_sha256"]
            source_peaks = e.findall("c:peakList/c:peak", NS)
            assert len(source_peaks) == len(s["assigned_shifts"])
            for original, normalized in zip(source_peaks, s["assigned_shifts"], strict=True):
                assert original.attrib == normalized["source_attributes"]
                assert float(original.get("xValue")) == normalized["shift_ppm"]
                assert original.get("atomRefs", "").split() == normalized["atom_refs"]
                assert set(normalized["atom_refs"]) <= atom_ids
                assert original.get("peakMultiplicity") == normalized["source_multiplicity"]
                assert normalized["measured_integral"] is None
                assert [ET.tostring(c, encoding="unicode") for c in original] == normalized["source_children_xml"]
                counts[s["nucleus"] + "_assigned_peaks"] += 1
        assert row["synthetic_augmentation"] == []
    assert dict(Counter(r["solvent"] for r in rows)) == manifest["solvent_counts"]
    print(json.dumps({"verified_records": len(rows), **counts}, indent=2))


if __name__ == "__main__":
    main()
