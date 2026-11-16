#!/usr/bin/env python3
"""Acquire auditable paired NMRShiftDB2 shifts; never manufacture integrations.

Run from repo: .venv/bin/python external_data/nmr/acquire_nmrshiftdb.py
The full exports stay in ignored raw/. Selected unmodified source fragments are kept.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import UTC, datetime
from pathlib import Path

BASE = Path(__file__).resolve().parent
DOWNLOAD = "https://downloads.sourceforge.net/project/nmrshiftdb2/data/"
LICENSE_URL = (
    "https://svn.code.sf.net/p/nmrshiftdb2/code/trunk/snapshots/nmrshiftdb2datalicense.txt"
)
LICENSE_NAME = "nmrshiftdb2 Database License"
ATTRIBUTION = (
    "Contains information from nmrshiftdb2 (www.nmrshiftdb.org), "
    "which is made available here under the nmrshiftdb2 Database License."
)
NS = {"c": "http://www.xml-cml.org/schema", "b": "http://bibtexml.sf.net/"}
SOLVENTS = {
    "Chloroform-D1 (CDCl3)": "CDCl3",
    "Dimethylsulphoxide-D6 (DMSO-D6, C2D6SO)": "DMSO-d6",
    "Methanol-D4 (CD3OD)": "CD3OD",
}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def file_sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def fetch(url, path, byte_limit=None):
    if path.exists():
        return
    temporary = path.with_suffix(path.suffix + ".part")
    with urllib.request.urlopen(url, timeout=120) as source, temporary.open("wb") as target:
        remaining = byte_limit
        while remaining is None or remaining > 0:
            block = source.read(
                min(1024 * 1024, remaining) if remaining is not None else 1024 * 1024
            )
            if not block:
                break
            target.write(block)
            if remaining is not None:
                remaining -= len(block)
    temporary.replace(path)


def fragments(path, allow_partial=False):
    """Yield exact original byte slices for top-level molecules and spectra."""
    current = []
    kind = None
    with path.open("rb") as handle:
        for line in handle:
            if line.startswith(b"<molecule ") or line.startswith(b"<spectrum "):
                if kind is not None:
                    raise ValueError("Nested or truncated source element")
                kind = "molecule" if line.startswith(b"<molecule ") else "spectrum"
                current = [line]
            elif kind:
                current.append(line)
                if line.strip() == ("</" + kind + ">").encode():
                    yield kind, b"".join(current)
                    current, kind = [], None
    if kind and not allow_partial:
        raise ValueError("Incomplete source download")


def parse_fragment(raw):
    return ET.fromstring(b'<cml xmlns="http://www.xml-cml.org/schema">' + raw + b"</cml>")[0]


def spectrum(raw):
    e = parse_fragment(raw)
    metadata = {m.get("name"): m.get("content") for m in e.findall("c:metadataList/c:metadata", NS)}
    nucleus = metadata.get("nmr:OBSERVENUCLEUS")
    solvents = e.findall('c:substanceList/c:substance[@dictRef="cml:solvent"]', NS)
    conditions = {
        m.get("dictRef"): (m.text or "").strip() for m in e.findall("c:conditionList/c:scalar", NS)
    }
    bibliography = [
        {
            item.tag.split("}")[-1]: (item.text or "").strip()
            for item in entry.iter()
            if not len(item)
        }
        for entry in e.findall("b:entry", NS)
    ]
    if nucleus not in ("1H", "13C"):
        return None
    if any(key in metadata for key in ("nmr:Program", "nmr:NMRModel", "nmr:GeomMethod")):
        return None
    if len(solvents) != 1 or solvents[0].get("title") not in SOLVENTS:
        return None
    try:
        field_mhz = float(conditions.get("cml:field", "0"))
    except ValueError:
        return None
    if not bibliography or field_mhz <= 0:
        return None
    peaks = []
    for p in e.findall("c:peakList/c:peak", NS):
        if p.get("xUnits") != "units:ppm" or p.get("xValue") is None:
            return None
        peaks.append(
            {
                "source_peak_id": p.get("id"),
                "shift_ppm": float(p.get("xValue")),
                "atom_refs": p.get("atomRefs", "").split(),
                "source_multiplicity": p.get("peakMultiplicity"),
                "source_attributes": dict(p.attrib),
                "source_children_xml": [ET.tostring(c, encoding="unicode") for c in p],
                "measured_integral": None,
            }
        )
    if not peaks:
        return None
    return {
        "source_spectrum_id": e.get("id"),
        "source_spectrum_url": "https://www.nmrshiftdb.org/spectrum/"
        + e.get("id").removeprefix("nmrshiftdb"),
        "source_molecule_ref": e.get("moleculeRef"),
        "nucleus": nucleus,
        "solvent": SOLVENTS[solvents[0].get("title")],
        "solvent_source_text": solvents[0].get("title"),
        "conditions_source": conditions,
        "metadata_source": metadata,
        "bibliography_source": bibliography,
        "assigned_shifts": peaks,
        "measured_integrations": [],
        "observed_unassigned_peak_list": None,
        "source_fragment_sha256": sha(raw),
        "_raw": raw,
    }


def candidates(xml_path, allow_partial=False):
    pools = defaultdict(list)
    stats = Counter()
    molecule_raw, entries = None, []

    def finish():
        if molecule_raw is None:
            return
        mol = parse_fragment(molecule_raw)
        by_solvent = defaultdict(dict)
        for entry in entries:
            if entry["source_molecule_ref"] != mol.get("id"):
                raise ValueError("Spectrum does not belong to preceding source molecule")
            by_solvent[entry["solvent"]].setdefault(entry["nucleus"], entry)
        for solvent, pair in by_solvent.items():
            if {"1H", "13C"} <= pair.keys():
                pools[solvent].append(
                    (mol.get("id"), mol.get("title"), molecule_raw, pair["1H"], pair["13C"])
                )

    for kind, raw in fragments(xml_path, allow_partial=allow_partial):
        stats[kind + "_elements_scanned"] += 1
        if kind == "molecule":
            finish()
            molecule_raw, entries = raw, []
        else:
            item = spectrum(raw)
            if item:
                entries.append(item)
                stats["eligible_" + item["nucleus"] + "_spectra"] += 1
    finish()
    return pools, stats


def structures(sdf_path, wanted):
    from rdkit import Chem

    result = {}
    record = []
    with sdf_path.open("rb") as handle:
        for line in handle:
            record.append(line)
            if line.strip() != b"$$$$":
                continue
            raw = b"".join(record)
            record = []
            match = re.search(rb"> <nmrshiftdb2 ID>\s*(\d+)", raw)
            if not match or "nmrshiftdb" + match[1].decode() not in wanted:
                continue
            key = "nmrshiftdb" + match[1].decode()
            text = raw.decode("utf-8").replace("\r\n", "\n")
            mol = Chem.MolFromMolBlock(text.split("M  END")[0] + "M  END\n")
            if mol is None:
                continue
            inchi = re.search(r"> <INChI>\s*(\S+)", text)
            result[key] = {
                "canonical_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
                "inchi": inchi[1] if inchi else None,
                "source_sdf_record_sha256": sha(raw),
                "_raw": raw,
            }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=100)
    parser.add_argument("--offline", action="store_true")
    parser.add_argument(
        "--xml-prefix-bytes",
        type=int,
        default=None,
        help=("Acquire/use this exact XML prefix; discard the incomplete final element"),
    )
    args = parser.parse_args()
    if args.limit < 1:
        parser.error("--limit must be positive")
    raw_dir = BASE / "raw"
    raw_dir.mkdir(exist_ok=True)
    paths = [raw_dir / "nmrshiftdb2.xml", raw_dir / "nmrshiftdb2withsignals.sd"]
    if not args.offline:
        for path in paths:
            fetch(
                DOWNLOAD + path.name, path, args.xml_prefix_bytes if path.suffix == ".xml" else None
            )
        fetch(LICENSE_URL, BASE / "DATABASE-LICENSE.txt")
    with paths[0].open("rb") as handle:
        handle.seek(-100, 2)
        xml_complete = handle.read().strip().endswith(b"</cml>")
    if not xml_complete and args.xml_prefix_bytes != paths[0].stat().st_size:
        raise ValueError(
            "Incomplete XML requires explicit --xml-prefix-bytes equal to its byte size"
        )
    pools, stats = candidates(paths[0], allow_partial=not xml_complete)
    eligible_counts = {key: len(value) for key, value in pools.items()}
    selected, seen = [], set()
    while any(pools.values()):
        for solvent in ("CDCl3", "DMSO-d6", "CD3OD"):
            while pools[solvent] and pools[solvent][0][0] in seen:
                pools[solvent].pop(0)
            if pools[solvent]:
                item = pools[solvent].pop(0)
                selected.append(item)
                seen.add(item[0])
    structs = structures(paths[1], seen)
    out = BASE / "selected_sources"
    out.mkdir(exist_ok=True)
    normalized = []
    for molid, title, mol_raw, proton, carbon in selected:
        if len(normalized) >= args.limit:
            break
        if molid not in structs:
            stats["selected_structure_parse_failures"] += 1
            continue
        structure = structs[molid]
        sdf_raw = structure.pop("_raw")
        (out / (molid + ".sdf")).write_bytes(sdf_raw)
        source_xml = (
            b'<?xml version="1.0" encoding="UTF-8"?>\n'
            b'<cml xmlns="http://www.xml-cml.org/schema">\n'
            + mol_raw
            + proton.pop("_raw")
            + carbon.pop("_raw")
            + b"</cml>\n"
        )
        (out / (molid + ".xml")).write_bytes(source_xml)
        record = {
            "record_id": molid,
            "name": title,
            "structure": structure,
            "source": {
                "database": "nmrshiftdb2",
                "molecule_id": molid.removeprefix("nmrshiftdb"),
                "record_url": "https://www.nmrshiftdb.org/molecule/"
                + molid.removeprefix("nmrshiftdb"),
                "download_url": DOWNLOAD + "nmrshiftdb2.xml",
                "license": LICENSE_NAME,
                "license_url": LICENSE_URL,
                "attribution": ATTRIBUTION,
                "molecule_fragment_sha256": sha(mol_raw),
                "selected_xml_sha256": sha(source_xml),
            },
            "evidence_kind": "experimental_assigned_shifts",
            "experimental_selection_basis": (
                "experimental metadata branch: no computational Program/GeomMethod/NMRModel; "
                "positive instrument frequency; literature reference; supported explicit solvent"
            ),
            "solvent": proton["solvent"],
            "spectra": [proton, carbon],
            "synthetic_augmentation": [],
            "limitations": [
                "Assigned shifts are not quantitative integrated peak lists.",
                "Measured integrals are absent; atom-reference counts are not measured integrals.",
                "Pairing is by structure and solvent, not a claim of the same physical sample.",
                "Source multiplicities are retained verbatim, never inferred.",
            ],
        }
        normalized.append(record)
    output = BASE / "nmrshiftdb2_100.jsonl"
    output.write_text(
        "".join(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n" for row in normalized)
    )
    manifest = {
        "schema_version": 1,
        "dataset": "OpenChemProof NMRShiftDB2 paired experimental assigned shifts",
        "retrieved_utc": datetime.now(UTC).isoformat(),
        "requested_count": args.limit,
        "normalized_count": len(normalized),
        "quantitatively_integrated_case_count": 0,
        "solvent_counts": dict(Counter(r["solvent"] for r in normalized)),
        "eligible_same_solvent_pairs_before_structure_validation": eligible_counts,
        "scan_counts": dict(stats),
        "selection": (
            "deterministic source-order round robin across CDCl3, DMSO-d6, CD3OD; "
            "unique molecule IDs; skip unparseable structures"
        ),
        "source_exports": [
            {
                "url": DOWNLOAD + p.name,
                "filename": p.name,
                "bytes": p.stat().st_size,
                "sha256": file_sha(p),
                "complete": xml_complete if p.suffix == ".xml" else True,
                "snapshot_kind": ("full_export" if xml_complete else "exact_byte_prefix")
                if p.suffix == ".xml"
                else "full_export",
            }
            for p in paths
        ],
        "xml_export_complete": xml_complete,
        "xml_prefix_bytes": None if xml_complete else paths[0].stat().st_size,
        "scope": (
            "Counts apply only to scanned complete elements; "
            "partial-export counts are not whole-database totals."
        ),
        "normalized_sha256": file_sha(output),
        "license": LICENSE_NAME,
        "license_url": LICENSE_URL,
        "license_sha256": file_sha(BASE / "DATABASE-LICENSE.txt"),
        "attribution": ATTRIBUTION,
        "blocker_to_integrated_peak_list_target": (
            "Selected source CML has assigned chemical shifts but no measured proton integration "
            "values. No integral or impurity measurement is fabricated."
        ),
    }
    (BASE / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "normalized_count": len(normalized),
                "solvent_counts": manifest["solvent_counts"],
                "quantitatively_integrated_case_count": 0,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
