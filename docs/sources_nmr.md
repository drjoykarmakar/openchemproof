# NMRShiftDB2 experimental data and its limits

Contains information from nmrshiftdb2 (www.nmrshiftdb.org), which is made available here under the [nmrshiftdb2 Database License](../external_data/nmr/DATABASE-LICENSE.txt).

## What is included

`external_data/nmr/nmrshiftdb2_100.jsonl` contains paired experimental **assigned chemical shifts** for distinct source molecule IDs. Each record has one 1H and one 13C spectrum in the same explicitly reported solvent, a canonical isomeric SMILES derived with RDKit from the original SDF, original InChI, molecule/spectrum IDs and stable source URLs, literature metadata, source conditions, and exact source-fragment SHA-256 hashes. `manifest.json` is authoritative for counts, selected solvents, retrieval time, source-byte counts and checksums.

The official SDF export was downloaded completely. The XML export was acquired as an explicitly identified byte prefix once enough eligible records were available. The manifest records the exact prefix length/hash and labels it `exact_byte_prefix`; it is **not** a hash or inventory of the full 648,631,898-byte upstream XML export. The unfinished final XML element is discarded. Counts of eligible records apply only to complete elements in that prefix.

`selected_sources/` retains each selected original SDF record and the original XML molecule plus the two selected spectrum elements. XML fragments retain their original bytes inside a new CML wrapper. Individual fragment hashes refer to the exact source element bytes, including line endings; the selected XML file also has its own hash. SDF records can contain other spectra, including predictions; only the two explicitly selected experimental spectra enter the normalized JSONL.

## Selection and reproducibility

Download sources, verified against the official [download documentation](https://nmrshiftdb.nmr.uni-koeln.de/nmrshiftdbhtml/t1.html) and [SourceForge data directory](https://sourceforge.net/projects/nmrshiftdb2/files/data/):

- https://downloads.sourceforge.net/project/nmrshiftdb2/data/nmrshiftdb2withsignals.sd
- https://downloads.sourceforge.net/project/nmrshiftdb2/data/nmrshiftdb2.xml
- License: https://svn.code.sf.net/p/nmrshiftdb2/code/trunk/snapshots/nmrshiftdb2datalicense.txt

Validate the compact subset without downloading anything using `.venv/bin/python external_data/nmr/verify_subset.py`; this compares every retained peak and hash against the original source fragments.

The importer requires an explicit supported solvent, positive reported spectrometer frequency, and bibliography on **both** selected spectra; it excludes computational `Program`, `NMRModel`, or `GeomMethod` metadata. The official [condition definitions](https://svn.code.sf.net/p/nmrshiftdb2/code/trunk/nmrshiftdb2/src/sql/populate-mysql.sql) classify solvent and field conditions as measured (`m`) and Program conditions as calculated (`c`) or predicted (`p`). The [exporter implementation](https://svn.code.sf.net/p/nmrshiftdb2/code/trunk/nmrshiftdb2/src/java/org/openscience/nmrshiftdb/om/DBSpectrum.java) maps these conditions into CML. Thus the classification is grounded in source measurement metadata; a literal experimental boolean is not present in this export.

Selection is deterministic: source order within each solvent, then round robin across CDCl3, DMSO-d6 and CD3OD, taking unique molecule IDs. This convenience subset is not a representative or randomized benchmark of organic chemistry. Pairing asserts a shared structure and solvent, not a shared physical sample or measurement session. Bibliographies are preserved as deposited; this workflow does not independently recheck every original article.

Run from the repository with its RDKit environment:

```sh
.venv/bin/python external_data/nmr/acquire_nmrshiftdb.py --offline --xml-prefix-bytes PREFIX_BYTES
```

For this snapshot `PREFIX_BYTES` is `398454784`, also recorded as `xml_prefix_bytes` in the manifest. Omit `--offline` to acquire missing inputs. The prefix option stops reading the XML download at that exact byte count; the full SDF is still acquired. Omit both options to process full downloads. Because upstream filenames can be replaced, compare source hashes with this manifest to verify snapshot identity. Full downloads and research cache files remain under ignored `external_data/nmr/raw/`; they are not needed to inspect the compact selected subset.

## Assigned shifts are not integrated peak lists

The [official CML convention](https://nmrshiftdb.sourceforge.net/nmrshiftdb-convention.html) specifies atom-assigned peaks and optional multiplicity and coupling data. The source exporter writes an intensity as `peakHeight` if present. A peak height is not a measured integral. No atom count, peak height, implicit hydrogen count, or SDF zero intensity is converted into a measured integral.

The normalized records preserve source peak positions, atom references, multiplicity strings, other peak attributes and child coupling elements. Missing multiplicity remains null. `measured_integrations` stays empty and each `measured_integral` is null. The data do not establish complete unassigned observed spectra, quantitative proton integration, sample purity, contaminant presence/absence, or raw-FID availability. An assignment can be incomplete, overlap with another signal, or cover multiple atoms; repeated atom references are preserved. Carbon multiplicity labels must not be interpreted as observed proton splitting patterns.

The verified subset contains **100 real experimental structure/shift pairs** (68 CDCl3, 14 DMSO-d6, 18 CD3OD), 1,649 assigned proton peaks and 1,871 assigned carbon peaks, but **zero quantitatively integrated proton peak-list cases**. One upstream structure with invalid carbon valence was skipped and replaced. Any evaluation requiring measured proton integrals must report that missing-evidence limitation instead of claiming to validate it with this dataset.

## Impurity augmentation provenance

`impurity_reference_positions.json` records only the selected 1H reference positions for ethyl acetate, dichloromethane and water in the three solvents. Source: Gottlieb, Kotlyar and Nudelman, *J. Org. Chem.* 1997, 62, 7512–7515, [DOI 10.1021/jo971176v](https://pubs.acs.org/doi/10.1021/jo971176v), Table 1; [university-hosted article copy](https://nmr.chem.ucsb.edu/docs/NMR_Impurities.pdf), PDF page 3. The JSON records the copy's hash. The combined PDF itself is a research cache, not redistributed as part of the selected data.

| Solvent | EtOAc CH3CO / OCH2 / CH3, ppm | DCM, ppm | Water, ppm |
|---|---|---|---|
| CDCl3 | 2.05 s / 4.12 q / 1.26 t | 5.30 s | 1.56 s |
| DMSO-d6 | 1.99 s / 4.03 q / 1.17 t | 5.76 s | 3.33 s |
| CD3OD | 2.01 s / 4.09 q / 1.24 t | 5.49 s | 4.87 s |

These values are not evidence that any selected NMRShiftDB2 sample contains those impurities. Injected signals must live in a separate `synthetic_augmentation` list or separately labeled derived cases, with their parent source record, reference, injected peak values and generation rule. Injected integrals and jitter are synthetic. The ethyl acetate 3:2:3 proton stoichiometry is not an experimentally measured integral for a selected sample. Water positions are particularly sensitive to conditions; HDO in DMSO can appear near 3.30 ppm, distinct from the tabulated H2O value. A position match alone is not a definitive contaminant identification.

## License boundary

The data are under the **nmrshiftdb2 Database License**, not Apache-2.0 and not an unmodified ODbL. The full upstream license is included; attribution is embedded in each normalized record and the manifest. It includes database share-alike and notice provisions, and section 4.5 adds an OSI-approved software licensing condition for software that depends on the database for its intended functionality. The project source-code license does not relicense upstream data. Retain the data license and notices with this derivative subset.
