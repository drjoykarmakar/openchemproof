# Real ChEMBL stereoisomer evidence

`external_data/structure/records.jsonl` contains **100 distinct ChEMBL compounds with measured activity**, arranged as **50 disjoint stereoisomer pairs** in `pairs.jsonl`. These are real retrieved records, not generated stereoisomers or synthetic potency values. The subset uses 21 assays against the human D(2) dopamine receptor (`CHEMBL217`). It is an evidence/identity fixture, not a gold standard for predicting a binary activity flip.

## Source and redistribution

- Source: [ChEMBL at EMBL-EBI](https://www.ebi.ac.uk/chembl/), public [REST API](https://www.ebi.ac.uk/chembl/api/data/docs).
- Retrieved release: **ChEMBL_37**, release date **2026-05-01**, as returned by the API status endpoint; the exact response and retrieval timestamp are retained in `api_status.json`.
- ChEMBL data and these adapted data records: **Creative Commons Attribution–ShareAlike 3.0 Unported**, [license](https://creativecommons.org/licenses/by-sa/3.0/). See the [official licensing FAQ](https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/general-questions) and [ChEMBL licensing notice](https://chembl.github.io/chembl-licensing/). The data license is distinct from the repository's software license. Preserve attribution and share alike when redistributing this data adaptation.
- Attribution: ChEMBL data is from https://www.ebi.ac.uk/chembl/; version ChEMBL_37. Cite Zdrazil et al., *The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods*, Nucleic Acids Research, [DOI 10.1093/nar/gkad1004](https://doi.org/10.1093/nar/gkad1004). API reference: Davies et al., [DOI 10.1093/nar/gkv352](https://doi.org/10.1093/nar/gkv352).
- Commercially calculated molecule properties and ligand-efficiency fields are not extracted. The retained values are structures, activity measurements, and their assay/document provenance.

## Acquisition and selection

Run from the repository root with RDKit installed:

```bash
python external_data/structure/acquire_chembl.py --max-pages 10 --pairs 50
```

The script requests the first 10,000 D2 activity records in ascending `activity_id` order with `standard_type` in `Ki,IC50`, `standard_relation` equal to `=`, and `standard_units` equal to `nM`. The raw cache contains the selected source fields, exact request URLs, timestamps and API pagination metadata. It excludes records carrying a ChEMBL data-validity warning or potential-duplicate flag and requires a positive numeric value.

For each molecule, RDKit parses the source canonical SMILES and requires at least one tetrahedral stereocenter with no unassigned tetrahedral centers. It generates a canonical isomeric SMILES, a canonical SMILES after removing stereochemistry, and a locally calculated InChIKey. Pairs must have identical stereo-stripped connectivity but distinct canonical isomeric SMILES. No stereoisomer is enumerated or invented. This definition includes enantiomeric or diastereomeric relationships and does not assert which class applies to each pair.

Activity comparisons additionally require identical ChEMBL assay, document, target, endpoint and standardized units, plus identical assay description, recorded assay variants, activity properties and target organism. Each ChEMBL molecule occurs in only one selected pair. The script selects the first 50 qualifying disjoint pairs in deterministic source activity order; it does not select on the size of the activity difference. For repeated measurements of one isomer in a group, the earliest activity ID is retained; values are not averaged and disagreements are not resolved.

The output preserves original and standardized value/relation/type/units, molecule and activity IDs, target identity, full available assay metadata and description, activity properties, source document metadata (including DOI/PMID where reported), and resolvable API URLs. Unreported conditions remain unknown; temperature, pH, duration, substrate or radioligand concentrations, and uncertainty are not guessed from common assay protocols.

For an offline deterministic rebuild from the included source cache:

```bash
python external_data/structure/acquire_chembl.py --offline --pairs 50
```

`manifest.json` records the source release, exact counts, RDKit version, selection rule and SHA-256 checksums. The live API changes releases; cache-backed reproduction is the version-preserving route. A new live acquisition needs a fresh cache directory if a different release is desired.

## What the evidence supports

A pair supports that ChEMBL records two distinct specified stereochemical structures with matched-assay numeric measurements. The `observed_potency_ratio` is simply the larger reported nM value divided by the smaller. Lower Ki/IC50 is ordinarily a more potent measurement within that same endpoint; the ratio is not a significance test. Every `activity_flip_label` is `null`: there is no user-specified active/inactive threshold, and the source generally does not provide enough replicate uncertainty to establish a statistically meaningful flip.

Identical assay IDs strengthen comparability but cannot guarantee that all experimental details were identical or completely reported. Source curation and structure assignment errors remain possible. This is a convenience sample from one target and literature-biased early activity IDs; it is not balanced, chemically independent of ChEMBL training corpora, or representative of all stereochemistry.

**Coverage limits:** 100 molecule/activity records, 50 paired comparison cases, 0 tautomer cases, and 0 experimentally confirmed aggregation-artifact cases. Aggregation and tautomer status are explicitly `not_assessed`. No PAINS match, aggregation-risk score, or structural similarity is treated as experimental confirmation of aggregation. Do not count the 100 rows as 100 independent stereoisomer comparisons.
