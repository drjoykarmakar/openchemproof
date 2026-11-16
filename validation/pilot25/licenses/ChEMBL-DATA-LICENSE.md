# ChEMBL data attribution and license

The source-field caches and derived `records.jsonl`, `pairs.jsonl` and manifest are adapted from ChEMBL, https://www.ebi.ac.uk/chembl/, version ChEMBL_37.

ChEMBL data and these adapted data records are distributed under the **Creative Commons Attribution–ShareAlike 3.0 Unported license**:
https://creativecommons.org/licenses/by-sa/3.0/

Official licensing statement:
https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/general-questions

Source citation: Zdrazil et al. *The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods*. Nucleic Acids Research. https://doi.org/10.1093/nar/gkad1004

Changes: subset extraction, exclusion of commercially calculated properties, RDKit structure normalization, same-assay stereoisomer pairing, provenance fields and observed potency-ratio calculation. No new activity measurements were generated.

This data license does not replace the repository's software license for the acquisition script. See `../../docs/sources_structure.md` for selection details, reproducibility instructions and scientific limitations.
