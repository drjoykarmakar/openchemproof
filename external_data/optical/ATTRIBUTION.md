# Attribution and redistribution

The raw CSV and the selected/normalized experimental values in this directory derive from:

Joung, Joonyoung F.; Han, Minhi; Jeong, Minseok; Park, Sungnam (2020). **DB for chromophore**, version 2. Figshare. [doi:10.6084/m9.figshare.12045567.v2](https://doi.org/10.6084/m9.figshare.12045567.v2).

Dataset license: [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/). The original license is preserved in `raw/figshare_v2_metadata.json`. The repository's software license does not replace this dataset license.

Associated article: Joung, J. F., Han, M., Jeong, M. & Park, S. **Experimental database of optical properties of organic compounds**. *Scientific Data* **7**, 295 (2020). [doi:10.1038/s41597-020-00634-8](https://doi.org/10.1038/s41597-020-00634-8).

OpenChemProof changes: a deterministic 100-row subset; source and RDKit-canonical SMILES; numeric parsing and explicit nulls; DOI parsing; descriptive solvent names; conversion of reported log10 molar extinction coefficient to linear units; provenance and missingness metadata. `source_rows.csv` preserves all original columns for the selected rows. Individual literature references are retained in every record. No endorsement by the source authors is implied.

Acquisition script is repository software. The data files retain CC BY 4.0 attribution. No PhotochemCAD files or private DyeMind/NarcoticSense tables are included.
