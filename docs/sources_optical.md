# Experimental optical source subset

`external_data/optical/records.jsonl` contains **100 real literature-compiled experimental records**, representing 100 RDKit-canonical chromophores, 100 original literature DOI references, and 18 solvents. Every selected record contains absorption maximum, emission maximum, photoluminescence quantum yield, and the logarithm of the molar extinction coefficient at the absorption maximum. There are **zero simulated records** in this file. The original studies have not been independently re-extracted or experimentally replicated.

## Official source and license

The source is the publication-pinned **Deep4Chem / DB for chromophore v2**, by Joonyoung F. Joung, Minhi Han, Minseok Jeong, and Sungnam Park:

- Dataset: [10.6084/m9.figshare.12045567.v2](https://doi.org/10.6084/m9.figshare.12045567.v2)
- Article and column definitions: [10.1038/s41597-020-00634-8](https://doi.org/10.1038/s41597-020-00634-8)
- Official metadata: <https://api.figshare.com/v2/articles/12045567/versions/2>
- Official CSV: <https://ndownloader.figshare.com/files/23637518>
- Dataset license: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), verified in official Figshare metadata and preserved locally. See `external_data/optical/ATTRIBUTION.md` for credit and changes.

The downloaded CSV has 20,236 rows and 3,406,886 bytes; its MD5 is `6619ca296ebbf14e8c80dda084b43ac5`. The live metadata description mentions a different total, but this subset's counts come from parsing the version-pinned CSV. This is the 2020 published release, not a claim to mirror the current Deep4Chem service. PhotochemCAD was unnecessary because this official release supplies sufficient complete records under an explicit redistribution license. Private DyeMind/NarcoticSense data were not accessed.

## Reproduction and selection

From the repository root, after installing the project dependencies:

```bash
python external_data/optical/acquire_deep4chem.py
# Re-download official inputs and verify the pinned CSV checksum:
python external_data/optical/acquire_deep4chem.py --refresh
```

The script requires RDKit and otherwise uses the Python standard library. This run used RDKit 2026.09.2. The raw CSV and official license metadata are cached in `external_data/optical/raw/`. Each normalized row retains the original Tag, CSV line number, original reference, source file SHA-256, original SMILES, and canonical SMILES.

The full CSV contains 6,565 rows with all four required optical properties; the other 13,671 are excluded for missing required values. The eligible rows also pass numeric ranges, DOI parsing, RDKit SMILES parsing, and single-component solvent filtering. Eligible rows are ranked by SHA-256 of `openchemproof-optical-v1:<Tag>`. The script takes the first 100 with distinct canonical chromophores and distinct reference DOIs, then writes them in numeric Tag order. Exact selected Tags, checksums, software version, and counts are in `manifest.json`.

This is a convenience subset favoring complete properties and reference diversity; it is not a statistically representative sample. The source compilation may contain extraction errors. Canonical SMILES validation establishes parseability, not correctness of the original structure assignment. Do not use this subset to claim measured artifact-detection accuracy.

## Fields and evidence boundaries

| Field | Meaning / evidence |
|---|---|
| `absorption_max_nm` | Reported first absorption maximum; not the experimental excitation wavelength. |
| `emission_max_nm` | Reported emission maximum. |
| `quantum_yield` | Source photoluminescence quantum yield, dimensionless fraction; no percent conversion or imputation. |
| `log10_epsilon_at_absorption_max` | Source base-10 logarithm of molar extinction coefficient at absorption maximum. |
| `epsilon_at_absorption_max_M_inv_cm_inv` | `10 ** log10_epsilon`, in M⁻¹ cm⁻¹. Derived unit representation of the source measurement, not a new measurement. |
| `excitation_wavelength_nm` | Null in all 100 rows; the CSV does not report excitation wavelength. |
| `epsilon_at_excitation_M_inv_cm_inv`, `epsilon_at_emission_M_inv_cm_inv` | Null; peak epsilon cannot establish epsilon at other wavelengths. |
| `absorbance_at_excitation`, `absorbance_at_emission` | Null; no measured OD supplied. |
| `concentration_M`, `path_length_cm`, `pH`, `temperature_K` | Null; these conditions are absent from the CSV. |
| `absorption_spectrum`, `emission_spectrum` | Null; no spectral traces in this CSV. Optional reported FWHMs are retained separately. |
| `fluorescence_lifetime_ns` | Present for 36 selected rows; missing for 64. |
| `record_kind`, `is_simulated` | `literature_compiled_experimental`, `false`. |

Missing source values are JSON null, never zero. Source rows remain in `source_rows.csv` so numeric normalization can be audited without inferring original precision from JSON floats.

## Use in synthetic artifact challenges

These records can anchor **separately labeled simulated scenarios**. For instance, choose a simulated concentration `c` and path length `l` and evaluate `A(lambda_abs_max) = epsilon_max * c * l`. That absorbance is a Beer–Lambert calculation under an assumed setup, not an observed experimental OD. A simulation may explicitly choose excitation at the absorption maximum, but this must remain a scenario assumption rather than populate the original excitation field.

To evaluate absorption at both excitation and emission wavelengths, supply actual spectral data (and concentration/path length or measured absorbance). A peak epsilon alone cannot establish reabsorption or `A_ex + A_em`; FWHM does not uniquely determine a spectrum. Any reconstructed Gaussian spectrum would itself be a model and must be labeled simulated. Do not report inner-filter effects, quenching, aggregation, or observed apparatus artifacts as experimental labels for these original records.

Keep generated perturbations in a different file with `record_kind: simulated_perturbation`, the original record ID, scenario assumptions, perturbation parameters, and a clear synthetic ground-truth label. Evaluation against such challenges measures performance on constructed scenarios, not validation against measured artifact cases.

## Files

- `external_data/optical/records.jsonl`: 100 normalized experimental records.
- `external_data/optical/source_rows.csv`: the same 100 rows, original source columns and strings.
- `external_data/optical/manifest.json`: source, license, checksums, selection, counts, and missingness.
- `external_data/optical/acquire_deep4chem.py`: reproducible acquisition and normalization.
- `external_data/optical/ATTRIBUTION.md`: source attribution, dataset license, and modification notice.
- `external_data/optical/raw/deep4chem_v2.csv`: official, checksum-verified source CSV.
- `external_data/optical/raw/figshare_v2_metadata.json`: official metadata including CC BY 4.0 license.
