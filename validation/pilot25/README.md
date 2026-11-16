# OpenChemProof: 25-case validation pilot

This package contains **25 source-linked cases, five per category**, with expected verdicts frozen before the first run. The unchanged engine agreed on **24/25 verdicts and required diagnostics**, with no execution errors. All cases still need human scientific review.

This is an **unblinded curator pilot**, not independent experimental validation or a 96% accuracy estimate. Ten cases retain source values without deliberate corruption, eight simulate observations using real source structures/properties, and seven perturb source evidence or identity assertions. The source-only category also includes explicit missing-data adapters; it does not imply complete experimental evidence.

## Use it on your Mac

1. Unzip the package.
2. Copy its `validation` folder into your existing `openchemproof` repository. Merge folders if asked.
3. Open `OpenChemProof-validation-25.xlsx` to review the cases. The amber columns on Results are for the reviewer, review status and notes.
4. In GitHub Desktop, choose **Repository → Open in Terminal**, then run:

```bash
uv run python validation/pilot25/run_pilot.py --output-dir validation/pilot25/results-rerun
```

No API key, model subscription or new experiment is needed. The dependencies must already be installed; if needed, run `uv sync --locked --python 3.12` once. If `uv` is unavailable but the project environment exists, use `.venv/bin/python` in place of `uv run python`.

The runner verifies the frozen file hashes, validates all input schemas, runs the engine and writes a report, CSV, per-case verdicts and content-addressed audit bundles. The original run remains in `results`. An execution/integrity failure returns a nonzero exit code. Scientific disagreement is recorded without aborting; add `--fail-on-disagreement` to make it fail a CI check.

To upload this addition, use GitHub Desktop to commit and push the new `validation` folder. It contains only the small pilot subset; the bulk external-data caches are unnecessary.

## First-run results

| Category | Cases | Expected status and diagnostics matched |
| --- | ---: | ---: |
| LCMS_ADDUCT_DIMER_TRAP | 5 | 5 |
| OPTICAL_IFE_ARTIFACT | 5 | 5 |
| NMR_SOLVENT_INTEGRAL_TRAP | 5 | 5 |
| STEREO_TAUTOMER_TRAP | 5 | 5 |
| VERIFIED_TRUE_AND_NEGATIVE | 5 | 4 |
| Total | 25 | 24 |

The frozen expectations contain 9 ARTIFACT_DETECTED, 5 VERIFIED and 11 INCONCLUSIVE cases. The engine produced 10, 5 and 10, respectively. There are **no qualified NEGATIVE_RESULT controls**, so sensitivity for that class cannot be estimated.

### TR03: unresolved ion metadata

ORD record `ord-7a17fdd932c54be285708fd64e8dcd7b` reports HRMS found **357.1027 m/z** for a product represented as C16H16BClF3NO2. Its retained characterization excerpt does not specify an ion/adduct. It also lacks a measured isotope envelope. The adapter explicitly assumes positive polarity and leaves the adduct missing.

The curator expected INCONCLUSIVE. The engine selects the closest option from its supported positive-adduct list, `[M+H]+`, calculates **358.098747492 m/z**, and emits a critical mass mismatch of approximately **−2781.49 ppm**. It therefore returns ARTIFACT_DETECTED even though the source ion identity is unresolved. This is **not a literal hard-coded proton-adduct default**: the closest supported adduct is chosen.

This disagreement needs adjudication. It does not establish that the source HRMS assignment is correct or that the observed peak is a particular unsupported ion. The next review should consult the original characterization and assess whether missing/unsupported ion identity should produce an inconclusive result. The engine and the frozen expectation were not changed to remove the disagreement.

## How the cases were constructed

- **HRMS:** one genuine processed MS2 surviving-precursor peak, two simulated adduct errors, one simulated chlorine M+2 error and one ambiguity control using two real same-formula constitutional isomers. MassBank precursor metadata were not presented as independently measured HRMS accuracy. The artificial isotope envelope is not an experimental spectrum. Same-formula isomers cannot be distinguished by exact mass or isotope composition alone.
- **Optical:** five Deep4Chem records provide absorption/emission maxima, quantum yields and extinction coefficients. One case retains the missing cuvette evidence. Four use explicit optical simulations. Excitation at the absorption maximum, a 1 cm path length and epsilon_em/epsilon_ex=0.1 are assumptions. The source does not report those experimental settings. Concentration is calculated from A=epsilon*c*l. QY is held as a reference property, not claimed as a measured high-concentration value.
- **NMR:** four ORD records provide reported 1H peak lists with explicit nH assignments. Values are transcribed, not rescaled to the molecular formula. Interval shifts use midpoints and complex multiplicities map to unknown. Two cases add synthetic EtOAc/DCM signals at literature positions; their intensities are stipulated. One removes integral calibration. A fifth case retains NMRShiftDB's missing measured integrals in CD3OD rather than turning atom counts into measured integrals. No raw-spectrum reintegration or qNMR calibration is claimed.
- **Stereo/tautomer:** ChEMBL structures support a matching control, swapped stereoisomer reference, omitted stereochemistry, inconsistent full InChIKey and a generated tautomer comparison. These test identity assertions. They do not validate activity flips, aggregation or experimentally populated tautomers.
- **True/negative:** two source-characterized products test NMR count consistency, one combines NMR with incomplete HRMS ion metadata, and two real reported failures retain missing detection limits and controls. For those failures, `controls_passed=False` means not established, rather than a source-observed failed control. The target product structure describes the intended product and does not assert that it was isolated.

The optical illustration uses `Fcorrected=Fobserved*10^((Aex+Aem)/2)` and the corresponding differential sample/reference factor. Ratios are per-concentration/absorbed-photon response ratios, not raw fluorescence intensities between differently concentrated samples. At sample Aex=1.5 the engine abstains under its declared correction policy. A formula returning 1 in the workbook does not validate that high-OD approximation physically.

## Files and review workflow

- `expected_cases.json`: frozen expectations, rationale, provenance, assumptions and limitations.
- `cases/`: 25 runnable VerificationRecord JSON inputs.
- `sources/`: normalized source snapshots, with original reference IDs and available acquisition hashes.
- `expectations.lock.json`: SHA-256 hashes of inputs, source snapshots and expected labels.
- `results/report.json`: all outcomes, confusion matrix, engine fingerprint and timestamps.
- `results/case_results.csv`: flat comparison table for other analysis tools.
- `results/verdicts/` and `results/audits/`: detailed engine findings and reproducibility bundles.
- `results/remote_comparison.json`: local engine/dependency file comparison with GitHub commit `8e3642c641df37c701bc3887c286b44031ad97c8`; all 16 compared files matched.
- `build_cases.py`: deterministic curation logic. It refuses to overwrite an existing frozen pilot. Rebuilding requires the original repository's external-data corpus.
- `run_pilot.py`: reruns the existing frozen cases without accessing the source repositories or using an LLM.
- `licenses/` and `ATTRIBUTION.md`: source-specific terms and attribution.

For review, read each source snapshot and its cited source, assess the transformation and expected verdict, then fill the workbook's reviewer fields. Preserve this first-run package. If adjudication changes a label or case, create a new version with the reason and reviewer recorded, rather than silently editing the frozen files. Workbook review entries are not automatically imported into the JSON labels. Changing the optical worksheet assumptions recalculates that illustration only; it does not rerun the engine or modify the frozen cases.

VERIFIED has the engine's bounded meaning: no implemented consistency check failed on the supplied evidence. It does not independently establish synthesis success, identity, purity, stereochemical assignment or activity. The small hand-selected sample, shared curator/engine-development context, simulated traps, dependent sources and incomplete modalities prevent a population-level accuracy claim.

The next scientific step is human adjudication of TR03 and the other 24 cases, followed by an independently curated holdout with complete measured ion/isotope data, measured optical settings, integrated spectra and qualified negative controls.
