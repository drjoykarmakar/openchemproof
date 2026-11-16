# Public experimental LC–MS/MS source records

`external_data/lcms/massbank_100.jsonl` is a normalized acquisition subset of real
MassBank experimental reference records. `manifest.json` supplies exact counts,
source version, data hashes, selection rules and limitations. The unmodified source
bytes are in `external_data/lcms/raw/`, and each JSONL row records its original
file SHA-256, immutable GitHub URL, author attribution and source license.

## Source and reproducibility

- Official repository: <https://github.com/MassBank/MassBank-data>.
- Release: **2026.11**, archived at <https://doi.org/10.5281/zenodo.24082491>.
- Pinned commit: `72980b3af504cfd1e89ac54d7698d3da45d85a58`.
- Contributor: Eawag; EA series with the source comment `CONFIDENCE standard compound`.
- Record-format reference: <https://github.com/MassBank/MassBank-documentation/blob/dev/MassBankRecordFormat.md>.

Acquisition is restricted to source records with an `LC-` instrument type, MS2
acquisition, valid RDKit-parsable SMILES, a positive source-reported precursor m/z,
ion mode, adduct, and a nonempty peak list matching `PK$NUM_PEAK`. Selection uses
ascending source compound IDs and first available spectra per polarity block,
preferring alternating polarity when available. Canonical isomeric SMILES are
de-duplicated. This convenience subset does not represent all of MassBank, and it
is neither a random sample nor a blind challenge.

From the repository root, using an environment with RDKit installed:

```bash
python external_data/lcms/acquire_massbank.py --count 100
```

Downloads use a pinned commit and a local ignored cache. `--offline` regenerates
the JSONL from that cache when present. Source bytes are never rewritten. SMILES
normalization uses RDKit canonical isomeric SMILES without salt stripping,
neutralization or tautomer canonicalization; the RDKit version and original
SMILES are recorded. Canonicalization does not establish chemical identity.

## What is experimental, reported, inferred, or missing

| Content | Status | Permitted interpretation |
| --- | --- | --- |
| `PK$PEAK` m/z and intensities | Published processed experimental MS2 spectrum | Fragment-spectrum measurements, with source processing retained |
| SMILES, formula and compound identity | Source reference assignment | Standard-compound attribution, not independent re-identification |
| `MS$FOCUSED_ION: PRECURSOR_M/Z` | Source-reported precursor metadata | May be calculated; not independent measured mass accuracy |
| Adduct, polarity and instrument | Source acquisition metadata | Preserved verbatim where applicable and normalized separately |
| Canonical SMILES | Computed representation | Deterministic parsing/serialization only |
| M, M+1, M+2 precursor isotope envelope | **Unavailable in this subset** | All three are null; no isotope-agreement claim is supported |
| Raw instrument files or replicates | Not supplied | Cannot estimate raw-data quality, uncertainty, or repeatability |
| Synthetic perturbations | None | Must be kept in a separate dataset if generated later |
| Expected verifier verdict | Unassigned | No model output or prediction-derived oracle labels |

The MassBank format explicitly prefers **calculated exact mass** for its
`PRECURSOR_M/Z` field (section 2.5.1). Treating agreement between that field and a
formula mass calculation as experimental validation would be circular. Use this
subset for ingestion, provenance and metadata-consistency checks; it does not
establish experimental LCMS verifier accuracy. An independently annotated measured
precursor peak or an MS1 scan is needed for measured mass-error evaluation.

MS2 fragment peaks, including a surviving precursor-like peak, are not an MS1
isotope envelope. No closest-peak heuristic, formula-derived isotope distribution
or synthetic M+1/M+2 intensity is inserted into the missing observations. The
source processing fields are retained because these reference peaks may have been
centroided, recalibrated or filtered by RMassBank. Automated fragment formulas,
where present in the original records, are tentative source annotations rather
than independently verified labels.

## Attribution and licensing

Every selected source record states **CC BY**. The original records do not state
a Creative Commons license version, so this subset preserves the literal license
and leaves its version unspecified. Source data and their normalized
representations are not relicensed under OpenChemProof's software license. Preserve
the source authors, copyright, record URLs, original license notice and the fact
of normalization when redistributing this subset. Attribution is included per row
and in every unmodified raw record.

MassBank's record-level licensing approach is documented in the record format,
section 2.1.5, and confirmed by a maintainer at
<https://github.com/MassBank/MassBank-data/issues/394>.

## Evaluation traps to keep separate

- Missing isotope data must produce an explicit unperformed isotope check, not a
  fabricated match or a negative finding.
- Source precursor-metadata agreement cannot be counted as independent measured
  PASS accuracy.
- Multiple collision energies for one structure are correlated observations;
  split by structure/source before any learned evaluation.
- Polarity or adduct edits are synthetic tests, even when their parent row came
  from a real record. Keep the original accession and the exact perturbation.
- A source-standard record is not a negative ground-truth example. Real
  contradictions require documented independent adjudication.
- `LC-ESI-*` records must not be mixed silently with direct-infusion ESI records.

No MoNA records were needed for this subset. There are no invented accessions,
generated spectra or inferred isotope intensities in the acquired cohort.
