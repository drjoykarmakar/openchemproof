# Evidence collection and benchmark dataset card

## Two different collections

`benchmarks/chem_artifact_bench/data/chem_artifact_bench_v 1.*` contains 500 **synthetic integration
fixtures**. `external_data/` contains 500 **real public-source records**. They are separate datasets,
not interchangeable versions of the same ground truth. The public collection has no engine-derived
gold verdicts. A uniform index records native row pointers, category, source license, source URL,
row hash, available evidence and missing evidence; the original normalized measurements remain
in their source-specific files.

The initial public collection was acquired on 2026-11-16. It is a curated feasibility corpus,
not a randomized representative sample of all chemistry. Source records may contain errors,
uncertain annotations or incomplete conditions. Each source document describes acquisition,
selection and changes from upstream values.

| Collection | Unit of record | Count | Independent comparison units |
|---|---|---:|---|
| MassBank | Unique compound with one selected processed MS2 spectrum |100|100 compounds; one selected contributor/instrument family |
| Deep4Chem | Unique chromophore/solvent literature row |100|100 distinct underlying literature DOIs |
| NMRShiftDB2 | Structure with same-solvent proton/carbon assignment pair |100|100 structures; physical-sample equivalence is unconfirmed |
| ChEMBL | Compound/activity row from same-assay stereoisomer pair |100|50 pairs,21 assays,19 documents |
| ORD | Reaction/product outcome |100|60 successes across 9 publications plus 40 failures from another dataset |

## Measurement lineage

Preserve four separate concepts: an original observed value; a source assertion; a deterministic
derived value; a simulated or adversarially modified value. Examples:

- MassBank's processed MS2 peaks are experimental. Its precursor metadata may be calculated;
  it cannot serve as independent observed mass accuracy. The selected set contains 78 positive
  `[M+H]+` and 22 negative `[M-H]-` records, with no measured MS1 isotope envelopes.
- Deep4Chem's source log ε is measured/compiled at the absorption maximum. Its conversion to
  linear ε is derived. Absorption maximum is not a recorded excitation setting. A simulated
  Beer–Lambert absorbance needs an explicitly chosen path length/concentration; absorption at
  emission remains unknown unless measured or explicitly modeled.
- NMRShiftDB2 atom assignments and peak multiplicity annotations are preserved. Counting atom
  references is not quantitative integration; no measured integrals were fabricated.
- ChEMBL potency ratios are descriptive comparisons. Same assay/endpoint/units/document is a
  useful matching criterion but does not establish statistically significant activity inversion.
- ORD product-linked characterization supports a source claim. Zero yield plus a failure narrative
  is retained as an outcome, not automatically converted to a strict negative with successful
  analytical controls and a known detection limit.

## Constructing future challenges

Each derived challenge should keep the unchanged source row, immutable source ID, transformation
version, random seed, exact changed fields, injected values, assumptions and a new derivative ID.
Store expected behavior separately from model inputs. The ground truth for an injected error is
the recorded transformation; it is not a new experimental artifact observation.

| Proposed transformation | Defensible interpretation | Requirement before labeling a real artifact |
|---|---|---|
| Relabel an adduct | Corrupted source annotation challenge | Independent measured precursor and correct ion/adduct evidence |
| Substitute an equal-formula isomer | Ambiguity/abstention challenge if no orthogonal identity data | Chiral/structural evidence able to distinguish the isomers |
| Increase OD using ε and chosen concentration | Simulated IFE stress test | Measured sample/reference absorbance, geometry and fluorescence response |
| Add solvent NMR peaks | Synthetic contamination challenge | Real integration/calibration or independently confirmed contaminant identity |
| Swap a stereo identity assertion | Deterministic structural-identity contradiction | Independent experimental stereo assignment for a real misidentification |
| Treat a zero yield as a negative | Source-outcome annotation task | Detection limit, positive/negative controls and target-specific assay validation |

Do not fill missing MS1 envelopes with calculated isotope ratios and call them experimental.
Do not turn all PAINS/aggregator-risk molecules into negative labels. Keep valid controls with
these alerts, and preserve true tautomers and isobaric ambiguities to test specificity.

## Evaluation and leakage

The synthetic oracle currently checks 500/500 construction labels and their expected flags. This
is a regression check against shared chemistry primitives, not validation of general diagnostic
accuracy. Parameter jitter and repeated structures are not 500 independent chemical systems.
Live baseline/MCP comparison is implemented and tested with mocked provider responses; no paid
model evaluation has been run or reported.

Before scientific evaluation, curate an external holdout independently of rule development.
Group split by molecule/scaffold, stereoisomer pair, source document and experiment as appropriate.
Keep all perturbations and controls from one parent record in the same split. The ORD positive
and negative subsets differ in chemistry family/source, providing an obvious shortcut for a
classifier. Report source-stratified performance and uncertainty, not just a pooled accuracy.
No model-training split is claimed for this first acquisition.

## Licensing and redistribution

Data remain separate subcollections with their native licenses. The software Apache license
does not override CC-BY-SA or NMRShiftDB2 terms. Exact attribution, transformations and notices
are retained. See `NOTICE`, each source directory, and `docs/sources_*.md`. Bulk NMR exports are
excluded; selected original fragments, checksums and reproducible acquisition scripts are included.
