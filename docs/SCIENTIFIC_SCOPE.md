# Scientific scope and policy

OpenChemProof tests internal consistency of supplied observations and assertions. It cannot
authenticate an instrument export, establish identity from formula alone, or certify that a
record accurately represents the laboratory experiment. File hashes detect changes relative to
the saved snapshot; they do not prove the original evidence is true.

## Structural standardization

RDKit SaltRemover removes recognized fragments and the software records actual removed SMILES,
including multiplicity. Remaining unresolved components are not silently reduced to the largest
fragment. Neutralizable salts are uncharged to the parent; permanent ions are outside the mass
model. `mass_basis` states whether a reported formula/mass refers to the neutral parent or the
as-drawn salt. `canonical_smiles`, `inchikey` and related fields are assertions to check, not trusted
computed values.

The first 14 InChIKey characters encode the connectivity hash. A different remainder can reflect
stereo, isotope or protonation differences; the software does not call every full-key difference
a stereochemical inversion. R/S and E/Z identity comparisons use isomeric graph representations.
Canonical tautomer equivalence does not establish the tautomer distribution under assay conditions.
Enumeration is capped at 256 tautomers/512 transforms; reaching a cap yields an inconclusive flag.
PAINS/Brenk and `cLogP >4.5` plus ≥3 aromatic rings are risk screens, not experimental conclusions.

## HRMS

Eight neutral-parent singly charged ESI adducts are supported, including ammonium, potassium,
chloride, formate and protonated dimers. Adduct mass increments are the requested rounded constants;
they are included in the MCP reference resource. Ppm is signed `(observed-theoretical)/theoretical*1e6`.
The 5 ppm acceptance and 3 ppm alternative gates are configurable policies.

The isotope model convolves natural-abundance polynomials, retaining all terms up to nominal M+2
including cross terms. Adduct atoms and dimers change the composition. Explicit isotope labels
remain fixed. Supported elements: C/H/N/O/S/Cl/Br/F/P/I/Na/K. The returned ratios are relative to
M=1; M, M+1, M+2 are not normalized to sum=1 and high-order isotope peaks are not renormalized in.
Natural abundance is a representative assumption, not a sample-specific measurement. Other
elements/enrichment, resolved fine structure, mixtures and multiply charged ions need a richer model.

Isotope residual tolerances are max(0.05 absolute M-relative ratio, 25% of expected ratio).
They are screening defaults, not empirically calibrated instrument-error estimates. A sodium or
dimer match is a candidate explanation. If both target and starting material fit, the result is
inconclusive. Equal-formula isomers can have identical exact masses and nominal isotope patterns.
Missing isotope evidence remains missing: MS2 fragment peaks are not reused as M+1/M+2 peaks.

## Fluorescence

The center-observation approximation multiplies fluorescence by
`10**((A_ex+A_em)/2)`. It assumes a homogeneous, non-scattering, right-angle geometry with
effective half-path lengths. Relative response requires correction of both sample and reference.
An apparent change is flagged as IFE-explained only when correction returns the ratio within 15%
of unity, the uncorrected change exceeds 15%, and at least one absorbance exceeds 0.15.
The default refuses application above OD 1.0 or outside the stated geometry. Even inside that
range, photobleaching, reabsorption, scattering and instrumentation can invalidate the approximation.
Dilution series and blanks remain essential experimental controls.

Ordinary single-photon fluorescence QY is checked against 0–1. The model does not address quantum
cutting or multiphoton yield conventions. Anti-Stokes emission without a declared alternative
mechanism conflicts with the ordinary model; it is not universally forbidden by physics. The 15 nm
excitation exclusion window is a conservative screening rule, not a complete Rayleigh/Raman model.
An explicitly supplied solvent Raman shift can trigger a tentative overlap flag. Solvatochromism
checks use a declared trend and ≥3 distinct polarities; specific solvation remains an alternative.

## NMR

The H budget is determined from the normalized RDKit structure, excluding deuterium/tritium.
H bound to O/N/S is considered exchangeable. D2O/CD3OD suppresses the expected exchangeable
contribution in the model. CDCl3 and DMSO-d6 allow a range unless exchange visibility is explicitly
specified; DMSO does not guarantee sharp OH/NH peaks for all samples. Carbon-bound H/D exchange,
kinetics and unusual nuclei require additional modeling.

Integrals must carry an independent `integral_scale` (reported units per proton). Scaling the
total to the expected formula would hide contamination and is never done. Default H tolerance
is max(0.35 H, 5% of the upper expected budget). Solvent-matching shifts flag candidates without
automatic subtraction. The analyte itself may be ethyl acetate or overlap those shifts. Carbon
counts are an upper-bound consistency aid; symmetry/overlap reduce peaks and rotamers can increase
them, so extra carbon signals produce an inconclusive warning rather than definitive failure.

## Verdict and audit semantics

Precedence: a critical inconsistency → ARTIFACT_DETECTED; otherwise a sufficiently documented
negative with successful controls → NEGATIVE_RESULT; incomplete/unsupported evidence or structure
metadata alone → INCONCLUSIVE; otherwise → VERIFIED. “Artifact” includes inconsistent identity
assertions; it is not always an experimentally confirmed instrumental artifact. The
`confidence_score` is a deterministic rule-completion indicator, not a calibrated probability.
Warnings can coexist with VERIFIED; the flags must be reviewed.

Negative evidence references are retained, not dereferenced/authenticated automatically. Ordinarily,
external evidence review would be required before using a result for regulatory or publication claims.
ORD ID fields support identifier linkage only: this is not a full ORD schema or bidirectional adapter.

## Sources

- [RDKit SaltRemover](https://www.rdkit.org/docs/source/rdkit.Chem.SaltRemover.html),
  [descriptors](https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html),
  [MolStandardize](https://www.rdkit.org/docs/source/rdkit.Chem.MolStandardize.rdMolStandardize.html).
- Gottlieb, Kotlyar and Nudelman. NMR Chemical Shifts of Common Laboratory Solvents as Trace
  Impurities. J. Org. Chem. 1997, 62, 7512–7515. [DOI 10.1021/jo971176v](https://doi.org/10.1021/jo971176v).
- On the origin and correction for inner filter effects in fluorescence, Parts I and II:
  [PubMed 32428893](https://pubmed.ncbi.nlm.nih.gov/32428893/) and
  [PubMed 34032610](https://pubmed.ncbi.nlm.nih.gov/34032610/).
- [MassBank record format](https://github.com/MassBank/MassBank-documentation/blob/dev/MassBankRecordFormat.md):
  precursor metadata may be calculated, which prevents independent mass-accuracy validation.
- [FastMCP client and testing](https://gofastmcp.com/docs/clients/client).
- [OpenAI function calling](https://developers.openai.com/api/docs/guides/function-calling).

Quantitative limits not attributed above are explicit v0.1 engineering policies. Publication-level
validation requires independently curated, instrument-aware measured datasets and blinded evaluation.
