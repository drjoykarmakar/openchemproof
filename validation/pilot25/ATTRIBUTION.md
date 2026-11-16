# Sources and changes

The collection preserves source-specific data terms. The repository's Apache-2.0 software license does not replace the licenses attached to source records or their adaptations. Each case identifies its sources, changes and limitations in `expected_cases.json` and the workbook. This subset is neither endorsed nor experimentally verified by the source providers.

## MassBank

MassBank-data release 2026.11, commit `72980b3af504cfd1e89ac54d7698d3da45d85a58`. Selected records are attributed to Stravs M, Schymanski E, Singer H, Department of Environmental Chemistry, Eawag. Each snapshot preserves its author, copyright and literal `CC BY` license statement, source URL and raw-file SHA-256. No version of the license is inferred where the source does not specify one.

Changes: subset selection, normalized structures and peak fields, one surviving-precursor extraction, simulated mass/adduct/isotope perturbations and an isomer ambiguity control. Experimental MS2 fragment lists remain in the snapshots; they are not recast as MS1 isotope envelopes.

## Deep4Chem

Joung, Joonyoung F.; Han, Minhi; Jeong, Minseok; Park, Sungnam (2020). DB for chromophore, version 2. https://doi.org/10.6084/m9.figshare.12045567.v2. CC BY 4.0: https://creativecommons.org/licenses/by/4.0/.

Associated publication: Experimental database of optical properties of organic compounds. Scientific Data 7, 295 (2020). https://doi.org/10.1038/s41597-020-00634-8.

Individual literature references are preserved per record. Changes: select five source records, convert log10 epsilon to linear units and add explicitly simulated cuvette conditions/response ratios. The source excitation setting, measured OD and epsilon at emission remain unknown. No private DyeMind/NarcoticSense-AI data or PhotochemCAD records are included.

## NMRShiftDB2 and impurity reference

NMRShiftDB2: https://nmrshiftdb.nmr.uni-koeln.de/. The selected CD3OD structure/assignment record retains its record URL, spectrum IDs, source hashes and bibliography. The complete downloaded database license is preserved in `licenses/NMRShiftDB2-DATABASE-LICENSE.txt`.

Changes: subset one normalized experimental assignment record. No measured integral is imputed from atom assignments. The engine receives an explicitly incomplete integration list because its paired-array schema cannot represent a shift with a missing integral.

Synthetic impurity reference positions: Gottlieb HE, Kotlyar V, Nudelman A. NMR Chemical Shifts of Common Laboratory Solvents as Trace Impurities. J Org Chem. 1997;62:7512–7515. https://doi.org/10.1021/jo971176v. The injected EtOAc and DCM amounts/integrals are stipulated simulations, not measurements from that paper or the target ORD samples.

## ChEMBL

ChEMBL 37: https://www.ebi.ac.uk/chembl/. CC BY-SA 3.0: https://creativecommons.org/licenses/by-sa/3.0/. The source-specific attribution statement is retained in `licenses/ChEMBL-DATA-LICENSE.md`. Source documents and individual activity IDs are retained in every snapshot.

Changes: select records, compare stereochemical assertions, use one paired isomer's full InChIKey, omit required stereo in one case, and generate a comparison tautomer using RDKit. No new biological measurements, activity-flip labels or aggregator assignments are created.

## Open Reaction Database

Open Reaction Database: https://open-reaction-database.org/ and https://github.com/open-reaction-database/ord-data. The selected data retain CC BY-SA 4.0; the full license is in `licenses/ORD-LICENSE-CC-BY-SA-4.0.txt`.

Selected records derive from the Borylation dataset from Roche and 750 AstraZeneca ELN dataset. Source commit: `398d7349d4218eccf8782a350a13588c7c7fb88f`; normalized snapshots preserve original reaction IDs, dataset and reaction hashes, publication references, and the pinned download location.

Changes: extract reported 1H peak lists and nH assignments, represent reported proton equivalents with scale 1, use midpoints for shift intervals and unknown for complex multiplicities, introduce two explicit contaminant stress tests, remove calibration in one case, retain unspecified ion/adduct data and retain missing negative-control qualifications. Reported yields are provenance only.

## Optical correction reference

The right-angle center-path approximation and its geometry limitations are discussed in the source manuscript available from the Technical University of Denmark: https://backend.orbit.dtu.dk/ws/portalfiles/portal/298034559/Corrections_of_changes_remove_Fluorescence_Spectra_of_oil_and_asphaltenes_OM_12_Simon_Andersen_July_2022._corrected.pdf. The workbook uses a deterministic algebraic simulation; it does not supply measured correction performance.

## Newly authored material

Python scripts follow the repository's Apache-2.0 license, included in `licenses/SOFTWARE-APACHE-2.0.txt`. Source-derived content retains the respective source terms above. Novel explanations, review labels and table organization do not remove those terms. Human reviewers should preserve attribution when reusing or revising this pilot.
