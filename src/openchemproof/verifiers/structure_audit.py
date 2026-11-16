"""Structural consistency, identity-layer comparisons, and non-diagnostic risk alerts."""

from functools import lru_cache
from typing import Any, cast

from rdkit import Chem
from rdkit.Chem import Crippen, FilterCatalog, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

from openchemproof.schema.models import CompoundSpec, ModalityResult
from openchemproof.verifiers.common import flag, prepare


@lru_cache(maxsize=1)
def _catalog() -> FilterCatalog.FilterCatalog:
    parameters = FilterCatalog.FilterCatalogParams()
    parameters.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    parameters.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog.FilterCatalog(parameters)


def audit_structure(compound: CompoundSpec) -> ModalityResult:
    flags = []
    metrics: dict[str, float | str] = {}
    checks = ["smiles_parse"]
    try:
        prepared = prepare(compound.proposed_smiles)
    except ValueError as exc:
        return ModalityResult(
            modality="structure",
            flags=[flag("INVALID_SMILES", "CRITICAL_ARTIFACT", "structure", str(exc))],
            checks_performed=checks,
        )
    parent = prepared.parent
    canonical = str(Chem.MolToSmiles(parent))
    key = str(cast(Any, Chem).MolToInchiKey(parent))
    mass = float(rdMolDescriptors.CalcExactMolWt(parent))
    basis = prepared.original if compound.mass_basis == "as_drawn" else parent
    basis_mass = float(rdMolDescriptors.CalcExactMolWt(basis))
    basis_formula = str(rdMolDescriptors.CalcMolFormula(basis))
    metrics.update(
        canonical_smiles=canonical,
        inchikey=key,
        inchi_connectivity_layer=key[:14],
        molecular_formula=str(rdMolDescriptors.CalcMolFormula(parent)),
        parent_monoisotopic_mass=mass,
        as_drawn_monoisotopic_mass=float(rdMolDescriptors.CalcExactMolWt(prepared.original)),
    )
    inconclusive = prepared.unresolved_fragments
    if prepared.unresolved_fragments:
        flags.append(
            flag(
                "UNRESOLVED_MIXTURE",
                "WARNING",
                "structure",
                "Multiple unrecognized components remain; parent identity is ambiguous",
            )
        )
    if prepared.stripped_fragments:
        flags.append(
            flag(
                "SALT_STRIPPED",
                "INFO",
                "structure",
                "Counterions removed; neutral parent and as-drawn masses retained",
                removed_fragments=".".join(prepared.stripped_fragments),
            )
        )
    if Chem.GetFormalCharge(parent) != 0:
        inconclusive = True
        flags.append(
            flag(
                "PERMANENT_ION",
                "WARNING",
                "structure",
                "Permanent ions require an explicit ion model beyond neutral-parent v0.1",
            )
        )
    if compound.molecular_formula is not None and compound.molecular_formula != basis_formula:
        flags.append(
            flag(
                "FORMULA_MISMATCH",
                "CRITICAL_ARTIFACT",
                "structure",
                "Reported formula differs from the declared mass basis",
                reported=compound.molecular_formula,
                calculated=basis_formula,
            )
        )
    if compound.monoisotopic_mass is not None and abs(compound.monoisotopic_mass - basis_mass) > 0.001:
        code = "SALT_MASS_MISMATCH" if prepared.stripped_fragments else "EXACT_MASS_MISMATCH"
        flags.append(
            flag(
                code,
                "CRITICAL_ARTIFACT",
                "structure",
                "Reported neutral mass disagrees with the declared molecular basis",
                reported=compound.monoisotopic_mass,
                calculated=basis_mass,
            )
        )
    if compound.inchikey is not None and compound.inchikey != key:
        code = "INCHI_CONNECTIVITY_MISMATCH" if compound.inchikey[:14] != key[:14] else "INCHI_LAYER_MISMATCH"
        flags.append(
            flag(
                code,
                "CRITICAL_ARTIFACT",
                "structure",
                "InChIKey assertion conflicts with calculated parent identity; full-key "
                "differences can include stereo, isotope, or protonation layers",
                reported=compound.inchikey,
                calculated=key,
            )
        )
    if compound.inchi_connectivity_layer is not None and compound.inchi_connectivity_layer != key[:14]:
        flags.append(
            flag(
                "CONNECTIVITY_LAYER_MISMATCH",
                "CRITICAL_ARTIFACT",
                "structure",
                "Reported connectivity block disagrees with calculated InChIKey",
            )
        )
    enumerator = rdMolStandardize.TautomerEnumerator()
    enumerator.SetMaxTautomers(256)
    enumerator.SetMaxTransforms(512)
    tautomers = enumerator.Enumerate(parent)
    metrics["tautomer_count"] = float(len(tautomers))
    if str(tautomers.status) != "Completed":
        inconclusive = True
        flags.append(flag("TAUTOMER_LIMIT", "WARNING", "structure", "Enumeration limit reached"))
    for reference in [compound.reference_smiles, compound.canonical_smiles]:
        if reference is None:
            continue
        try:
            reference_mol = prepare(reference).parent
        except ValueError:
            flags.append(
                flag(
                    "INVALID_REFERENCE",
                    "CRITICAL_ARTIFACT",
                    "structure",
                    "Reference or asserted canonical SMILES is invalid",
                )
            )
            continue
        if Chem.MolToSmiles(reference_mol) == canonical:
            continue
        # Compare stereo first: tautomer canonicalization can erase stereo at affected atoms.
        same_graph = Chem.MolToSmiles(reference_mol, isomericSmiles=False) == Chem.MolToSmiles(
            parent, isomericSmiles=False
        )
        has_isotopes = any(a.GetIsotope() for a in parent.GetAtoms()) or any(
            a.GetIsotope() for a in reference_mol.GetAtoms()
        )
        if same_graph:
            flags.append(
                flag(
                    "ISOTOPE_IDENTITY_MISMATCH" if has_isotopes else "STEREO_MISMATCH",
                    "CRITICAL_ARTIFACT",
                    "structure",
                    "Same connectivity but different or missing stereochemical/isotope specification",
                )
            )
        elif Chem.MolToSmiles(enumerator.Canonicalize(reference_mol)) == Chem.MolToSmiles(
            enumerator.Canonicalize(parent)
        ):
            flags.append(
                flag(
                    "TAUTOMER_EQUIVALENT",
                    "INFO",
                    "structure",
                    "Structures belong to the same enumerated tautomer family",
                )
            )
        else:
            flags.append(
                flag(
                    "STRUCTURE_IDENTITY_MISMATCH",
                    "CRITICAL_ARTIFACT",
                    "structure",
                    "Reference and proposed structures disagree",
                )
            )
    centers = cast(Any, Chem).FindMolChiralCenters(parent, includeUnassigned=True)
    assigned = sum(label != "?" for _, label in centers)
    metrics["assigned_stereocenters"] = float(assigned)
    if compound.expected_stereocenters > assigned:
        flags.append(
            flag(
                "MISSING_STEREOCHEMISTRY",
                "CRITICAL_ARTIFACT",
                "structure",
                "Fewer assigned tetrahedral centers than explicitly required",
                expected=float(compound.expected_stereocenters),
                assigned=float(assigned),
            )
        )
    elif any(label == "?" for _, label in centers):
        inconclusive = True
        flags.append(
            flag(
                "UNSPECIFIED_STEREO",
                "WARNING",
                "structure",
                "Potential tetrahedral stereocenters lack assignments",
            )
        )
    for entry in _catalog().GetMatches(parent):
        flags.append(
            flag(
                "STRUCTURAL_ALERT",
                "WARNING",
                "structure",
                "PAINS/Brenk substructure alert is a risk screen, not proof of assay interference",
                alert=str(entry.GetDescription()),
            )
        )
    logp = float(cast(Any, Crippen).MolLogP(parent))
    rings = int(rdMolDescriptors.CalcNumAromaticRings(parent))
    metrics.update(clogp=logp, aromatic_rings=float(rings))
    if logp > 4.5 and rings >= 3:
        flags.append(
            flag(
                "AGGREGATION_RISK",
                "WARNING",
                "structure",
                "Lipophilicity/aromaticity heuristic warrants detergent and concentration controls",
                clogp=logp,
                aromatic_rings=float(rings),
            )
        )
    checks.extend(
        [
            "salt_reconciliation",
            "identity_assertions",
            "stereochemistry",
            "tautomer_family",
            "pains_brenk",
            "aggregation_heuristic",
        ]
    )
    return ModalityResult(
        modality="structure",
        flags=flags,
        corrected_metrics=metrics,
        checks_performed=checks,
        inconclusive=inconclusive,
    )
