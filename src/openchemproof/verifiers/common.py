"""RDKit normalization and reporting shared by deterministic engines."""

from collections import Counter
from dataclasses import dataclass
from typing import Any, cast

from rdkit import Chem
from rdkit.Chem import SaltRemover, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

from openchemproof.schema.models import AuditFlag, CompoundSpec, Modality


@dataclass(frozen=True)
class PreparedMolecule:
    original: Chem.Mol
    parent: Chem.Mol
    stripped_fragments: tuple[str, ...]
    unresolved_fragments: bool


def prepare(smiles: str) -> PreparedMolecule:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError("Invalid or empty SMILES")
    if mol.GetNumAtoms() > 1000:
        raise ValueError("v0.1 supports at most 1000 atoms")
    stripped = cast(Any, SaltRemover).SaltRemover().StripMol(mol, dontRemoveEverything=True)
    original_fragments = Counter(Chem.MolToSmiles(m) for m in Chem.GetMolFrags(mol, asMols=True))
    remaining = Counter(Chem.MolToSmiles(m) for m in Chem.GetMolFrags(stripped, asMols=True))
    removed = tuple(sorted((original_fragments - remaining).elements()))
    # Do not silently pick the largest of two unknown organic components.
    unresolved = len(Chem.GetMolFrags(stripped)) != 1
    parent = rdMolStandardize.Uncharger().uncharge(stripped)
    return PreparedMolecule(mol, parent, removed, unresolved)


def exact_mass(smiles: str) -> float:
    return float(rdMolDescriptors.CalcExactMolWt(prepare(smiles).parent))


def compound_from_smiles(compound_id: str, smiles: str) -> CompoundSpec:
    prepared = prepare(smiles)
    parent = prepared.parent
    key = str(cast(Any, Chem).MolToInchiKey(parent))
    return CompoundSpec(
        compound_id=compound_id,
        proposed_smiles=smiles,
        canonical_smiles=str(Chem.MolToSmiles(parent)),
        inchikey=key,
        inchi_connectivity_layer=key[:14],
        molecular_formula=str(rdMolDescriptors.CalcMolFormula(parent)),
        monoisotopic_mass=float(rdMolDescriptors.CalcExactMolWt(parent)),
        expected_stereocenters=len(cast(Any, Chem).FindMolChiralCenters(parent, includeUnassigned=False)),
    )


def flag(code: str, severity: str, modality: Modality, message: str, **evidence: float | str) -> AuditFlag:
    # model_validate is used so even dynamically assembled flags obey the wire schema.
    return AuditFlag.model_validate(
        dict(
            code=code,
            severity=severity,
            modality=modality,
            message=message,
            quantitative_evidence=evidence,
        )
    )
