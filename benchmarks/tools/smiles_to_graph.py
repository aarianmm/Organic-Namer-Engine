#!/usr/bin/env python3
"""Convert SMILES strings to the Organic-Namer-Engine's molecule JSON format.

Engine format (see examples/ethanol.json in the repo root):

    {"atoms": [ {"element": "C", "bonds": [{"to": 1, "order": 1}, ...]}, ... ]}

Rules enforced here:
  - Bonds are listed symmetrically on both endpoints.
  - `to` is a 0-based index into the atoms array.
  - `order` is 1, 2 or 3.
  - Hydrogens are omitted entirely (the engine derives them from valence).
  - Only elements the engine understands are allowed: C, H, O, N, F, Cl, Br, I, S.
  - Aromatic rings are written Kekule-style (alternating single/double bonds).
  - Charge-separated nitro groups ([N+](=O)[O-]) are rewritten to the neutral
    hypervalent form N(=O)=O, which is the only form the engine accepts since
    it has no formal-charge field.
  - Any other charge, any isotope, any radical, or any unsupported element
    makes the molecule unconvertible.
  - Multi-fragment SMILES (containing '.') are unconvertible.
  - Stereochemistry is stripped/ignored (the engine has no stereo model).

Usage as a CLI:

    python smiles_to_graph.py --in in.jsonl --out out.jsonl --fail-out unconvertible.jsonl
    python smiles_to_graph.py --selftest

Usage as a library:

    from smiles_to_graph import smiles_to_atoms, atoms_to_smiles
    atoms, err = smiles_to_atoms("CCO")
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any, Optional

from rdkit import Chem
from rdkit.Chem import rdmolops

# Elements the naming engine knows how to handle.
SUPPORTED_ELEMENTS = {"C", "H", "O", "N", "F", "Cl", "Br", "I", "S"}

# SMARTS for the charge-separated nitro group: a nitrogen bearing a formal
# +1 charge, double-bonded to one O and single-bonded to an O- .
NITRO_CHARGE_SEPARATED = Chem.MolFromSmarts("[N+](=O)[O-]")


class ConversionError(Exception):
    """Raised internally to signal a structured, reported failure reason."""

    def __init__(self, reason: str):
        super().__init__(reason)
        self.reason = reason


def _find_nitro_atoms(mol: Chem.Mol) -> tuple[set[int], dict[int, int]]:
    """Find charge-separated nitro groups [N+](=O)[O-] in `mol`.

    Returns (all_nitro_atom_indices, n_idx_to_o_minus_idx) where the latter
    maps each nitro nitrogen's atom index to the index of its single-bonded
    (charge -1) oxygen -- the bond that needs to be promoted to a double
    bond in the engine's neutral hypervalent output representation.

    Note: RDKit cannot represent the neutral hypervalent form
    (N with bonds 1,2,2 and no charge) as a valid sanitized mol -- its
    valence model only allows this pattern via the charge-separated form.
    So rather than rewriting the RDKit mol itself, we detect the pattern
    here and apply the bond-order override only when emitting the final
    JSON atom/bond list (see smiles_to_atoms).
    """
    all_idxs: set[int] = set()
    n_to_ominus: dict[int, int] = {}
    for match in mol.GetSubstructMatches(NITRO_CHARGE_SEPARATED):
        n_idx, o_double_idx, o_minus_idx = match
        all_idxs.update(match)
        n_to_ominus[n_idx] = o_minus_idx
    return all_idxs, n_to_ominus


def smiles_to_atoms(smiles: str) -> tuple[Optional[list[dict[str, Any]]], Optional[str]]:
    """Convert a SMILES string to the engine's heavy-atom-only JSON atom list.

    Returns (atoms, None) on success, or (None, reason) on failure, where
    `reason` is a short machine-readable string describing why the molecule
    could not be converted.
    """
    if smiles is None:
        return None, "empty_smiles"
    smiles = smiles.strip()
    if not smiles:
        return None, "empty_smiles"

    if "." in smiles:
        return None, "multi_fragment"

    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
    except Exception:
        mol = None
    if mol is None:
        return None, "rdkit_parse_failed"

    # Identify charge-separated nitro groups; their N/O atoms are exempt
    # from the generic charge check below because we neutralise them in
    # the output (as the neutral hypervalent form) rather than rejecting.
    nitro_atom_idxs, _ = _find_nitro_atoms(mol)

    # Reject any charges/isotopes/radicals not accounted for by a
    # recognised nitro group.
    for atom in mol.GetAtoms():
        if atom.GetIdx() in nitro_atom_idxs:
            continue
        if atom.GetFormalCharge() != 0:
            return None, "charged_atom"
        if atom.GetIsotope() != 0:
            return None, "isotope"
        if atom.GetNumRadicalElectrons() != 0:
            return None, "radical"

    # Reject unsupported elements (heavy atoms only; hydrogens are fine but
    # will be stripped below anyway).
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in SUPPORTED_ELEMENTS:
            return None, f"unsupported_element:{symbol}"

    # Kekulize: turn aromatic rings into explicit alternating single/double
    # bonds and clear aromatic flags, so a benzene ring emits 1,2,1,2,1,2
    # rather than all-single bonds.
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception as exc:
        return None, f"kekulize_failed: {exc}"

    # Strip explicit hydrogens; the engine derives them from valence.
    try:
        mol = rdmolops.RemoveHs(mol, updateExplicitCount=False, sanitize=False)
    except Exception:
        pass

    # Re-locate nitro groups on the (possibly re-indexed) H-stripped mol so
    # we can override their N-O bond order to 2 (the neutral hypervalent
    # form) when building the output below.
    _, nitro_n_to_ominus = _find_nitro_atoms(mol)
    nitro_bond_overrides: set[tuple[int, int]] = set()
    for n_idx, o_minus_idx in nitro_n_to_ominus.items():
        key = (min(n_idx, o_minus_idx), max(n_idx, o_minus_idx))
        nitro_bond_overrides.add(key)

    bond_order_map = {
        Chem.BondType.SINGLE: 1,
        Chem.BondType.DOUBLE: 2,
        Chem.BondType.TRIPLE: 3,
    }

    atoms: list[dict[str, Any]] = []
    heavy_atoms = [a for a in mol.GetAtoms() if a.GetSymbol() != "H"]
    # Build an index remap in case RemoveHs left any H atoms in place
    # (shouldn't normally happen, but guard against it) or reordered indices.
    old_to_new: dict[int, int] = {}
    for new_idx, atom in enumerate(heavy_atoms):
        old_to_new[atom.GetIdx()] = new_idx

    for atom in heavy_atoms:
        element = atom.GetSymbol()
        bonds_out = []
        for bond in atom.GetBonds():
            other = bond.GetOtherAtom(atom)
            if other.GetSymbol() == "H":
                continue  # heavy atoms only
            order = bond_order_map.get(bond.GetBondType())
            if order is None:
                return None, f"unsupported_bond_type:{bond.GetBondType()}"
            other_new_idx = old_to_new.get(other.GetIdx())
            if other_new_idx is None:
                return None, "internal_index_error"
            # Apply the nitro neutral-hypervalent override: the N-O(-) bond
            # (currently a single bond) becomes a double bond.
            bond_key = (min(atom.GetIdx(), other.GetIdx()), max(atom.GetIdx(), other.GetIdx()))
            if bond_key in nitro_bond_overrides:
                order = 2
            bonds_out.append({"to": other_new_idx, "order": order})
        atoms.append({"element": element, "bonds": bonds_out})

    if not atoms:
        return None, "no_heavy_atoms"

    return atoms, None


def atoms_to_mol(atoms: list[dict[str, Any]]) -> Chem.Mol:
    """Rebuild an RDKit mol (heavy atoms + explicit bonds) from the engine's
    atom list, WITHOUT adding hydrogens. Symmetric duplicate bonds (i lists
    j, j lists i) are collapsed to a single bond.

    The engine's neutral hypervalent nitro form (N with bonds 1,2,2 to a
    substituent and two terminal O atoms, no charge field) cannot be
    sanitised by RDKit directly -- its valence model requires the
    charge-separated form for this pattern. So any such pattern found here
    is converted back to charge-separated ([N+](=O)[O-]) before returning,
    mirroring the inverse of the override applied in smiles_to_atoms.
    """
    rw = Chem.RWMol()
    for atom_dict in atoms:
        a = Chem.Atom(atom_dict["element"])
        rw.AddAtom(a)

    bond_type_map = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
    seen: set[tuple[int, int]] = set()
    for i, atom_dict in enumerate(atoms):
        for bond in atom_dict.get("bonds", []):
            j = bond["to"]
            key = (min(i, j), max(i, j))
            if key in seen:
                continue
            seen.add(key)
            rw.AddBond(i, j, bond_type_map[bond["order"]])

    # Detect and fix neutral hypervalent nitro groups: an N with exactly
    # three bonds, orders [1, 2, 2], where both order-2 neighbours are
    # terminal O atoms (degree 1).
    for i, atom_dict in enumerate(atoms):
        if atom_dict["element"] != "N":
            continue
        bonds = atom_dict.get("bonds", [])
        if len(bonds) != 3:
            continue
        if sorted(b["order"] for b in bonds) != [1, 2, 2]:
            continue
        o_candidates = [b["to"] for b in bonds if b["order"] == 2]
        if any(atoms[o_idx]["element"] != "O" or len(atoms[o_idx].get("bonds", [])) != 1
               for o_idx in o_candidates):
            continue

        n_atom = rw.GetAtomWithIdx(i)
        n_atom.SetFormalCharge(1)
        o_idx = o_candidates[0]
        bond = rw.GetBondBetweenAtoms(i, o_idx)
        bond.SetBondType(Chem.BondType.SINGLE)
        rw.GetAtomWithIdx(o_idx).SetFormalCharge(-1)

    return rw.GetMol()


def atoms_to_smiles(atoms: list[dict[str, Any]]) -> str:
    """Rebuild a molecule from the engine-format atom list and return its
    canonical SMILES, adding explicit hydrogens via sanitisation (the
    inverse of the heavy-atom-only conversion done by smiles_to_atoms).
    """
    mol = atoms_to_mol(atoms)
    mol = mol.GetMol() if hasattr(mol, "GetMol") else mol
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


# ---------------------------------------------------------------------------
# CLI plumbing
# ---------------------------------------------------------------------------


def _run_conversion(in_path: str, out_path: str, fail_out_path: str) -> None:
    n_ok = 0
    n_fail = 0
    with open(in_path, "r", encoding="utf-8") as fin, \
            open(out_path, "w", encoding="utf-8") as fout, \
            open(fail_out_path, "w", encoding="utf-8") as ffail:
        for line_no, line in enumerate(fin, 1):
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                sys.stderr.write(f"line {line_no}: invalid JSON ({exc}), skipping\n")
                continue

            # Accept the lowercase "smiles" key from the generic pipeline
            # contract, as well as PubChem's raw property casings ("SMILES",
            # falling back to the legacy "CanonicalSMILES").
            smiles = record.get("smiles") or record.get("SMILES") or record.get("CanonicalSMILES")
            atoms, err = smiles_to_atoms(smiles)
            if err is None:
                out_record = dict(record)
                out_record["atoms"] = atoms
                fout.write(json.dumps(out_record) + "\n")
                n_ok += 1
            else:
                fail_record = dict(record)
                fail_record["convert_error"] = err
                ffail.write(json.dumps(fail_record) + "\n")
                n_fail += 1

    print(f"Converted: {n_ok} ok, {n_fail} failed -> {out_path} / {fail_out_path}")


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

SELFTEST_MOLECULES: list[tuple[str, str]] = [
    ("methane", "C"),
    ("ethane", "CC"),
    ("propane", "CCC"),
    ("butane (linear alkane)", "CCCC"),
    ("2-methylpropane (branched alkane)", "CC(C)C"),
    ("2,2-dimethylpropane (branched alkane)", "CC(C)(C)C"),
    ("ethene (alkene)", "C=C"),
    ("propene (alkene)", "C=CC"),
    ("2-butene (alkene)", "CC=CC"),
    ("1,3-butadiene (diene)", "C=CC=C"),
    ("ethyne (alkyne)", "C#C"),
    ("propyne (alkyne)", "C#CC"),
    ("ethanol (alcohol)", "CCO"),
    ("propan-2-ol (alcohol)", "CC(O)C"),
    ("propan-1-ol (alcohol)", "CCCO"),
    ("propanone (ketone)", "CC(=O)C"),
    ("butanone (ketone)", "CCC(=O)C"),
    ("ethanal (aldehyde)", "CC=O"),
    ("propanal (aldehyde)", "CCC=O"),
    ("ethanoic acid (carboxylic acid)", "CC(=O)O"),
    ("propanoic acid (carboxylic acid)", "CCC(=O)O"),
    ("methyl ethanoate (ester)", "CC(=O)OC"),
    ("ethyl ethanoate (ester)", "CC(=O)OCC"),
    ("methoxymethane (ether)", "COC"),
    ("diethyl ether (ether)", "CCOCC"),
    ("methanamine (amine)", "CN"),
    ("ethanamine (amine)", "CCN"),
    ("dimethylamine (secondary amine)", "CNC"),
    ("trimethylamine (tertiary amine)", "CN(C)C"),
    ("ethanamide (amide)", "CC(=O)N"),
    ("N-methylethanamide (amide)", "CC(=O)NC"),
    ("acetonitrile (nitrile)", "CC#N"),
    ("propanenitrile (nitrile)", "CCC#N"),
    ("chloromethane (halide)", "CCl"),
    ("2-chloropropane (halide)", "CC(Cl)C"),
    ("bromoethane (halide)", "CCBr"),
    ("1-fluorobutane (halide)", "CCCCF"),
    ("iodoethane (halide)", "CCI"),
    ("nitromethane (nitro compound)", "C[N+](=O)[O-]"),
    ("nitrobenzene (aromatic nitro compound)", "c1ccccc1[N+](=O)[O-]"),
    ("benzene (aromatic)", "c1ccccc1"),
    ("toluene (substituted benzene)", "Cc1ccccc1"),
    ("phenol (aromatic alcohol)", "Oc1ccccc1"),
    ("aniline (aromatic amine)", "Nc1ccccc1"),
    ("chlorobenzene (substituted benzene)", "Clc1ccccc1"),
    ("cyclohexane (cycloalkane)", "C1CCCCC1"),
    ("cyclopentane (cycloalkane)", "C1CCCC1"),
    ("methylcyclohexane (cycloalkane)", "CC1CCCCC1"),
    ("cyclohexanol (cyclic alcohol)", "OC1CCCCC1"),
    ("thioethanol / ethanethiol (thiol)", "CCS"),
    ("dimethyl sulfide (sulfide)", "CSC"),
]


def _run_selftest() -> bool:
    passed = 0
    failed: list[tuple[str, str, str]] = []
    for name, smiles in SELFTEST_MOLECULES:
        try:
            canon_in = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
        except Exception as exc:
            failed.append((name, smiles, f"could not canonicalise input SMILES: {exc}"))
            continue

        atoms, err = smiles_to_atoms(smiles)
        if err is not None:
            failed.append((name, smiles, f"smiles_to_atoms failed: {err}"))
            continue

        try:
            canon_out = atoms_to_smiles(atoms)
        except Exception as exc:
            failed.append((name, smiles, f"atoms_to_smiles raised: {exc}"))
            continue

        if canon_out != canon_in:
            failed.append((name, smiles, f"round-trip mismatch: {canon_in!r} != {canon_out!r}"))
            continue

        passed += 1

    total = len(SELFTEST_MOLECULES)
    print(f"Self-test: {passed}/{total} round-trips passed")
    if failed:
        print("Failures:")
        for name, smiles, reason in failed:
            print(f"  - {name} ({smiles}): {reason}")
    return len(failed) == 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in", dest="in_path", help="Input JSONL file with {'smiles': ..., 'id': ...} records")
    parser.add_argument("--out", dest="out_path", help="Output JSONL file for convertible molecules")
    parser.add_argument("--fail-out", dest="fail_out_path", help="Output JSONL file for unconvertible molecules")
    parser.add_argument("--selftest", action="store_true", help="Run the built-in round-trip self-test and exit")
    args = parser.parse_args()

    if args.selftest:
        ok = _run_selftest()
        return 0 if ok else 1

    if not (args.in_path and args.out_path and args.fail_out_path):
        parser.error("--in, --out and --fail-out are all required unless --selftest is given")

    _run_conversion(args.in_path, args.out_path, args.fail_out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
