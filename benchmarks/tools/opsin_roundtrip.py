#!/usr/bin/env python3
"""
Structural round-trip validation of engine names, adjudicated by OPSIN.

The engine turns a molecular graph into an IUPAC name. OPSIN (Open Parser for
Systematic IUPAC Nomenclature, Lowe et al., University of Cambridge) turns an
IUPAC name back into a structure. Composing them closes a loop:

    graph  --engine-->  name  --OPSIN-->  graph'

If canonical(graph) == canonical(graph'), the name unambiguously denotes the
molecule it was generated from. This is a far stronger and more objective test
than string-matching against a reference name, because it is immune to
legitimate differences in naming convention -- "4-(1-methylethyl)benzaldehyde"
and "4-propan-2-ylbenzaldehyde" are the same name to this test, as they should
be, while a name that quietly denotes the wrong structure fails no matter how
plausible it reads.

Comparison is constitutional only: stereochemistry is stripped (the engine has
no stereo model) and nitro groups are normalised to one tautomeric form, since
OPSIN emits the charge-separated [N+](=O)[O-] and the engine's input model uses
the neutral hypervalent N(=O)=O.

Usage:
    python opsin_roundtrip.py --in <engine-output.jsonl> --out <results.jsonl>
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile

from rdkit import Chem, RDLogger

from smiles_to_graph import atoms_to_mol

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
OPSIN_JAR = os.path.join(HERE, "opsin-2.9.0.jar")


def run_opsin(names: list[str]) -> list[str | None]:
    """Batch-resolve names to SMILES. Returns None for names OPSIN can't parse."""
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        for n in names:
            fh.write(n.replace("\n", " ") + "\n")
        path = fh.name
    try:
        proc = subprocess.run(
            ["java", "-jar", OPSIN_JAR, "-o", "smi", path],
            capture_output=True, text=True, timeout=1800,
        )
        lines = proc.stdout.split("\n")
        # OPSIN emits exactly one line per input, blank where it failed.
        out = [ln.strip() or None for ln in lines[:len(names)]]
        while len(out) < len(names):
            out.append(None)
        return out
    finally:
        os.unlink(path)


def normalise(mol: Chem.Mol | None) -> str | None:
    """Canonical constitutional SMILES: no stereo, nitro in one fixed form."""
    if mol is None:
        return None
    mol = Chem.RWMol(mol)
    # Collapse charge-separated nitro to the neutral hypervalent form so both
    # sides of the comparison agree on one representation.
    patt = Chem.MolFromSmarts("[N+;$([N+](=O)[O-])]")
    for match in mol.GetSubstructMatches(patt):
        n = mol.GetAtomWithIdx(match[0])
        for nb in n.GetNeighbors():
            if nb.GetSymbol() == "O" and nb.GetFormalCharge() == -1:
                nb.SetFormalCharge(0)
                mol.GetBondBetweenAtoms(n.GetIdx(), nb.GetIdx()).SetBondType(
                    Chem.BondType.DOUBLE)
                n.SetFormalCharge(0)
                n.SetNoImplicit(True)
                break
    m = mol.GetMol()
    try:
        Chem.SanitizeMol(m, Chem.SanitizeFlags.SANITIZE_ALL ^
                         Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    except Exception:
        return None
    Chem.RemoveStereochemistry(m)
    try:
        return Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(m)) or m)
    except Exception:
        return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    args = ap.parse_args()

    rows = []
    for line in open(args.inp):
        d = json.loads(line)
        if d.get("status") == "named" and d.get("names"):
            rows.append(d)

    # The engine may emit several equally-valid names; the round-trip passes if
    # ANY of them resolves to the input structure.
    flat, owner = [], []
    for i, d in enumerate(rows):
        for nm in d["names"]:
            flat.append(nm)
            owner.append(i)

    print(f"resolving {len(flat)} names for {len(rows)} molecules through OPSIN...",
          file=sys.stderr)
    smis = run_opsin(flat)

    per_mol: dict[int, list[str | None]] = {}
    for idx, smi in zip(owner, smis):
        per_mol.setdefault(idx, []).append(smi)

    n_pass = n_fail = n_unparsed = 0
    with open(args.out, "w") as fh:
        for i, d in enumerate(rows):
            target = normalise(atoms_to_mol(d["atoms"]))
            got = per_mol.get(i, [])
            cands = [normalise(Chem.MolFromSmiles(s)) if s else None for s in got]
            if target is not None and target in cands:
                verdict, n_pass = "match", n_pass + 1
            elif all(c is None for c in cands):
                verdict, n_unparsed = "opsin_could_not_parse", n_unparsed + 1
            else:
                verdict, n_fail = "structure_mismatch", n_fail + 1
            fh.write(json.dumps({
                "id": d.get("id"), "names": d["names"],
                "reference_name": d.get("reference_name"),
                "target_smiles": target, "roundtrip_smiles": cands,
                "verdict": verdict,
            }) + "\n")

    total = len(rows)
    print(f"\nround-trip over {total} named molecules")
    print(f"  match                 {n_pass:6d}  {100*n_pass/total:6.2f}%")
    print(f"  structure mismatch    {n_fail:6d}  {100*n_fail/total:6.2f}%")
    print(f"  OPSIN could not parse {n_unparsed:6d}  {100*n_unparsed/total:6.2f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
