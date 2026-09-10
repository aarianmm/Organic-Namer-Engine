# benchmarks/tools

Python tooling for building a reference dataset to benchmark the Organic-Namer-Engine's
IUPAC naming against real compound data from PubChem.

## Setup

RDKit has no wheel for Python 3.14, so these tools use a Python 3.13 virtualenv.

```bash
/opt/homebrew/bin/python3.13 -m venv benchmarks/tools/.venv
benchmarks/tools/.venv/bin/pip install --quiet --upgrade pip rdkit requests
```

(If the 3.13 wheel is ever unavailable, fall back to `/opt/homebrew/bin/python3.11`.)

Verify it works:

```bash
benchmarks/tools/.venv/bin/python -c "from rdkit import Chem; print(Chem.MolFromSmiles('CCO'))"
```

## Scripts

### `smiles_to_graph.py`

Converts SMILES to the engine's molecule JSON format (`{"atoms": [...]}`, heavy atoms
only, hydrogens omitted, bonds listed symmetrically on both endpoints, aromatic rings
Kekulized to alternating single/double bonds, charge-separated nitro groups
`[N+](=O)[O-]` rewritten to the neutral hypervalent form the engine expects).

Molecules that can't be represented in the engine's format (any formal charge other
than a recognised nitro group, isotopes, radicals, unsupported elements, or
multi-fragment/salt SMILES) are reported as unconvertible with a reason, not guessed at.

Library use:

```python
from smiles_to_graph import smiles_to_atoms, atoms_to_smiles
atoms, err = smiles_to_atoms("CCO")
```

CLI:

```bash
# Convert a JSONL file of {"smiles": ..., "id": ...} (+ passthrough fields) records.
# Also accepts PubChem's raw "SMILES"/"CanonicalSMILES" field casing directly.
benchmarks/tools/.venv/bin/python benchmarks/tools/smiles_to_graph.py \
  --in benchmarks/data/pubchem_raw.jsonl \
  --out benchmarks/data/pubchem_graphs.jsonl \
  --fail-out benchmarks/data/pubchem_unconvertible.jsonl

# Run the built-in round-trip self-test (SMILES -> atoms -> SMILES, canonical
# forms must match) over 51 molecules spanning the major functional groups.
benchmarks/tools/.venv/bin/python benchmarks/tools/smiles_to_graph.py --selftest
```

### `fetch_pubchem.py`

Fetches compound records (IUPACName, SMILES, MolecularFormula, Title) from PubChem's
PUG-REST API by CID, in batches of 100, respecting PubChem's 5 requests/second rate
limit. Caches raw responses to a JSONL file so re-runs with `--resume` skip CIDs
already fetched. Records with no usable IUPAC name or SMILES are skipped.

```bash
# By CID range:
benchmarks/tools/.venv/bin/python benchmarks/tools/fetch_pubchem.py \
  --cids 1-6000 --out benchmarks/data/pubchem_raw.jsonl

# By CID file (one CID per line), resuming a previous run:
benchmarks/tools/.venv/bin/python benchmarks/tools/fetch_pubchem.py \
  --cid-file cids.txt --out benchmarks/data/pubchem_raw.jsonl --resume
```

## End-to-end pipeline

```bash
# 1. Set up the venv (see Setup above).

# 2. Fetch raw PubChem records for CIDs 1-6000 (low CIDs skew toward small,
#    well-characterised organic molecules -- a good benchmark population).
benchmarks/tools/.venv/bin/python benchmarks/tools/fetch_pubchem.py \
  --cids 1-6000 --out benchmarks/data/pubchem_raw.jsonl

# 3. Convert to the engine's molecule JSON format, splitting out anything
#    that can't be represented (salts, most charged species, unsupported
#    elements like P/B/Se/Si, etc).
benchmarks/tools/.venv/bin/python benchmarks/tools/smiles_to_graph.py \
  --in benchmarks/data/pubchem_raw.jsonl \
  --out benchmarks/data/pubchem_graphs.jsonl \
  --fail-out benchmarks/data/pubchem_unconvertible.jsonl
```

`benchmarks/data/pubchem_graphs.jsonl` is the dataset consumed by the benchmark runner:
each line is a PubChem property record plus an `"atoms"` field in the engine's format.
`benchmarks/data/pubchem_unconvertible.jsonl` holds everything that was skipped, each
with a `"convert_error"` reason for auditing.
