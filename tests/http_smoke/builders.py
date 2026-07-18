"""
Molecule-graph builders for the HTTP smoke tests.

Builds atom/bond JSON payloads in the project's explicit-hydrogen convention
(see examples/README.md for the format). These are intentionally small,
composable helpers rather than a general SMILES-style parser - just enough
to construct rings, chains, and simple substituents by hand.
"""

VALENCY = {"C": 4, "O": 2, "N": 3, "F": 1, "Cl": 1, "Br": 1, "I": 1, "S": 2, "H": 1}


def new_atoms():
    return []


def add_atom(atoms, element):
    idx = len(atoms)
    atoms.append({"element": element, "bonds": []})
    return idx


def bond(atoms, a, b, order=1):
    atoms[a]["bonds"].append({"to": b, "order": order})
    atoms[b]["bonds"].append({"to": a, "order": order})


def carbon_ring(size, bond_orders=None):
    """Build an all-carbon ring of `size` atoms (indices 0..size-1).
    bond_orders: optional list of length `size`, order of bond ring[i]->ring[i+1].
    Defaults to all single bonds (cycloalkane)."""
    atoms = [{"element": "C", "bonds": []} for _ in range(size)]
    orders = bond_orders or [1] * size
    for i in range(size):
        j = (i + 1) % size
        bond(atoms, i, j, orders[i])
    return atoms


def benzene_ring(pattern="A"):
    """6-carbon ring with alternating Kekule bond orders.
    pattern A: 2,1,2,1,2,1 ; pattern B: 1,2,1,2,1,2 (equivalent structures)."""
    orders = [2, 1, 2, 1, 2, 1] if pattern == "A" else [1, 2, 1, 2, 1, 2]
    return carbon_ring(6, orders)


def attach_chain(atoms, attach_to, length):
    """Attach a straight (unbranched) alkyl chain of `length` carbons to `attach_to`.
    Returns the list of new carbon indices, chain[0] is bonded to attach_to."""
    prev = attach_to
    chain = []
    for _ in range(length):
        c = add_atom(atoms, "C")
        bond(atoms, prev, c)
        chain.append(c)
        prev = c
    return chain


def attach_atom(atoms, attach_to, element, order=1):
    """Attach a single heteroatom (Cl, Br, F, I, or O for =O) to `attach_to`."""
    idx = add_atom(atoms, element)
    bond(atoms, attach_to, idx, order)
    return idx


def attach_branch(atoms, attach_to, kind):
    """Attach a branched alkyl substituent: 'isopropyl' or 'tertbutyl'."""
    c1 = add_atom(atoms, "C")
    bond(atoms, attach_to, c1)
    n_methyls = 2 if kind == "isopropyl" else 3
    for _ in range(n_methyls):
        m = add_atom(atoms, "C")
        bond(atoms, c1, m)
    return c1


def attach_hydroxymethyl(atoms, attach_to):
    """Attach -CH2OH (branch-tip functional group, Phase 4)."""
    c1 = add_atom(atoms, "C")
    bond(atoms, attach_to, c1)
    o = add_atom(atoms, "O")
    bond(atoms, c1, o)
    return c1


def attach_chloromethyl(atoms, attach_to):
    """Attach -CH2Cl (branch-tip functional group, Phase 4)."""
    c1 = add_atom(atoms, "C")
    bond(atoms, attach_to, c1)
    cl = add_atom(atoms, "Cl")
    bond(atoms, c1, cl)
    return c1


def attach_carboxyl(atoms, attach_to):
    """Attach -COOH."""
    c = add_atom(atoms, "C")
    bond(atoms, attach_to, c)
    attach_atom(atoms, c, "O", order=2)
    attach_atom(atoms, c, "O")
    return c


def attach_acyl(atoms, attach_to, tail_len=0, terminal=None, terminal_order=1):
    """Attach -C(=O)R : tail_len carbons after the carbonyl C, and/or a
    terminal heteroatom on the carbonyl C (Cl for -COCl, N for -CONH2)."""
    c = add_atom(atoms, "C")
    bond(atoms, attach_to, c)
    attach_atom(atoms, c, "O", order=2)
    if tail_len:
        attach_chain(atoms, c, tail_len)
    if terminal:
        attach_atom(atoms, c, terminal, order=terminal_order)
    return c


def attach_nitrile(atoms, attach_to):
    """Attach -C≡N."""
    c = add_atom(atoms, "C")
    bond(atoms, attach_to, c)
    attach_atom(atoms, c, "N", order=3)
    return c


def attach_vinyl(atoms, attach_to):
    """Attach -CH=CH2."""
    c1 = add_atom(atoms, "C")
    bond(atoms, attach_to, c1)
    c2 = add_atom(atoms, "C")
    bond(atoms, c1, c2, 2)
    return c1


def ring_with_chain_listed_first(ring_size, chain_length, ring_bond_orders=None):
    """Build a ring with a straight alkyl chain attached at ring position 0, but with
    the chain's atoms placed BEFORE the ring atoms in the atom list. Regression-tests
    FindRing()'s anchor selection: a chain carbon can have the same C-neighbour count
    as a true ring carbon, so the anchor search must not just take the first match by
    index order (see Algorithm-Extension-Plan.md, Phase 2, Prerequisite fix 1)."""
    atoms = []
    chain_atoms = []
    prev = None
    for _ in range(chain_length):
        c = add_atom(atoms, "C")
        if prev is not None:
            bond(atoms, prev, c)
        chain_atoms.append(c)
        prev = c

    ring_offset = len(atoms)
    ring = carbon_ring(ring_size, ring_bond_orders)
    for a in ring:
        a["bonds"] = [{"to": b["to"] + ring_offset, "order": b["order"]} for b in a["bonds"]]
    atoms.extend(ring)

    if chain_atoms:
        bond(atoms, chain_atoms[-1], ring_offset)  # attach chain tip to ring position 0
    return atoms


def fill_hydrogens(atoms):
    """Append explicit H atoms so every atom reaches its valency. Mutates in place."""
    n = len(atoms)
    for i in range(n):
        element = atoms[i]["element"]
        used = sum(b["order"] for b in atoms[i]["bonds"])
        needed = VALENCY[element] - used
        for _ in range(needed):
            h = add_atom(atoms, "H")
            bond(atoms, i, h)


def payload(atoms, spec="AllGroups"):
    """Fill hydrogens and wrap into the API request body."""
    fill_hydrogens(atoms)
    return {"atoms": atoms, "specificationSet": spec}
