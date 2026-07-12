"""
Named test cases for the HTTP smoke harness. Each case maps to:
  {"payload": <request body>, "expect": {...}}

expect["type"]:
  "exact"    - HTTP 200 and names == [expect["name"]]
  "reject"   - HTTP 4xx (a clean ChemistryError/ValidationError, not a crash)
  "no_crash" - any well-formed JSON response (200 or error); just proves the
               server didn't fall over on adversarial/malformed input

These cases track what's actually implemented as of Phase 2 (aromatic compounds).
See Algorithm-Extension-Plan.md for the phase-by-phase design this exercises.
"""

from builders import (
    new_atoms, add_atom, bond, carbon_ring, benzene_ring,
    attach_chain, attach_atom, attach_branch, attach_hydroxymethyl,
    attach_chloromethyl, ring_with_chain_listed_first, payload,
)


def exact(p, name):
    return {"payload": p, "expect": {"type": "exact", "name": name}}


def reject(p):
    return {"payload": p, "expect": {"type": "reject"}}


def no_crash(p):
    return {"payload": p, "expect": {"type": "no_crash"}}


CASES = {}

# ── Regression: existing chain-based naming (unaffected by cyclic/aromatic work) ──
def _chain(n):
    a = new_atoms()
    c0 = add_atom(a, "C")
    attach_chain(a, c0, n - 1)
    return a

CASES["ethane"] = exact(payload(_chain(2)), "ethane")
CASES["propane"] = exact(payload(_chain(3)), "propane")

def _chloropropane_2():
    a = new_atoms()
    c0 = add_atom(a, "C")
    c1 = add_atom(a, "C")
    bond(a, c0, c1)
    c2 = add_atom(a, "C")
    bond(a, c1, c2)
    attach_atom(a, c1, "Cl")
    return a

CASES["2_chloropropane"] = exact(payload(_chloropropane_2()), "2-chloropropane")

# ── Phase 1: cyclic (cycloalkanes) ──
CASES["cyclopropane"] = exact(payload(carbon_ring(3)), "cyclopropane")
CASES["cyclohexane"] = exact(payload(carbon_ring(6)), "cyclohexane")

def _methylcyclohexane():
    a = carbon_ring(6)
    attach_chain(a, 0, 1)
    return a

CASES["methylcyclohexane"] = exact(payload(_methylcyclohexane()), "methylcyclohexane")

CASES["ethylcyclohexane_chain_listed_first"] = exact(
    payload(ring_with_chain_listed_first(6, 2)), "ethylcyclohexane"
)  # regression for the FindRing() anchor-selection bug (Phase 2, prereq fix 1)

def _dimethylcyclohexane(pos_a, pos_b):
    a = carbon_ring(6)
    attach_chain(a, pos_a, 1)
    attach_chain(a, pos_b, 1)
    return a

CASES["12_dimethylcyclohexane"] = exact(payload(_dimethylcyclohexane(0, 1)), "1,2-dimethylcyclohexane")
CASES["14_dimethylcyclohexane"] = exact(payload(_dimethylcyclohexane(0, 3)), "1,4-dimethylcyclohexane")

def _chlorocyclohexane():
    a = carbon_ring(6)
    attach_atom(a, 0, "Cl")
    return a

CASES["chlorocyclohexane"] = exact(payload(_chlorocyclohexane()), "chlorocyclohexane")

# Branched/functionalised ring substituents are explicitly rejected as of Phase 2
# (see ValidateRingSubstituents in Algorithm-Extension-Plan.md) - previously this
# silently misnamed as "propylcyclohexane", ignoring the branching.
def _isopropylcyclohexane():
    a = carbon_ring(6)
    attach_branch(a, 0, "isopropyl")
    return a

CASES["isopropylcyclohexane_branched_sub"] = reject(payload(_isopropylcyclohexane()))

# ── Phase 4: branch-tip functional groups ──
def _hydroxymethylpropane():
    a = new_atoms()
    c0 = add_atom(a, "C")
    c1 = add_atom(a, "C")
    bond(a, c0, c1)
    c2 = add_atom(a, "C")
    bond(a, c1, c2)
    attach_hydroxymethyl(a, c1)
    return a

CASES["2_hydroxymethylpropane"] = exact(payload(_hydroxymethylpropane()), "2-(hydroxymethyl)propane")

def _chloromethylbutane():
    a = new_atoms()
    c0 = add_atom(a, "C")
    c1 = add_atom(a, "C")
    bond(a, c0, c1)
    c2 = add_atom(a, "C")
    bond(a, c1, c2)
    c3 = add_atom(a, "C")
    bond(a, c2, c3)
    attach_chloromethyl(a, c2)
    return a

# Locant is 2, not 3: numbering from the other end of butane gives the lower
# locant, per IUPAC rules - verified against the live API (see Algorithm-Extension-Plan.md).
CASES["2_chloromethylbutane"] = exact(payload(_chloromethylbutane()), "2-(chloromethyl)butane")

# ── Phase 2: aromatic (benzene) ──
CASES["benzene_pattern_a"] = exact(payload(benzene_ring("A")), "benzene")
CASES["benzene_pattern_b"] = exact(payload(benzene_ring("B")), "benzene")

def _methylbenzene():
    a = benzene_ring("A")
    attach_chain(a, 0, 1)
    return a

CASES["methylbenzene"] = exact(payload(_methylbenzene()), "methylbenzene")

CASES["ethylbenzene_chain_listed_first"] = exact(
    payload(ring_with_chain_listed_first(6, 2, [2, 1, 2, 1, 2, 1])), "ethylbenzene"
)

def _chlorobenzene():
    a = benzene_ring("A")
    attach_atom(a, 0, "Cl")
    return a

CASES["chlorobenzene"] = exact(payload(_chlorobenzene()), "chlorobenzene")

def _dimethylbenzene(pos_a, pos_b):
    a = benzene_ring("A")
    attach_chain(a, pos_a, 1)
    attach_chain(a, pos_b, 1)
    return a

CASES["12_dimethylbenzene"] = exact(payload(_dimethylbenzene(0, 1)), "1,2-dimethylbenzene")
CASES["14_dimethylbenzene"] = exact(payload(_dimethylbenzene(0, 3)), "1,4-dimethylbenzene")

def _chloromethylbenzene():
    a = benzene_ring("A")
    attach_atom(a, 0, "Cl")
    attach_chain(a, 1, 1)
    return a

CASES["1chloro2methylbenzene"] = exact(payload(_chloromethylbenzene()), "1-chloro-2-methylbenzene")

# ── Phase 2 negative cases: must reject, not silently misname ──
def _trimethylbenzene():
    a = benzene_ring("A")
    attach_chain(a, 0, 1)
    attach_chain(a, 2, 1)
    attach_chain(a, 4, 1)
    return a

CASES["trimethylbenzene_should_reject"] = reject(payload(_trimethylbenzene()))
CASES["cyclohexene_should_reject"] = reject(payload(carbon_ring(6, [2, 1, 1, 1, 1, 1])))

def _benzyl_alcohol():
    a = benzene_ring("A")
    attach_hydroxymethyl(a, 0)
    return a

CASES["benzyl_alcohol_should_reject"] = reject(payload(_benzyl_alcohol()))

def _styrene():
    a = benzene_ring("A")
    c6 = add_atom(a, "C")
    bond(a, 0, c6)
    c7 = add_atom(a, "C")
    bond(a, c6, c7, 2)
    return a

CASES["styrene_should_reject"] = reject(payload(_styrene()))

def _isopropylbenzene():
    a = benzene_ring("A")
    attach_branch(a, 0, "isopropyl")
    return a

CASES["isopropylbenzene_should_reject"] = reject(payload(_isopropylbenzene()))

def _cyclohexanone():
    a = carbon_ring(6)
    attach_atom(a, 0, "O", order=2)
    return a

CASES["cyclohexanone_should_reject"] = reject(payload(_cyclohexanone()))

def _fused_ring():
    a = benzene_ring("A")
    c6 = add_atom(a, "C")
    bond(a, 0, c6)
    c7 = add_atom(a, "C")
    bond(a, c6, c7)
    bond(a, c7, 1)  # closes a second ring fused to the benzene ring
    return a

CASES["fused_ring_should_reject"] = reject(payload(_fused_ring()))

# ── Adversarial / malformed input: must not crash the server ──
CASES["empty_atoms_list"] = no_crash({"atoms": [], "specificationSet": "AllGroups"})
CASES["out_of_range_bond_index"] = no_crash({
    "atoms": [{"element": "C", "bonds": [{"to": 999, "order": 1}]}],
    "specificationSet": "AllGroups",
})
CASES["overvalent_carbon"] = no_crash({
    "atoms": [
        {"element": "C", "bonds": [{"to": i, "order": 1} for i in range(1, 6)]},
        *[{"element": "H", "bonds": [{"to": 0, "order": 1}]} for _ in range(5)],
    ],
    "specificationSet": "AllGroups",
})
