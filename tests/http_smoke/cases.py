"""
Named test cases for the HTTP smoke harness. Each case maps to:
  {"payload": <request body>, "expect": {...}}

expect["type"]:
  "exact"    - HTTP 200 and names == [expect["name"]]
  "reject"   - HTTP 4xx (a clean ChemistryError/ValidationError, not a crash)
  "no_crash" - any well-formed JSON response (200 or error); just proves the
               server didn't fall over on adversarial/malformed input

These cases track what's actually implemented as of Phase 3 (esters/ethers).
See Algorithm-Extension-Plan.md for the phase-by-phase design this exercises.
"""

from builders import (
    new_atoms, add_atom, bond, carbon_ring, benzene_ring,
    attach_chain, attach_atom, attach_branch, attach_hydroxymethyl,
    attach_chloromethyl, attach_carboxyl, attach_acyl, attach_nitrile,
    attach_vinyl, attach_nitro, attach_bridged_arm, ester_bridge,
    ring_with_chain_listed_first, payload,
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

def _isopropylcyclohexane():
    a = carbon_ring(6)
    attach_branch(a, 0, "isopropyl")
    return a

CASES["isopropylcyclohexane"] = exact(payload(_isopropylcyclohexane()), "(1-methylethyl)cyclohexane")

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

CASES["135_trimethylbenzene"] = exact(payload(_trimethylbenzene()), "1,3,5-trimethylbenzene")
CASES["cyclohexene_should_reject"] = reject(payload(carbon_ring(6, [2, 1, 1, 1, 1, 1])))

def _benzyl_alcohol():
    a = benzene_ring("A")
    attach_hydroxymethyl(a, 0)
    return a

CASES["phenylmethanol"] = exact(payload(_benzyl_alcohol()), "phenylmethanol")

def _styrene():
    a = benzene_ring("A")
    c6 = add_atom(a, "C")
    bond(a, 0, c6)
    c7 = add_atom(a, "C")
    bond(a, c6, c7, 2)
    return a

CASES["ethenylbenzene"] = exact(payload(_styrene()), "ethenylbenzene")

def _isopropylbenzene():
    a = benzene_ring("A")
    attach_branch(a, 0, "isopropyl")
    return a

CASES["isopropylbenzene"] = exact(payload(_isopropylbenzene()), "(1-methylethyl)benzene")

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

# ── Wave 2.b: G1 unified substituent namer (halo-tip / branched-alkyl subs) ──
def _tertbutylbenzene():
    a = benzene_ring("A")
    attach_branch(a, 0, "tertbutyl")
    return a

CASES["tertbutylbenzene"] = exact(payload(_tertbutylbenzene()), "(1,1-dimethylethyl)benzene")

def _hydroxymethylcyclohexane():
    a = carbon_ring(6)
    attach_hydroxymethyl(a, 0)
    return a

CASES["hydroxymethylcyclohexane"] = exact(payload(_hydroxymethylcyclohexane()), "(hydroxymethyl)cyclohexane")

def _chloromethylbenzene_mono():
    a = benzene_ring("A")
    attach_chloromethyl(a, 0)
    return a

CASES["chloromethylbenzene_mono"] = exact(payload(_chloromethylbenzene_mono()), "(chloromethyl)benzene")

def _chloro_chloromethyl_benzene():   # composite-name alphabetisation + paren dashes
    a = benzene_ring("A")
    attach_atom(a, 0, "Cl")
    attach_chloromethyl(a, 3)
    return a

CASES["1chloro4chloromethylbenzene"] = exact(
    payload(_chloro_chloromethyl_benzene()), "1-chloro-4-(chloromethyl)benzene")

def _chloromethyl_hydroxymethyl_cyclohexane():   # pins the )digit dash rule
    a = carbon_ring(6)
    attach_chloromethyl(a, 0)
    attach_hydroxymethyl(a, 1)
    return a

CASES["chloromethyl_hydroxymethyl_cyclohexane"] = exact(
    payload(_chloromethyl_hydroxymethyl_cyclohexane()),
    "1-(chloromethyl)-2-(hydroxymethyl)cyclohexane")

def _2_methylpropanal():              # Q4 fix: was "(methyl)" garbage territory
    a = _chain(3)
    attach_acyl(a, 1)
    return a

CASES["2_methylpropanal"] = exact(payload(_2_methylpropanal()), "2-methylpropan-1-al")

# ── Wave 2.b negatives: G1 guard gates ──
def _dichloromethylbenzene():         # -CHCl2: multi-group tip (Q2 gate)
    a = benzene_ring("A")
    c = add_atom(a, "C")
    bond(a, 0, c)
    attach_atom(a, c, "Cl")
    attach_atom(a, c, "Cl")
    return a

CASES["dichloromethylbenzene_should_reject"] = reject(payload(_dichloromethylbenzene()))

def _cyclohexylbenzene():             # cyclic substituent (tree guard)
    a = benzene_ring("A")
    ring2 = [add_atom(a, "C") for _ in range(6)]
    for i in range(6):
        bond(a, ring2[i], ring2[(i + 1) % 6])
    bond(a, 0, ring2[0])
    return a

CASES["cyclohexylbenzene_should_reject"] = reject(payload(_cyclohexylbenzene()))

def _vinylcyclohexane():              # unsaturated sub on a NON-aromatic ring
    a = carbon_ring(6)
    attach_vinyl(a, 0)
    return a

CASES["vinylcyclohexane_should_reject"] = reject(payload(_vinylcyclohexane()))

def _aminomethylbenzene():            # suffix-capable tip on benzene (D5 gate)
    a = benzene_ring("A")
    c = add_atom(a, "C")
    bond(a, 0, c)
    attach_atom(a, c, "N")
    return a

CASES["aminomethylbenzene_should_reject"] = reject(payload(_aminomethylbenzene()))

def _2_phenylethanol():               # -CH2CH2OH: suffix-capable tip, 2-carbon spine
    a = benzene_ring("A")
    chain = attach_chain(a, 0, 2)
    attach_atom(a, chain[-1], "O")
    return a

CASES["2_phenylethanol_should_reject"] = reject(payload(_2_phenylethanol()))

def _1_phenylethanol():               # OH on a NON-tip substituent carbon
    a = benzene_ring("A")
    c1 = add_atom(a, "C")
    bond(a, 0, c1)
    attach_chain(a, c1, 1)
    attach_atom(a, c1, "O")
    return a

CASES["1_phenylethanol_should_reject"] = reject(payload(_1_phenylethanol()))

def _cyanomethylbenzene():            # carbon-subsuming tip (Q3 gate)
    a = benzene_ring("A")
    c = add_atom(a, "C")
    bond(a, 0, c)
    attach_nitrile(a, c)
    return a

CASES["cyanomethylbenzene_should_reject"] = reject(payload(_cyanomethylbenzene()))

# ── Wave 2.c: retained aromatic parent names (G4/E1) ──
def _mono_benzene(attacher, *args, **kw):
    a = benzene_ring("A")
    attacher(a, 0, *args, **kw)
    return a

CASES["phenol"]           = exact(payload(_mono_benzene(attach_atom, "O")), "phenol")
CASES["phenylamine"]      = exact(payload(_mono_benzene(attach_atom, "N")), "phenylamine")
CASES["benzoic_acid"]     = exact(payload(_mono_benzene(attach_carboxyl)), "benzoic acid")
CASES["benzaldehyde"]     = exact(payload(_mono_benzene(attach_acyl)), "benzaldehyde")
CASES["benzamide"]        = exact(payload(_mono_benzene(attach_acyl, terminal="N")), "benzamide")
CASES["benzoyl_chloride"] = exact(payload(_mono_benzene(attach_acyl, terminal="Cl")), "benzoyl chloride")
CASES["benzonitrile"]     = exact(payload(_mono_benzene(attach_nitrile)), "benzonitrile")
CASES["phenylethanone"]   = exact(payload(_mono_benzene(attach_acyl, tail_len=1)), "phenylethanone")

def _phenol_plus(pos, attacher, *args):
    a = benzene_ring("A")
    attach_atom(a, 0, "O")
    attacher(a, pos, *args)
    return a

CASES["2_chlorophenol"]         = exact(payload(_phenol_plus(1, attach_atom, "Cl")), "2-chlorophenol")
CASES["2_chlorophenol_reverse"] = exact(payload(_phenol_plus(5, attach_atom, "Cl")), "2-chlorophenol")  # anchor direction choice
CASES["4_bromophenol"]          = exact(payload(_phenol_plus(3, attach_atom, "Br")), "4-bromophenol")
CASES["4_aminophenol"]          = exact(payload(_phenol_plus(3, attach_atom, "N")), "4-aminophenol")    # OH outranks NH2

def _2_hydroxybenzoic_acid():
    a = benzene_ring("A")
    attach_carboxyl(a, 0)
    attach_atom(a, 1, "O")
    return a

CASES["2_hydroxybenzoic_acid"] = exact(payload(_2_hydroxybenzoic_acid()), "2-hydroxybenzoic acid")  # COOH outranks OH; OH demotes

def _3_methylbenzaldehyde():
    a = benzene_ring("A")
    attach_acyl(a, 0)
    attach_chain(a, 2, 1)
    return a

CASES["3_methylbenzaldehyde"] = exact(payload(_3_methylbenzaldehyde()), "3-methylbenzaldehyde")

def _benzene_12_diol():
    a = benzene_ring("A")
    attach_atom(a, 0, "O")
    attach_atom(a, 1, "O")
    return a

CASES["benzene_12_diol_should_reject"] = reject(payload(_benzene_12_diol()))  # D7

def _benzene_14_dicarboxylic():
    a = benzene_ring("A")
    attach_carboxyl(a, 0)
    attach_carboxyl(a, 3)
    return a

CASES["benzene_dicarboxylic_should_reject"] = reject(payload(_benzene_14_dicarboxylic()))

def _chloro_benzyl_alcohol():         # chain-parent pattern + extra sub → gate
    a = benzene_ring("A")
    attach_hydroxymethyl(a, 0)
    attach_atom(a, 3, "Cl")
    return a

CASES["chloro_benzyl_alcohol_should_reject"] = reject(payload(_chloro_benzyl_alcohol()))

def _ethynylbenzene():                # C≡C is not the styrene pattern
    a = benzene_ring("A")
    c1 = add_atom(a, "C")
    bond(a, 0, c1)
    c2 = add_atom(a, "C")
    bond(a, c1, c2, 3)
    return a

CASES["ethynylbenzene_should_reject"] = reject(payload(_ethynylbenzene()))

# ── Wave 2.d: no ring-substituent-count limit (E3) — tie-break battery ──
def _polysub_benzene(spec_list):     # [(pos, attacher, args...), ...]
    a = benzene_ring("A")
    for pos, attacher, *args in spec_list:
        attacher(a, pos, *args)
    return a

CASES["4_chloro_12_dimethylbenzene"] = exact(
    payload(_polysub_benzene([(0, attach_chain, 1), (1, attach_chain, 1), (3, attach_atom, "Cl")])),
    "4-chloro-1,2-dimethylbenzene")               # lowest locant SET beats alphabetical
CASES["1_bromo_3_chloro_5_methylbenzene"] = exact(
    payload(_polysub_benzene([(0, attach_atom, "Br"), (2, attach_atom, "Cl"), (4, attach_chain, 1)])),
    "1-bromo-3-chloro-5-methylbenzene")           # set tie → first-alphabetical gets lowest
CASES["1245_tetramethylbenzene"] = exact(
    payload(_polysub_benzene([(0, attach_chain, 1), (1, attach_chain, 1),
                              (3, attach_chain, 1), (4, attach_chain, 1)])),
    "1,2,4,5-tetramethylbenzene")

def _246_tribromophenol():                        # end-state flagship (E1+E3)
    a = benzene_ring("A")
    attach_atom(a, 0, "O")
    for p in (1, 3, 5):
        attach_atom(a, p, "Br")
    return a

CASES["246_tribromophenol"] = exact(payload(_246_tribromophenol()), "2,4,6-tribromophenol")

def _24_dibromophenol():                          # anchored direction choice at n=3
    a = benzene_ring("A")
    attach_atom(a, 0, "O")
    attach_atom(a, 1, "Br")
    attach_atom(a, 3, "Br")
    return a

CASES["24_dibromophenol"] = exact(payload(_24_dibromophenol()), "2,4-dibromophenol")

def _124_trimethylcyclohexane():                  # cycloalkane n≥3 regression (never limited)
    a = carbon_ring(6)
    for p in (0, 1, 3):
        attach_chain(a, p, 1)
    return a

CASES["124_trimethylcyclohexane"] = exact(payload(_124_trimethylcyclohexane()), "1,2,4-trimethylcyclohexane")

# ── Phase 3: bridging heteroatoms (ethers, esters, secondary amines) ──
def _bridged(len_a, len_b, bridge_element):
    """Two straight chains joined end-to-end through a single bridging atom
    (R-O-R' ether, R-NH-R' secondary amine, R-S-R' sulfide)."""
    a = new_atoms()
    c_a = add_atom(a, "C")
    attach_chain(a, c_a, len_a - 1)
    c_b = add_atom(a, "C")
    attach_chain(a, c_b, len_b - 1)
    x = add_atom(a, bridge_element)
    bond(a, c_a, x)
    bond(a, c_b, x)
    return a

def _ester(acid_len, alkyl_len):
    """R-CO-O-R': carbonyl carbon is the first of `acid_len` acid-side carbons.
    acid_len=1 gives a methanoate (H-CO-O-R')."""
    a = new_atoms()
    acid_c = add_atom(a, "C")
    attach_chain(a, acid_c, acid_len - 1)
    attach_atom(a, acid_c, "O", order=2)
    alkyl_c = add_atom(a, "C")
    attach_chain(a, alkyl_c, alkyl_len - 1)
    o = add_atom(a, "O")
    bond(a, acid_c, o)
    bond(a, alkyl_c, o)
    return a

CASES["methoxymethane"] = exact(payload(_bridged(1, 1, "O")), "methoxymethane")
CASES["methoxyethane"] = exact(payload(_bridged(1, 2, "O")), "methoxyethane")
CASES["ethoxypropane"] = exact(payload(_bridged(2, 3, "O")), "ethoxypropane")
# side order must not matter: longer side listed first still names the same
CASES["methoxyethane_long_side_first"] = exact(payload(_bridged(2, 1, "O")), "methoxyethane")

CASES["methyl_ethanoate"] = exact(payload(_ester(2, 1)), "methyl ethanoate")
CASES["ethyl_ethanoate"] = exact(payload(_ester(2, 2)), "ethyl ethanoate")
CASES["methyl_propanoate"] = exact(payload(_ester(3, 1)), "methyl propanoate")
CASES["methyl_methanoate"] = exact(payload(_ester(1, 1)), "methyl methanoate")

# verifies the "N-" prefix stays uppercase (applied after FormatName's lowercasing)
CASES["n_methylmethanamine"] = exact(payload(_bridged(1, 1, "N")), "N-methylmethanamine")
CASES["n_methylethanamine"] = exact(payload(_bridged(1, 2, "N")), "N-methylethanamine")

# ── Phase 3 negative cases: must reject, not silently misname ──
def _anisole():
    a = benzene_ring("A")
    o = add_atom(a, "O")
    bond(a, 0, o)
    c = add_atom(a, "C")
    bond(a, o, c)
    return a

# Wave 4 (E6 stage a) flagship: the O-bridge-on-ring shape routes to the cyclic
# path and names as a substituted benzene, not the misnamed "methoxyhexane"
# (bridged exit) or "hydroxybenzene" (cyclic exit dropping the methyl) traps.
CASES["methoxybenzene"] = exact(payload(_anisole()), "methoxybenzene")

def _heteroatom_ring(carbon_count, hetero_element):
    """Ring closed THROUGH the heteroatom (THF, ethylene oxide, pyrrolidine).
    Invisible to IsCyclic() (carbon subgraph only) - must be caught by
    SplitAtBridgingAtom's containment guard, else THF names 'butoxybutane'."""
    a = new_atoms()
    c0 = add_atom(a, "C")
    chain = attach_chain(a, c0, carbon_count - 1)
    x = add_atom(a, hetero_element)
    bond(a, c0, x)
    bond(a, chain[-1] if chain else c0, x)
    return a

CASES["thf_should_reject"] = reject(payload(_heteroatom_ring(4, "O")))
CASES["ethylene_oxide_should_reject"] = reject(payload(_heteroatom_ring(2, "O")))
CASES["pyrrolidine_should_reject"] = reject(payload(_heteroatom_ring(4, "N")))

def _2_methoxypropane():
    a = _chain(3)
    o = add_atom(a, "O")
    bond(a, 1, o)  # bridge attached at propane's MIDDLE carbon
    c = add_atom(a, "C")
    bond(a, o, c)
    return a

CASES["2_methoxypropane_should_reject"] = reject(payload(_2_methoxypropane()))

def _2_methoxyethanol():
    a = _bridged(1, 2, "O")
    attach_atom(a, 2, "O")  # OH on the far carbon of the ethyl side
    return a

CASES["2_methoxyethanol_should_reject"] = reject(payload(_2_methoxyethanol()))

def _methyl_vinyl_ether():
    a = new_atoms()
    c0 = add_atom(a, "C")
    c1 = add_atom(a, "C")
    bond(a, c0, c1, 2)  # CH2=CH-
    c2 = add_atom(a, "C")
    o = add_atom(a, "O")
    bond(a, c1, o)
    bond(a, c2, o)
    return a

CASES["methyl_vinyl_ether_should_reject"] = reject(payload(_methyl_vinyl_ether()))

def _dimethoxymethane():
    a = _bridged(1, 1, "O")  # CH3-O-CH3
    o2 = add_atom(a, "O")
    bond(a, 1, o2)  # second bridge off the base carbon
    c = add_atom(a, "C")
    bond(a, o2, c)
    return a

CASES["dimethoxymethane_should_reject"] = reject(payload(_dimethoxymethane()))

def _sym_anhydride(acyl_len):
    """R-CO-O-CO-R symmetric acid anhydride."""
    a = new_atoms()
    c_a = add_atom(a, "C")
    attach_chain(a, c_a, acyl_len - 1)
    attach_atom(a, c_a, "O", order=2)
    c_b = add_atom(a, "C")
    attach_chain(a, c_b, acyl_len - 1)
    attach_atom(a, c_b, "O", order=2)
    o = add_atom(a, "O")
    bond(a, c_a, o)
    bond(a, c_b, o)
    return a

CASES["ethanoic_anhydride"] = exact(payload(_sym_anhydride(2)), "ethanoic anhydride")
CASES["propanoic_anhydride"] = exact(payload(_sym_anhydride(3)), "propanoic anhydride")

def _n_methylethanamide():
    a = _bridged(2, 1, "N")  # CH3CH2-NH-CH3, then make the ethyl side's first C a carbonyl
    attach_atom(a, 0, "O", order=2)
    return a

CASES["n_methylethanamide"] = exact(payload(_n_methylethanamide()), "N-methylethanamide")

def _n_ethylethanamide():
    a = _bridged(2, 2, "N")    # CH3CH2-NH-CH2CH3
    attach_atom(a, 0, "O", order=2)  # carbonyl on side-A's first C -> acid side = ethanoyl
    return a

CASES["n_ethylethanamide"] = exact(payload(_n_ethylethanamide()), "N-ethylethanamide")

def _n_methylpropanamide():
    a = _bridged(3, 1, "N")    # CH3CH2CH2-NH-CH3
    attach_atom(a, 0, "O", order=2)  # acid side = propanoyl (3 C), N-side = methyl
    return a

CASES["n_methylpropanamide"] = exact(payload(_n_methylpropanamide()), "N-methylpropanamide")

def _nn_dimethylethanamide():
    a = new_atoms()
    acid_c = add_atom(a, "C")
    attach_chain(a, acid_c, 1)
    attach_atom(a, acid_c, "O", order=2)
    n = add_atom(a, "N")
    bond(a, acid_c, n)
    bond(a, n, add_atom(a, "C"))      # N-methyl
    bond(a, n, add_atom(a, "C"))      # N-methyl
    return a

CASES["nn_dimethylethanamide"] = exact(payload(_nn_dimethylethanamide()), "N,N-dimethylethanamide")

def _n_ethyl_n_methylethanamide():
    a = new_atoms()
    acid_c = add_atom(a, "C")
    attach_chain(a, acid_c, 1)
    attach_atom(a, acid_c, "O", order=2)
    n = add_atom(a, "N")
    bond(a, acid_c, n)
    bond(a, n, add_atom(a, "C"))                                   # N-methyl
    et = add_atom(a, "C")
    bond(a, n, et)
    attach_chain(a, et, 1)                                         # N-ethyl
    return a

CASES["n_ethyl_n_methylethanamide"] = exact(
    payload(_n_ethyl_n_methylethanamide()), "N-ethyl-N-methylethanamide")

def _diacetimide():                    # imide: N with two acyl carbons
    a = new_atoms()
    c_a = add_atom(a, "C")
    attach_chain(a, c_a, 1)
    attach_atom(a, c_a, "O", order=2)
    c_b = add_atom(a, "C")
    attach_chain(a, c_b, 1)
    attach_atom(a, c_b, "O", order=2)
    n = add_atom(a, "N")
    bond(a, c_a, n)
    bond(a, c_b, n)
    return a

CASES["diacetimide_should_reject"] = reject(payload(_diacetimide()))

def _n_phenylethanamide():
    a = benzene_ring("A")
    n = add_atom(a, "N")
    bond(a, 0, n)
    c = add_atom(a, "C")
    bond(a, n, c)
    attach_chain(a, c, 1)
    attach_atom(a, c, "O", order=2)
    return a

# Wave 4 (E6 stage e) flagship — arrives early, as Phase C fallout: IsAcidSide
# (Phase A) + AlkylSideName's AromaticRing case (Phase C) already compose
# through NameSecondaryAmide with zero additional code (R7/§11).
CASES["n_phenylethanamide"] = exact(payload(_n_phenylethanamide()), "N-phenylethanamide")

def _trimethylamine():
    a = new_atoms()
    n = add_atom(a, "N")
    for _ in range(3):
        c = add_atom(a, "C")
        bond(a, n, c)
    return a

CASES["trimethylamine_should_reject"] = reject(payload(_trimethylamine()))
CASES["dimethyl_sulfide_should_reject"] = reject(payload(_bridged(1, 1, "S")))

# ── Wave 3 / G2: pattern strictness - everything near-nitro still rejects ──
def _charge_separated_nitro():
    a = new_atoms()
    c = add_atom(a, "C")
    n = add_atom(a, "N")
    bond(a, c, n)
    attach_atom(a, n, "O", order=2)
    attach_atom(a, n, "O", order=1)   # single-bonded O: NOT the accepted form
    return a

CASES["charge_separated_nitro_should_reject"] = reject(payload(_charge_separated_nitro()))

def _nitrosomethane():
    a = new_atoms()
    c = add_atom(a, "C")
    n = add_atom(a, "N")
    bond(a, c, n)
    attach_atom(a, n, "O", order=2)
    return a

CASES["nitroso_should_reject"] = reject(payload(_nitrosomethane()))

def _methyl_nitrite():                 # C-O-N=O
    a = new_atoms()
    c = add_atom(a, "C")
    o = add_atom(a, "O")
    bond(a, c, o)
    n = add_atom(a, "N")
    bond(a, o, n)
    attach_atom(a, n, "O", order=2)
    return a

CASES["methyl_nitrite_should_reject"] = reject(payload(_methyl_nitrite()))

def _methyl_nitrate():                 # C-O-NO2 (N has no C neighbour)
    a = new_atoms()
    c = add_atom(a, "C")
    o = add_atom(a, "O")
    bond(a, c, o)
    n = add_atom(a, "N")
    bond(a, o, n)
    attach_atom(a, n, "O", order=2)
    attach_atom(a, n, "O", order=2)
    return a

CASES["methyl_nitrate_should_reject"] = reject(payload(_methyl_nitrate()))

def _dimethyl_peroxide():              # heteroatom-heteroatom throw intact: O-O
    a = new_atoms()
    c1 = add_atom(a, "C")
    o1 = add_atom(a, "O")
    bond(a, c1, o1)
    o2 = add_atom(a, "O")
    bond(a, o1, o2)
    c2 = add_atom(a, "C")
    bond(a, o2, c2)
    return a

CASES["dimethyl_peroxide_should_reject"] = reject(payload(_dimethyl_peroxide()))

def _dimethylhydrazine():              # heteroatom-heteroatom throw intact: N-N
    a = new_atoms()
    c1 = add_atom(a, "C")
    n1 = add_atom(a, "N")
    bond(a, c1, n1)
    n2 = add_atom(a, "N")
    bond(a, n1, n2)
    c2 = add_atom(a, "C")
    bond(a, n2, c2)
    return a

CASES["nn_dimethylhydrazine_should_reject"] = reject(payload(_dimethylhydrazine()))

def _nitrogen_dioxide():               # NO2 with no carbon
    a = new_atoms()
    n = add_atom(a, "N")
    attach_atom(a, n, "O", order=2)
    attach_atom(a, n, "O", order=2)
    return a

CASES["nitrogen_dioxide_should_reject"] = reject(payload(_nitrogen_dioxide()))

def _pentavalent_carbon():             # ValidateValences: deliberate rejection now
    a = new_atoms()
    c = add_atom(a, "C")
    for _ in range(5):
        attach_chain(a, c, 1)
    return a

CASES["pentavalent_carbon_should_reject"] = reject(payload(_pentavalent_carbon()))

def _overbonded_nitro_oxygen():        # terminal-O check: O bridging onward
    a = new_atoms()
    c = add_atom(a, "C")
    n = add_atom(a, "N")
    bond(a, c, n)
    attach_atom(a, n, "O", order=2)
    o2 = add_atom(a, "O")
    bond(a, n, o2, 2)
    c2 = add_atom(a, "C")
    bond(a, o2, c2)                    # over-valent O - must NOT be consumed
    return a

CASES["overbonded_nitro_oxygen_should_reject"] = reject(payload(_overbonded_nitro_oxygen()))

# ── Wave 3 / E4: nitro positives ──
def _nitrobenzene(pattern="A"):
    a = benzene_ring(pattern)
    attach_nitro(a, 0)
    return a

CASES["nitrobenzene"] = exact(payload(_nitrobenzene()), "nitrobenzene")
CASES["nitrobenzene_kekule_b"] = exact(payload(_nitrobenzene("B")), "nitrobenzene")

def _nitrophenol(pos):
    a = benzene_ring()
    attach_atom(a, 0, "O")
    attach_nitro(a, pos)
    return a

CASES["2_nitrophenol"] = exact(payload(_nitrophenol(1)), "2-nitrophenol")
CASES["4_nitrophenol"] = exact(payload(_nitrophenol(3)), "4-nitrophenol")

def _246_trinitrophenol():
    a = benzene_ring()
    attach_atom(a, 0, "O")
    for p in (1, 3, 5):
        attach_nitro(a, p)
    return a

CASES["246_trinitrophenol"] = exact(payload(_246_trinitrophenol()), "2,4,6-trinitrophenol")

def _3_nitrobenzoic_acid():
    a = benzene_ring()
    attach_carboxyl(a, 0)
    attach_nitro(a, 2)
    return a

CASES["3_nitrobenzoic_acid"] = exact(payload(_3_nitrobenzoic_acid()), "3-nitrobenzoic acid")

def _3_nitrobenzaldehyde():
    a = benzene_ring()
    attach_acyl(a, 0)                  # -CHO
    attach_nitro(a, 2)
    return a

CASES["3_nitrobenzaldehyde"] = exact(payload(_3_nitrobenzaldehyde()), "3-nitrobenzaldehyde")

def _4_nitrophenylamine():             # real NH2 = parent; nitro N = prefix (the
    a = benzene_ring()                 # degree-aware discrimination proof)
    attach_atom(a, 0, "N")
    attach_nitro(a, 3)
    return a

CASES["4_nitrophenylamine"] = exact(payload(_4_nitrophenylamine()), "4-nitrophenylamine")

def _dinitrobenzene(pos):
    a = benzene_ring()
    attach_nitro(a, 0)
    attach_nitro(a, pos)
    return a

CASES["12_dinitrobenzene"] = exact(payload(_dinitrobenzene(1)), "1,2-dinitrobenzene")
CASES["14_dinitrobenzene"] = exact(payload(_dinitrobenzene(3)), "1,4-dinitrobenzene")

def _trinitrotoluene():
    a = benzene_ring()
    attach_chain(a, 0, 1)
    for p in (1, 3, 5):
        attach_nitro(a, p)
    return a

# Systematic form, NOT "2,4,6-trinitromethylbenzene" - design ruling R6.
# Lowest-locants-as-a-set picks {1,2,3,5} over {1,2,4,6}, so methyl lands at
# 2 rather than 1 - matching TNT's actual systematic name, not the plan's
# predicted "1-methyl-2,4,6-trinitrobenzene" (a locant tie-break error in the
# design doc, verified against real IUPAC nomenclature on implementation).
CASES["trinitrotoluene_systematic"] = exact(
    payload(_trinitrotoluene()), "2-methyl-1,3,5-trinitrobenzene")

def _nitrocyclohexane():
    a = carbon_ring(6)
    attach_nitro(a, 0)
    return a

CASES["nitrocyclohexane"] = exact(payload(_nitrocyclohexane()), "nitrocyclohexane")

def _1_methyl_2_nitrocyclohexane():
    a = carbon_ring(6)
    attach_chain(a, 0, 1)
    attach_nitro(a, 1)
    return a

CASES["1_methyl_2_nitrocyclohexane"] = exact(
    payload(_1_methyl_2_nitrocyclohexane()), "1-methyl-2-nitrocyclohexane")

def _2_nitropropane():
    a = new_atoms()
    c0 = add_atom(a, "C")
    c1 = add_atom(a, "C")
    bond(a, c0, c1)
    c2 = add_atom(a, "C")
    bond(a, c1, c2)
    attach_nitro(a, c1)
    return a

CASES["2_nitropropane"] = exact(payload(_2_nitropropane()), "2-nitropropane")

def _22_dinitropropane():
    a = new_atoms()
    c0 = add_atom(a, "C")
    c1 = add_atom(a, "C")
    bond(a, c0, c1)
    c2 = add_atom(a, "C")
    bond(a, c1, c2)
    attach_nitro(a, c1)
    attach_nitro(a, c1)
    return a

CASES["22_dinitropropane"] = exact(payload(_22_dinitropropane()), "2,2-dinitropropane")

def _nitromethane():
    a = new_atoms()
    c = add_atom(a, "C")
    attach_nitro(a, c)
    return a

# House redundant-locant style (cf. 2-methylpropan-1-al).
CASES["nitromethane"] = exact(payload(_nitromethane()), "1-nitromethane")

# Hypervalent-N FillImplicitHydrogens pin-down (master-plan requirement): no
# client-side H fill - the API's own filler must add 3 H to C and NONE to N.
def _nitromethane_implicit_h():
    a = new_atoms()
    c = add_atom(a, "C")
    attach_nitro(a, c)
    return {"atoms": a, "specificationSet": "AllGroups"}   # note: no payload()/fill

CASES["nitromethane_implicit_h"] = exact(_nitromethane_implicit_h(), "1-nitromethane")

def _3_nitropropanoic_acid():          # merge-integrity proof: COOH merge with a
    a = new_atoms()                    # PolyatomicGroup elsewhere in the molecule
    c0 = add_atom(a, "C")
    c1 = add_atom(a, "C")
    bond(a, c0, c1)
    c2 = add_atom(a, "C")
    bond(a, c1, c2)
    attach_atom(a, c0, "O", order=2)
    attach_atom(a, c0, "O")
    attach_nitro(a, c2)
    return a

CASES["3_nitropropanoic_acid"] = exact(
    payload(_3_nitropropanoic_acid()), "3-nitropropan-1-oic acid")

# ── Wave 3 / E4: nitro negatives (capability gates that must hold with naming enabled) ──
CASES["nitrobenzene_hydrocarbons_spec_should_reject"] = reject(
    payload(_nitrobenzene(), spec="Hydrocarbons"))         # spec gate, ring path

CASES["2_nitropropane_hydrocarbons_spec_should_reject"] = reject(
    payload(_2_nitropropane(), spec="Hydrocarbons"))       # spec gate, chain path (CheckGroups)

def _nitromethylbenzene():             # ring-CH2-NO2 composite substituent
    a = benzene_ring()
    c = add_atom(a, "C")
    bond(a, 0, c)
    attach_nitro(a, c)
    return a

CASES["nitromethylbenzene_should_reject"] = reject(payload(_nitromethylbenzene()))

def _methyl_nitroethanoate():          # nitro on a bridged (ester) side
    a = new_atoms()
    c_acid = add_atom(a, "C")
    attach_atom(a, c_acid, "O", order=2)
    o_bridge = add_atom(a, "O")
    bond(a, c_acid, o_bridge)
    c_methyl = add_atom(a, "C")
    bond(a, o_bridge, c_methyl)
    c_alpha = add_atom(a, "C")
    bond(a, c_acid, c_alpha)
    attach_nitro(a, c_alpha)
    return a

CASES["methyl_nitroethanoate_should_reject"] = reject(payload(_methyl_nitroethanoate()))

# ── Wave 4 / E6 stage (a): aromatic and ring ethers (alkoxy substituents) ──
def _ethoxybenzene():
    a = benzene_ring("A")
    attach_bridged_arm(a, 0, "O", 2)
    return a

CASES["ethoxybenzene"] = exact(payload(_ethoxybenzene()), "ethoxybenzene")

def _methoxycyclohexane():
    a = carbon_ring(6)
    attach_bridged_arm(a, 0, "O", 1)
    return a

CASES["methoxycyclohexane"] = exact(payload(_methoxycyclohexane()), "methoxycyclohexane")

def _1_methoxy_2_methylcyclohexane():
    a = carbon_ring(6)
    attach_bridged_arm(a, 0, "O", 1)
    attach_chain(a, 1, 1)
    return a

CASES["1_methoxy_2_methylcyclohexane"] = exact(
    payload(_1_methoxy_2_methylcyclohexane()), "1-methoxy-2-methylcyclohexane")  # alphabetical: methoxy < methyl

def _2_methoxyphenol():
    a = benzene_ring("A")
    attach_atom(a, 0, "O")             # phenol OH — the retained-parent anchor
    attach_bridged_arm(a, 1, "O", 1)
    return a

CASES["2_methoxyphenol"] = exact(payload(_2_methoxyphenol()), "2-methoxyphenol")

def _4_methoxybenzaldehyde():
    a = benzene_ring("A")
    attach_acyl(a, 0)                  # -CHO — the retained carbon-parent anchor
    attach_bridged_arm(a, 3, "O", 1)
    return a

CASES["4_methoxybenzaldehyde"] = exact(payload(_4_methoxybenzaldehyde()), "4-methoxybenzaldehyde")

def _1_methoxy_4_nitrobenzene():
    a = benzene_ring("A")
    attach_bridged_arm(a, 0, "O", 1)
    attach_nitro(a, 3)
    return a

CASES["1_methoxy_4_nitrobenzene"] = exact(payload(_1_methoxy_4_nitrobenzene()), "1-methoxy-4-nitrobenzene")

# ── Wave 4 / E6 stage (a) negatives: ether shapes that must stay rejected ──
def _diphenyl_ether():
    a = benzene_ring("A")
    o = add_atom(a, "O")
    bond(a, 0, o)
    ring2_offset = len(a)
    ring2 = benzene_ring("A")
    for atom in ring2:
        atom["bonds"] = [{"to": b["to"] + ring2_offset, "order": b["order"]} for b in atom["bonds"]]
    a.extend(ring2)
    bond(a, o, ring2_offset)
    return a

CASES["diphenyl_ether_should_reject"] = reject(payload(_diphenyl_ether()))

def _benzyl_methyl_ether():             # Ph-CH2-O-CH3 — O not on a ring carbon
    a = benzene_ring("A")
    ch2 = attach_chain(a, 0, 1)[0]
    attach_bridged_arm(a, ch2, "O", 1)
    return a

CASES["benzyl_methyl_ether_should_reject"] = reject(payload(_benzyl_methyl_ether()))

def _dihydrobenzofuran():
    # 2,3-dihydrobenzofuran: benzene C1-C6; O bonded to C1 and CH2a; CH2a-CH2b;
    # CH2b bonded to C2. The R5 fused-guard pin (§8): without the
    # adjacent-to-blocked check, this arm looks like a clean unbranched
    # saturated "ethoxy" and would silently amputate the fusion.
    a = benzene_ring("A")
    o = add_atom(a, "O")
    bond(a, 0, o)
    ch2a = add_atom(a, "C")
    bond(a, o, ch2a)
    ch2b = add_atom(a, "C")
    bond(a, ch2a, ch2b)
    bond(a, ch2b, 1)                    # closes the fused ring back onto C2
    return a

CASES["dihydrobenzofuran_should_reject"] = reject(payload(_dihydrobenzofuran()))

def _phenyl_vinyl_ether():              # unsaturated arm
    a = benzene_ring("A")
    o = add_atom(a, "O")
    bond(a, 0, o)
    c1 = add_atom(a, "C")
    bond(a, o, c1)
    c2 = add_atom(a, "C")
    bond(a, c1, c2, 2)
    return a

CASES["phenyl_vinyl_ether_should_reject"] = reject(payload(_phenyl_vinyl_ether()))

def _isopropoxybenzene():               # branched arm
    a = benzene_ring("A")
    o = add_atom(a, "O")
    bond(a, 0, o)
    attach_branch(a, o, "isopropyl")
    return a

CASES["isopropoxybenzene_should_reject"] = reject(payload(_isopropoxybenzene()))

def _2_chloroethoxybenzene():           # heteroatom in the arm
    a = benzene_ring("A")
    _, arm = attach_bridged_arm(a, 0, "O", 2)
    attach_atom(a, arm[-1], "Cl")
    return a

CASES["2_chloroethoxybenzene_should_reject"] = reject(payload(_2_chloroethoxybenzene()))

def _12_dimethoxybenzene():             # multi-bridge (R12) — message changes, status doesn't
    a = benzene_ring("A")
    attach_bridged_arm(a, 0, "O", 1)
    attach_bridged_arm(a, 1, "O", 1)
    return a

CASES["12_dimethoxybenzene_should_reject"] = reject(payload(_12_dimethoxybenzene()))

# ── Wave 4 / E6 stages (b)+(c): esters with ring sides, unsubstituted ──
def _methyl_benzoate():
    a = benzene_ring("A")
    methyl_c = add_atom(a, "C")
    ester_bridge(a, 0, methyl_c)
    return a

CASES["methyl_benzoate"] = exact(payload(_methyl_benzoate()), "methyl benzoate")

def _ethyl_benzoate():
    a = benzene_ring("A")
    ethyl_c = add_atom(a, "C")
    attach_chain(a, ethyl_c, 1)
    ester_bridge(a, 0, ethyl_c)
    return a

CASES["ethyl_benzoate"] = exact(payload(_ethyl_benzoate()), "ethyl benzoate")

def _phenyl_ethanoate():
    a = benzene_ring("A")
    acyl_c, bridge_o = ester_bridge(a, None, 0)
    attach_chain(a, acyl_c, 1)          # ethanoyl: acyl_c + 1 tail carbon
    return a

CASES["phenyl_ethanoate"] = exact(payload(_phenyl_ethanoate()), "phenyl ethanoate")

def _phenyl_propanoate():
    a = benzene_ring("A")
    acyl_c, bridge_o = ester_bridge(a, None, 0)
    attach_chain(a, acyl_c, 2)          # propanoyl: acyl_c + 2 tail carbons
    return a

CASES["phenyl_propanoate"] = exact(payload(_phenyl_propanoate()), "phenyl propanoate")

def _phenyl_methanoate():
    a = benzene_ring("A")
    ester_bridge(a, None, 0)            # bare acyl C -> methanoate (R7)
    return a

CASES["phenyl_methanoate"] = exact(payload(_phenyl_methanoate()), "phenyl methanoate")

def _phenyl_benzoate():                 # both sides aromatic (R7)
    a = benzene_ring("A")               # acid-side ring
    ring2_offset = len(a)
    ring2 = benzene_ring("A")           # alkyl-side ring
    for atom in ring2:
        atom["bonds"] = [{"to": b["to"] + ring2_offset, "order": b["order"]} for b in atom["bonds"]]
    a.extend(ring2)
    ester_bridge(a, 0, ring2_offset)
    return a

CASES["phenyl_benzoate"] = exact(payload(_phenyl_benzoate()), "phenyl benzoate")

# ── Wave 4 / E6 stages (b)+(c) negatives ──
def _benzoic_anhydride():               # E6xE7: aromatic acid sides stay unsupported
    a = benzene_ring("A")
    c1 = add_atom(a, "C")
    bond(a, 0, c1)
    attach_atom(a, c1, "O", order=2)
    bridge_o = add_atom(a, "O")
    bond(a, c1, bridge_o)
    c2 = add_atom(a, "C")
    bond(a, bridge_o, c2)
    attach_atom(a, c2, "O", order=2)
    ring2_offset = len(a)
    ring2 = benzene_ring("A")
    for atom in ring2:
        atom["bonds"] = [{"to": b["to"] + ring2_offset, "order": b["order"]} for b in atom["bonds"]]
    a.extend(ring2)
    bond(a, c2, ring2_offset)
    return a

CASES["benzoic_anhydride_should_reject"] = reject(payload(_benzoic_anhydride()))

def _benzyl_ethanoate():                # ring deeper in the alkyl side
    a = benzene_ring("A")
    ch2 = attach_chain(a, 0, 1)[0]
    acyl_c, bridge_o = ester_bridge(a, None, ch2)
    attach_chain(a, acyl_c, 1)
    return a

CASES["benzyl_ethanoate_should_reject"] = reject(payload(_benzyl_ethanoate()))

def _cyclohexyl_ethanoate():            # non-aromatic ring alkyl side
    a = carbon_ring(6)
    acyl_c, bridge_o = ester_bridge(a, None, 0)
    attach_chain(a, acyl_c, 1)
    return a

CASES["cyclohexyl_ethanoate_should_reject"] = reject(payload(_cyclohexyl_ethanoate()))

def _4_nitrophenyl_ethanoate():         # substituted alkyl-side ring (R4 gate)
    a = benzene_ring("A")
    acyl_c, bridge_o = ester_bridge(a, None, 0)
    attach_chain(a, acyl_c, 1)
    attach_nitro(a, 3)
    return a

CASES["4_nitrophenyl_ethanoate_should_reject"] = reject(payload(_4_nitrophenyl_ethanoate()))

def _aspirin():                         # 2-(ethanoyloxy)benzoic acid — recorded non-goal (§19)
    a = benzene_ring("A")
    attach_carboxyl(a, 0)
    acyl_c, bridge_o = ester_bridge(a, None, 1)
    attach_chain(a, acyl_c, 1)
    return a

CASES["aspirin_should_reject"] = reject(payload(_aspirin()))

# ── Wave 4 / E6 stage (d): aromatic secondary amines ──
def _n_methylphenylamine():
    a = benzene_ring("A")
    n = add_atom(a, "N")
    bond(a, 0, n)
    methyl = add_atom(a, "C")
    bond(a, n, methyl)
    return a

# The master plan's trap test: assert exactly "N-methylphenylamine", never the
# misnamed "phenylamine" (which would silently drop the N-methyl).
CASES["n_methylphenylamine"] = exact(payload(_n_methylphenylamine()), "N-methylphenylamine")

def _n_ethylphenylamine():
    a = benzene_ring("A")
    n = add_atom(a, "N")
    bond(a, 0, n)
    ethyl_c = add_atom(a, "C")
    bond(a, n, ethyl_c)
    attach_chain(a, ethyl_c, 1)
    return a

CASES["n_ethylphenylamine"] = exact(payload(_n_ethylphenylamine()), "N-ethylphenylamine")

def _n_butylphenylamine():
    a = benzene_ring("A")
    n = add_atom(a, "N")
    bond(a, 0, n)
    butyl_c = add_atom(a, "C")
    bond(a, n, butyl_c)
    attach_chain(a, butyl_c, 3)
    return a

# Forced-base rule (R6/§10): the ring wins as the base even though the butyl
# side has more carbons than the ring — never Min/Max length selection.
CASES["n_butylphenylamine"] = exact(payload(_n_butylphenylamine()), "N-butylphenylamine")

# ── Wave 4 / E6 stage (d) negatives ──
def _diphenylamine():
    a = benzene_ring("A")
    n = add_atom(a, "N")
    bond(a, 0, n)
    ring2_offset = len(a)
    ring2 = benzene_ring("A")
    for atom in ring2:
        atom["bonds"] = [{"to": b["to"] + ring2_offset, "order": b["order"]} for b in atom["bonds"]]
    a.extend(ring2)
    bond(a, n, ring2_offset)
    return a

CASES["diphenylamine_should_reject"] = reject(payload(_diphenylamine()))

def _4_methyl_n_methylphenylamine():    # substituted ring — R4's side.Count==6 gate
    a = benzene_ring("A")
    n = add_atom(a, "N")
    bond(a, 0, n)
    methyl_n = add_atom(a, "C")
    bond(a, n, methyl_n)
    attach_chain(a, 3, 1)                # ring methyl substituent, para position
    return a

CASES["4_methyl_n_methylphenylamine_should_reject"] = reject(payload(_4_methyl_n_methylphenylamine()))

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
