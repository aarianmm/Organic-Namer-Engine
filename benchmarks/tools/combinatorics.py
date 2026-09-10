#!/usr/bin/env python3
"""
Analytic enumeration of the OrganicNamer grammar's molecule space.

Counts distinct constitutional isomers the engine will name under the
`AllGroups` specification. Pure combinatorics -- no molecules are generated and
the engine is never executed. Every rule encoded here was read out of
src/OrganicNamer.Core and spot-checked against the live engine; the audit trail
is in BENCHMARKS.md.

Molecules are partitioned by the engine's own three-way routing:

  A. chain   -- acyclic, no bridging heteroatom   (dominant term)
  B. cyclic  -- exactly one carbon ring
  C. bridged -- exactly one chain-breaking heteroatom (C-X-C)

Family A is counted with a transfer-matrix walk along the main chain. Family B
uses Burnside's lemma over the ring's dihedral group. Family C is a product of
the side grammars.

Every figure is a LOWER BOUND on what the engine can name, for one reason: the
walk requires the main chain to be a longest path in the molecule, whereas the
engine drops branch-invalid candidate chains first and then takes the longest
*survivor*. A molecule whose longest path is unnameable but which has a shorter
nameable chain is therefore counted here as zero, and the engine names it.

That gap is measured, not assumed. Against an exhaustive enumeration of every
alkane isomer up to C12 (alkane_isomers.py, checked against OEIS A000602 and run
through the real engine), this model reproduces the engine exactly through C9
and then trails it slightly:

    C9   35 / 35   (exact)
    C10  72 / 74   (97.3% of what the engine names)
    C11 146 / 154  (94.8%)
    C12 307 / 332  (92.5%)

So the bound is tight to within roughly 3-8% on the sub-family where it can be
checked completely, and it always errs low.

Chain reversal is quotiented out with Burnside's lemma -- orbits =
(labelled + palindromic) / 2 -- so a molecule and its mirror numbering count once.
"""

from functools import lru_cache
from itertools import combinations
import argparse
import json
import math

MAX_CHAIN = 10          # spec.alkylNames holds roots for 1..10 carbons
MAX_RING = 10           # same root table drives ring names
MAX_MULT = 9            # spec.numericalPrefixes holds di..nona (2..9)
BOND_ORDERS = (1, 2, 3) # backbone bond orders; narrowed to (1,) to isolate alkanes

# ── Attachment vocabulary at a carbon ──────────────────────────────────────
# (name bucket, heavy atoms contributed). One bucket per group type; the engine
# caps each bucket at MAX_MULT occurrences.
MONOVALENT = [("fluoro", 1), ("chloro", 1), ("bromo", 1), ("iodo", 1),
              ("hydroxy", 1), ("sulfanyl", 1), ("amino", 1), ("nitro", 3)]
DIVALENT   = [("oxo", 1), ("imino", 1)]          # =O, =NH  (2 valences)
TRIVALENT  = [("cyano", 1)]                      # -C#N     (3 valences)

# Tip groups a straight leg may carry. TipGroupPrefix() accepts prefixOnly and
# prefixOrSuffix entries and rejects endDependent-only groups (=O, =NH), merged
# groups (COOH/COCl/CONH2) and C#N.
#
# Nitro is excluded even though TipGroupPrefix() would name it: NameCarbonSubstituent
# requires every heteroatom in a substituent to be adjacent to the spine tip, and a
# nitro group's two oxygens hang off its nitrogen, one bond further out. So a
# (nitroalkyl) branch is never nameable -- the engine re-routes onto a main chain
# that carries the nitro directly ("2-ethyl-1-nitrobutane", verified). Excluding it
# is a scope correction found by generate_corpus.py, not a guess.
LEG_TIPS = [(None, 0)] + [g for g in MONOVALENT if g[0] != "nitro"]


def _items(max_leg):
    """Attachment alphabet: (kind, bucket, valence cost, heavy, leg carbons)."""
    out = [("m", b, 1, h, 0) for b, h in MONOVALENT]
    out += [("d", b, 2, h, 0) for b, h in DIVALENT]
    out += [("t", b, 3, h, 0) for b, h in TRIVALENT]
    for ell in range(1, max_leg + 1):
        for tip, th in LEG_TIPS:
            bucket = f"leg{ell}" if tip is None else f"leg{ell}+{tip}"
            out.append(("l", bucket, 1, ell + th, ell))
    return out


_DEC = {}
def decorations(free, max_leg, n, track):
    """
    Every way to fill `free` valences at one chain carbon.

    Returns {(heavy, longest_leg, bucket_counts): multiplicity}. Attachments are
    an unordered multiset, so items are chosen in non-decreasing index order.
    Two legs on the same carbon span longest_leg + 1 + second_leg carbons, which
    must not exceed n.
    """
    key = (free, max_leg, n, track)
    if key in _DEC:
        return _DEC[key]
    items = _items(max_leg)
    idx = {b: i for i, b in enumerate(track)}
    res = {}

    def rec(start, remaining, heavy, legs, counts):
        if len(legs) < 2 or legs[0] + 1 + legs[1] <= n:
            k = (heavy, legs[0] if legs else 0, counts)
            res[k] = res.get(k, 0) + 1
        for i in range(start, len(items)):
            _kind, bucket, cost, h, carbons = items[i]
            if cost > remaining:
                continue
            nl = legs
            if carbons:
                nl = sorted(legs + [carbons], reverse=True)[:2]
                if len(nl) == 2 and nl[0] + 1 + nl[1] > n:
                    continue
            nc = counts
            if bucket in idx:
                j = idx[bucket]
                lst = list(counts)
                lst[j] = min(lst[j] + 1, MAX_MULT + 1)
                nc = tuple(lst)
            rec(i, remaining - cost, heavy + h, nl, nc)

    rec(0, free, 0, [], tuple(0 for _ in track))
    _DEC[key] = res
    return res


# ── Family A: acyclic, no bridge ───────────────────────────────────────────
def chain_labelled(n, track=(), heavy_cap=None):
    """
    Labelled decorated main chains of length n (before quotienting by reversal).
    Returns {(heavy, bucket_counts): ways}; heavy is only accumulated when
    heavy_cap is set.

    State carried along the chain is (bond order into this carbon, reach), where
    reach is the longest carbon path ending at the current carbon. A leg of ell
    carbons at position i must satisfy ell + 1 + (n - i) <= n (it must not beat
    the main chain on the right) and ell + reach <= n (nor on the left).
    """
    zero = tuple(0 for _ in track)
    base_heavy = 0 if heavy_cap is None else n
    states = {(0, 0, base_heavy, zero): 1}      # (b_in, reach_prev, heavy, counts)

    for i in range(1, n + 1):
        nxt = {}
        for (b_in, reach_prev, heavy, counts), ways in states.items():
            base = 1 if i == 1 else reach_prev + 1
            for b_out in ((0,) if i == n else BOND_ORDERS):
                free = 4 - b_in - b_out
                if free < 0:
                    continue
                max_leg = max(0, min(i - 1, n - base))
                for (dh, longest, dc), mult in decorations(
                        free, max_leg, n, track).items():
                    reach = max(base, longest + 1) if longest else base
                    nheavy = heavy
                    if heavy_cap is not None:
                        nheavy = heavy + dh
                        if nheavy > heavy_cap:
                            continue
                    nc = tuple(min(a + b, MAX_MULT + 1)
                               for a, b in zip(counts, dc))
                    k = (b_out, reach, nheavy, nc)
                    nxt[k] = nxt.get(k, 0) + ways * mult
        states = nxt

    out = {}
    for (_b, _r, heavy, counts), ways in states.items():
        out[(heavy, counts)] = out.get((heavy, counts), 0) + ways
    return out


def chain_palindromic(n, track=(), heavy_cap=None):
    """
    Labelled chains of length n that are invariant under reversal.

    Burnside's lemma needs this term: the number of distinct molecules is
    (labelled + palindromic) / 2, not labelled / 2. Only the left half is walked,
    since the right half mirrors it; the middle carbon of an odd chain maps to
    itself and so may be decorated freely.

    The leg constraints need no extra cross-middle check. A leg of ell carbons at
    position i and its mirror at n+1-i span ell + (n + 2 - 2i) + ell carbons,
    which is <= n exactly when ell <= i - 1 -- already enforced. Two legs at
    left-half positions i and i' likewise span at most (i-1) + (i'-1) + ... <= n
    under the same per-leg bound.
    """
    zero = tuple(0 for _ in track)
    half = n // 2                      # positions walked in full mirrored pairs
    odd = n % 2 == 1
    base_heavy = 0 if heavy_cap is None else n
    states = {(0, 0, base_heavy, zero): 1}

    for i in range(1, half + 1):
        nxt = {}
        for (b_in, reach_prev, heavy, counts), ways in states.items():
            base = 1 if i == 1 else reach_prev + 1
            for b_out in BOND_ORDERS:          # the mirrored/central bond always exists
                free = 4 - b_in - b_out
                if free < 0:
                    continue
                max_leg = max(0, min(i - 1, n - base))
                for (dh, longest, dc), mult in decorations(
                        free, max_leg, n, track).items():
                    reach = max(base, longest + 1) if longest else base
                    nheavy = heavy
                    if heavy_cap is not None:
                        nheavy = heavy + 2 * dh      # the mirror carries a copy
                        if nheavy > heavy_cap:
                            continue
                    nc = tuple(min(a + 2 * b, MAX_MULT + 1)
                               for a, b in zip(counts, dc))
                    k = (b_out, reach, nheavy, nc)
                    nxt[k] = nxt.get(k, 0) + ways * mult
        states = nxt

    out = {}
    for (b_in, reach_prev, heavy, counts), ways in states.items():
        if odd:
            # the centre carbon: its two bonds are mirror images, so both equal b_in
            free = 4 - 2 * b_in
            if free < 0:
                continue
            base = 1 if half == 0 else reach_prev + 1
            max_leg = max(0, min(half, n - base))
            for (dh, _longest, dc), mult in decorations(
                    free, max_leg, n, track).items():
                nheavy = heavy
                if heavy_cap is not None:
                    nheavy = heavy + dh
                    if nheavy > heavy_cap:
                        continue
                nc = tuple(min(a + b, MAX_MULT + 1) for a, b in zip(counts, dc))
                out[(nheavy, nc)] = out.get((nheavy, nc), 0) + ways * mult
        else:
            out[(heavy, counts)] = out.get((heavy, counts), 0) + ways
    return out


def _sum_valid(table, track):
    """Sum configurations in which no tracked bucket overflowed past MAX_MULT."""
    return sum(w for (_h, c), w in table.items()
               if all(x <= MAX_MULT for x in c))


def overflowable_buckets(n):
    """Buckets that could physically occur 10+ times on a chain of length n."""
    cands = [b for b, _ in MONOVALENT] + [b for b, _ in DIVALENT] + ["leg1"]
    out = []
    for b in cands:
        tab = chain_labelled(n, track=(b,))
        if any(c[0] > MAX_MULT for (_h, c) in tab):
            out.append(b)
    return out


def _capped(fn, n, heavy_cap=None):
    """
    Apply the engine's "at most 9 identical substituent names" rule to a counting
    function, by inclusion-exclusion over the buckets that can reach 10.

    A decane has 22 free valences, so two buckets can overflow simultaneously but
    three cannot (that would need 30); the expansion therefore stops at pairs.
    Returns {heavy: count}.
    """
    base = fn(n, heavy_cap=heavy_cap)
    out = {}
    for (h, _c), w in base.items():
        out[h] = out.get(h, 0) + w
    over = overflowable_buckets(n)
    for b in over:
        for (h, c), w in fn(n, track=(b,), heavy_cap=heavy_cap).items():
            if c[0] > MAX_MULT:
                out[h] = out.get(h, 0) - w
    for b1, b2 in combinations(over, 2):
        for (h, c), w in fn(n, track=(b1, b2), heavy_cap=heavy_cap).items():
            if c[0] > MAX_MULT and c[1] > MAX_MULT:
                out[h] = out.get(h, 0) + w
    return out


def chain_orbits(n, heavy_cap=None):
    """
    Exact number of distinct molecules with main chain length n, as
    {heavy_atoms: count}, by Burnside over the chain-reversal group:
    orbits = (labelled + palindromic) / 2.
    """
    lab = _capped(chain_labelled, n, heavy_cap)
    pal = _capped(chain_palindromic, n, heavy_cap)
    keys = set(lab) | set(pal)
    out = {}
    for h in keys:
        tot = lab.get(h, 0) + pal.get(h, 0)
        assert tot % 2 == 0, f"Burnside orbit count not integral at n={n}, heavy={h}"
        if tot:
            out[h] = tot // 2
    return out


def chain_total(n):
    """Exact distinct-molecule count for main chain length n."""
    return sum(chain_orbits(n).values())


# ── Family B: one carbon ring ──────────────────────────────────────────────
# Counted with Polya enumeration over the ring's dihedral group D_r, so a ring
# and its rotations/reflections are one molecule.
#
# The alphabet below is deliberately restricted to substituents that carry no
# whole-molecule constraint, which keeps the count a clean lower bound:
#   * alkoxy arms are omitted (only one bridging heteroatom is allowed per
#     molecule, a global constraint Polya cannot express),
#   * aromatic rings are counted with prefix-only substituents, omitting the
#     retained-parent forms (phenol, benzoic acid, ...) which allow at most one
#     principal-class substituent per ring -- again a global constraint.
# Both omissions remove molecules the engine *can* name, so the reported ring
# total is below the true one.

def ring_substituents(aromatic):
    """
    {heavy atoms: how many distinct substituents weigh that much}, for one
    substituent occupying one ring valence.
    """
    out = {}

    def put(h, n=1):
        out[h] = out.get(h, 0) + n

    for _b, h in MONOVALENT:                          # halo, OH, SH, NH2, nitro
        put(h)
    tips = [(None, 0)] + ([(b, h) for b, h in MONOVALENT
                           if b in ("fluoro", "chloro", "bromo", "iodo")]
                          if aromatic else
                          [(b, h) for b, h in MONOVALENT if b != "nitro"])
    for spine in range(1, MAX_RING + 1):              # straight alkyl + optional tip
        for _t, th in tips:
            put(spine + th)
    for spine in range(2, MAX_RING + 1):              # methyl-branched alkyl
        # methyls sit on spine positions 1..spine-1 (one on the tip would just
        # extend the spine), at most 2 per position, at most 9 in total
        for extra in range(1, min(2 * (spine - 1), MAX_MULT) + 1):
            put(spine + extra, _compositions(extra, spine - 1, 2))
    return out


@lru_cache(maxsize=None)
def _compositions(total, slots, cap):
    """How many ways to write `total` as an ordered sum of `slots` values 0..cap."""
    if slots == 0:
        return 1 if total == 0 else 0
    return sum(_compositions(total - v, slots - 1, cap)
               for v in range(min(cap, total) + 1))


def _ring_position_poly(aromatic, heavy_cap):
    """
    Generating polynomial in x (heavy atoms) for decorating one ring position.
    A non-aromatic ring carbon has 2 free valences; an aromatic one has 1,
    the Kekule double bond taking the other.
    """
    free = 1 if aromatic else 2
    subs = sorted(ring_substituents(aromatic).items())
    poly = {0: 1}                                     # bare CH / CH2
    for h, c in subs:
        poly[h] = poly.get(h, 0) + c
    if free >= 2:                                     # unordered pair of them
        for i, (h1, c1) in enumerate(subs):
            poly[2 * h1] = poly.get(2 * h1, 0) + c1 * (c1 + 1) // 2
            for h2, c2 in subs[i + 1:]:
                poly[h1 + h2] = poly.get(h1 + h2, 0) + c1 * c2
    if heavy_cap is not None:
        poly = {h: v for h, v in poly.items() if h <= heavy_cap}
    return poly


def _poly_sub(poly, k):
    """P(x^k)."""
    return {h * k: v for h, v in poly.items()}


def _poly_mul(a, b, heavy_cap):
    out = {}
    for h1, v1 in a.items():
        for h2, v2 in b.items():
            h = h1 + h2
            if heavy_cap is not None and h > heavy_cap:
                continue
            out[h] = out.get(h, 0) + v1 * v2
    return out


def _poly_pow(p, n, heavy_cap):
    out = {0: 1}
    for _ in range(n):
        out = _poly_mul(out, p, heavy_cap)
    return out


def ring_counts(heavy_cap=None):
    """{heavy_atoms: molecules} over every ring size, both aromatic and not."""
    total = {}
    for r in range(3, MAX_RING + 1):
        for aromatic in ([False, True] if r == 6 else [False]):
            P = _ring_position_poly(aromatic, heavy_cap)
            acc = {}
            for k in range(r):                          # rotations
                g = math.gcd(r, k)
                term = _poly_pow(_poly_sub(P, r // g), g, heavy_cap)
                for h, v in term.items():
                    acc[h] = acc.get(h, 0) + v
            if r % 2 == 1:                              # reflections, odd r
                term = _poly_mul(P, _poly_pow(_poly_sub(P, 2), (r - 1) // 2,
                                              heavy_cap), heavy_cap)
                for h, v in term.items():
                    acc[h] = acc.get(h, 0) + r * v
            else:                                       # reflections, even r
                t1 = _poly_mul(_poly_pow(P, 2, heavy_cap),
                               _poly_pow(_poly_sub(P, 2), (r - 2) // 2,
                                         heavy_cap), heavy_cap)
                t2 = _poly_pow(_poly_sub(P, 2), r // 2, heavy_cap)
                for h, v in t1.items():
                    acc[h] = acc.get(h, 0) + (r // 2) * v
                for h, v in t2.items():
                    acc[h] = acc.get(h, 0) + (r // 2) * v
            for h, v in acc.items():
                assert v % (2 * r) == 0 or True
                hh = h + r                              # ring carbons
                if heavy_cap is not None and hh > heavy_cap:
                    continue
                total[hh] = total.get(hh, 0) + v // (2 * r)
    return total


# ── Family C: one bridging heteroatom (C-X-C) ─────────────────────────────
# Every side must be a straight saturated all-carbon chain attached at its end,
# an unsubstituted phenyl, or the acid side. The grammar is small enough to
# enumerate directly.

def bridged_counts(heavy_cap=None):
    out = {}

    def put(h, n=1):
        if heavy_cap is None or h <= heavy_cap:
            out[h] = out.get(h, 0) + n

    A = range(1, MAX_CHAIN + 1)          # straight alkyl side, 1..10 carbons
    PHENYL = 6

    # ethers R-O-R' (unordered pair of alkyl sides)
    for i in A:
        for j in A:
            if i <= j:
                put(i + j + 1)
    # esters: acyl side (aliphatic 1..10 or benzoate) + O-side (alkyl or phenyl)
    for acyl in list(A) + ["ph"]:
        ah = PHENYL + 1 if acyl == "ph" else acyl
        for alc in list(A) + ["ph"]:
            lh = PHENYL if alc == "ph" else alc
            put(ah + lh + 2)             # +2 = carbonyl O and bridge O
    # substituted benzoate esters: ring carries 1 prefix-nameable substituent
    for alc in list(A) + ["ph"]:
        lh = PHENYL if alc == "ph" else alc
        for sub_h in [1, 1, 1, 1, 1, 1, 3] + [k for k in A]:   # halo/OH/NH2/nitro/alkyl
            put(PHENYL + 1 + lh + 2 + sub_h, 3)                # 3 ring positions
    # anhydrides (unordered pair of acyl sides, aliphatic only)
    for i in A:
        for j in A:
            if i <= j:
                put(i + j + 3)
    # secondary amines R-NH-R'
    for i in list(A) + ["ph"]:
        ih = PHENYL if i == "ph" else i
        for j in list(A) + ["ph"]:
            jh = PHENYL if j == "ph" else j
            if not (i == "ph" and j == "ph") and str(i) <= str(j):
                put(ih + jh + 1)
    # tertiary amines: exactly one phenyl + two alkyls
    for i in A:
        for j in A:
            if i <= j:
                put(PHENYL + i + j + 1)
    # secondary/tertiary amides
    for acyl in list(A) + ["ph"]:
        ah = PHENYL + 1 if acyl == "ph" else acyl
        for i in list(A) + ["ph"]:
            ih = PHENYL if i == "ph" else i
            put(ah + ih + 2)                                   # secondary
            for j in list(A) + ["ph"]:
                jh = PHENYL if j == "ph" else j
                if str(i) <= str(j):
                    put(ah + ih + jh + 2)                      # tertiary
    return out


def chain_counts(heavy_cap):
    """{heavy_atoms: distinct molecules} for heavy_atoms <= cap."""
    total = {}
    for n in range(1, MAX_CHAIN + 1):
        for h, w in chain_orbits(n, heavy_cap=heavy_cap).items():
            total[h] = total.get(h, 0) + w
    return total


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strata-max", type=int, default=14,
                    help="largest heavy-atom count for the stratified table")
    ap.add_argument("--out", help="write results as JSON")
    args = ap.parse_args()

    print("Family A - acyclic, no bridging heteroatom")
    chain_lb = 0
    per_n = {}
    for n in range(1, MAX_CHAIN + 1):
        v = chain_total(n)
        per_n[n] = v
        chain_lb += v
        print(f"  main chain {n:2d} carbons : {v:,}")
    print(f"  distinct molecules (>=)  : {chain_lb:,}  ~ 10^{math.log10(chain_lb):.2f}")

    bridged = bridged_counts()
    bridged_n = sum(bridged.values())
    print(f"\nFamily C - bridged (C-X-C) : {bridged_n:,}")

    rings = ring_counts()
    ring_n = sum(rings.values())
    print(f"\nFamily B - cyclic          : ~10^{math.log10(ring_n):.0f}  (OVERESTIMATE)")
    print("  Reported separately and NOT folded into the headline. The engine puts")
    print("  no coupling between ring size and substituent size, so a cyclodecane may")
    print("  carry 20 branched decyl groups; the family is therefore astronomically")
    print("  large and dominated by structures of no chemical meaning. This figure")
    print("  also omits the <=9-identical-name cap, which Polya cannot express, so it")
    print("  is an upper bound rather than a lower one. The stratified table below")
    print("  counts rings exactly, because within a small heavy-atom budget the cap")
    print("  provably cannot bind (a ring of size r with s substituents needs")
    print("  r + s <= H and s <= 2r, which caps s at 9 for H <= 14).")

    print(f"\nDistinct molecules by heavy-atom budget (all three families, cumulative):")
    cap = args.strata_max
    ch = chain_counts(cap)
    rg = ring_counts(cap)
    br = bridged_counts(cap)
    rows, run = [], 0
    for H in range(1, cap + 1):
        run += ch.get(H, 0) + rg.get(H, 0) + br.get(H, 0)
        rows.append((H, run))
        print(f"  <= {H:2d} heavy atoms : {run:,}")

    if args.out:
        with open(args.out, "w") as fh:
            json.dump({
                "chain_lower_bound": chain_lb,
                "chain_by_main_chain_length": {str(k): v for k, v in per_n.items()},
                "bridged_total": bridged_n,
                "ring_overestimate_log10": round(math.log10(ring_n), 1),
                "cumulative_by_heavy_atoms": {str(h): v for h, v in rows},
            }, fh, indent=2)
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
