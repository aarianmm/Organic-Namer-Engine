#!/usr/bin/env python3
"""
Sample molecules from the engine's supported grammar.

Two jobs at once:

1. It is the large-N corpus for the OPSIN round-trip accuracy measurement. The
   PubChem set only contains a few hundred molecules inside the engine's scope,
   which is too few to characterise the tail; sampling the grammar directly
   gives tens of thousands, spread evenly over every path the engine supports
   rather than clustered on whatever PubChem happens to hold at low CIDs.

2. It empirically checks the scope model used by combinatorics.py. That script
   counts the molecule space from rules read out of the source; this one builds
   molecules from the *same* rules and feeds them to the real engine. If the
   model is faithful, essentially everything generated should be named. A
   generated molecule the engine rejects means the analytic model is too
   generous somewhere, and the count would be an overestimate.

Molecules are emitted in the engine's heavy-atom JSON format (hydrogens are
filled in by the engine from valency).

Usage:
    python generate_corpus.py --count 50000 --seed 1 --out corpus.jsonl
"""

from __future__ import annotations

import argparse
import json
import random

MAX_CHAIN = 10
MAX_RING = 10
MAX_MULT = 9

MONO = ["F", "Cl", "Br", "I", "OH", "SH", "NH2", "NO2"]
DI = ["O", "NH"]          # =O, =NH

# A branch tip may not carry nitro: NameCarbonSubstituent demands every heteroatom
# in a substituent be adjacent to the spine tip, and nitro's oxygens sit one bond
# further out. Aromatic ring branches are stricter still -- TipGroupPrefix refuses
# suffix-capable groups there (allowSuffixCapableTip: false), leaving halogens only.
LEG_TIPS = ["F", "Cl", "Br", "I", "OH", "SH", "NH2"]
AROMATIC_LEG_TIPS = ["F", "Cl", "Br", "I"]


class Mol:
    def __init__(self):
        self.atoms = []

    def add(self, el):
        self.atoms.append({"element": el, "bonds": []})
        return len(self.atoms) - 1

    def bond(self, i, j, order=1):
        self.atoms[i]["bonds"].append({"to": j, "order": order})
        self.atoms[j]["bonds"].append({"to": i, "order": order})

    def used(self, i):
        return sum(b["order"] for b in self.atoms[i]["bonds"])

    def free(self, i):
        return 4 - self.used(i)


def put_mono(m, c, kind):
    """Attach a monovalent group, consuming exactly one valence on carbon c."""
    if kind in ("F", "Cl", "Br", "I"):
        m.bond(c, m.add(kind))
    elif kind == "OH":
        m.bond(c, m.add("O"))
    elif kind == "SH":
        m.bond(c, m.add("S"))
    elif kind == "NH2":
        m.bond(c, m.add("N"))
    elif kind == "NO2":
        n = m.add("N")
        m.bond(c, n)
        m.bond(n, m.add("O"), 2)
        m.bond(n, m.add("O"), 2)


def put_di(m, c, kind):
    m.bond(c, m.add("O" if kind == "O" else "N"), 2)


def decorate(m, c, rng, budget, counts, allow_multi=True):
    """Fill up to `budget` free valences on carbon c with functional groups."""
    while budget > 0:
        r = rng.random()
        if r < 0.55:
            return
        if allow_multi and budget >= 3 and rng.random() < 0.06:
            m.bond(c, m.add("N"), 3)          # nitrile
            return
        if allow_multi and budget >= 2 and rng.random() < 0.22:
            k = rng.choice(DI)
            if counts[k] >= MAX_MULT:
                return
            counts[k] += 1
            put_di(m, c, k)
            budget -= 2
            continue
        k = rng.choice(MONO)
        if counts[k] >= MAX_MULT:
            return
        counts[k] += 1
        put_mono(m, c, k)
        budget -= 1


def gen_chain(rng):
    m = Mol()
    n = rng.randint(1, MAX_CHAIN)
    carbons = [m.add("C") for _ in range(n)]
    # backbone bond orders
    for i in range(n - 1):
        r = rng.random()
        order = 1 if r < 0.80 else (2 if r < 0.94 else 3)
        if m.free(carbons[i]) < order or m.free(carbons[i + 1]) < order:
            order = 1
        m.bond(carbons[i], carbons[i + 1], order)

    counts: dict[str, int] = {k: 0 for k in MONO + DI}
    leg_counts: dict[int, int] = {}
    reach = 1
    for i in range(1, n + 1):
        c = carbons[i - 1]
        base = 1 if i == 1 else reach + 1
        longest_leg = 0
        # legs: a leg of L carbons must not create a path longer than the chain
        cap = max(0, min(i - 1, n - base))
        while cap > 0 and m.free(c) > 0 and rng.random() < 0.22:
            L = rng.randint(1, cap)
            if leg_counts.get(L, 0) >= MAX_MULT:
                break
            leg_counts[L] = leg_counts.get(L, 0) + 1
            prev = c
            for _ in range(L):
                nc = m.add("C")
                m.bond(prev, nc)
                prev = nc
            # a leg may carry one tip group, and only at its tip
            if rng.random() < 0.18:
                put_mono(m, prev, rng.choice(LEG_TIPS))
            longest_leg = max(longest_leg, L)
            cap = min(cap, n - 1 - L)      # two legs here must also fit
        reach = max(base, longest_leg + 1)
        decorate(m, c, rng, m.free(c), counts)
    return m


def gen_ring(rng):
    m = Mol()
    size = rng.randint(3, MAX_RING)
    aromatic = size == 6 and rng.random() < 0.45
    ring = [m.add("C") for _ in range(size)]
    for k in range(size):
        order = 1
        if aromatic:
            order = 2 if k % 2 == 0 else 1
        m.bond(ring[k], ring[(k + 1) % size], order)

    counts: dict[str, int] = {k: 0 for k in MONO + DI}
    n_sub = rng.randint(0, min(size, 4))
    used_alkoxy = False
    for c in rng.sample(ring, n_sub):
        if m.free(c) <= 0:
            continue
        r = rng.random()
        if r < 0.40:
            k = rng.choice(MONO)
            if counts[k] < MAX_MULT:
                counts[k] += 1
                put_mono(m, c, k)
        elif r < 0.55 and not used_alkoxy:
            used_alkoxy = True                     # only ONE bridge allowed
            o = m.add("O")
            m.bond(c, o)
            prev = o
            for _ in range(rng.randint(1, 4)):
                nc = m.add("C")
                m.bond(prev, nc)
                prev = nc
        elif r < 0.85:
            prev = c                               # straight alkyl
            for _ in range(rng.randint(1, 4)):
                nc = m.add("C")
                m.bond(prev, nc)
                prev = nc
            if rng.random() < 0.15:
                put_mono(m, prev, rng.choice(
                    AROMATIC_LEG_TIPS if aromatic else LEG_TIPS))
        else:
            a = m.add("C")                         # branched alkyl, methyls only
            m.bond(c, a)
            for _ in range(rng.randint(2, 3)):
                m.bond(a, m.add("C"))
    return m


def straight_alkyl(m, start, length):
    prev = start
    for _ in range(length):
        nc = m.add("C")
        m.bond(prev, nc)
        prev = nc
    return prev


def phenyl(m, attach_to):
    ring = [m.add("C") for _ in range(6)]
    for k in range(6):
        m.bond(ring[k], ring[(k + 1) % 6], 2 if k % 2 == 0 else 1)
    m.bond(attach_to, ring[0])
    return ring


def gen_bridged(rng):
    """One chain-breaking heteroatom: ether, ester, amide, amine, anhydride."""
    m = Mol()
    kind = rng.choice(["ether", "ester", "amide2", "amine2", "amine3", "anhydride"])
    if kind == "ether":
        o = m.add("O")
        straight_alkyl(m, o, rng.randint(1, 5))
        straight_alkyl(m, o, rng.randint(1, 5))
    elif kind == "ester":
        acid = m.add("C")
        m.bond(acid, m.add("O"), 2)
        if rng.random() < 0.30:
            phenyl(m, acid)                        # benzoate
        else:
            straight_alkyl(m, acid, rng.randint(0, 4))
        o = m.add("O")
        m.bond(acid, o)
        if rng.random() < 0.20:
            phenyl(m, o)
        else:
            straight_alkyl(m, o, rng.randint(1, 5))
    elif kind == "amide2":
        acid = m.add("C")
        m.bond(acid, m.add("O"), 2)
        straight_alkyl(m, acid, rng.randint(0, 4))
        n = m.add("N")
        m.bond(acid, n)
        if rng.random() < 0.25:
            phenyl(m, n)
        else:
            straight_alkyl(m, n, rng.randint(1, 4))
    elif kind == "amine2":
        n = m.add("N")
        straight_alkyl(m, n, rng.randint(1, 5))
        if rng.random() < 0.30:
            phenyl(m, n)
        else:
            straight_alkyl(m, n, rng.randint(1, 5))
    elif kind == "amine3":
        n = m.add("N")
        phenyl(m, n)                               # exactly one phenyl required
        straight_alkyl(m, n, rng.randint(1, 3))
        straight_alkyl(m, n, rng.randint(1, 3))
    else:                                          # anhydride
        o = m.add("O")
        for _ in range(2):
            c = m.add("C")
            m.bond(o, c)
            m.bond(c, m.add("O"), 2)
            straight_alkyl(m, c, rng.randint(0, 3))
    return m


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--count", type=int, default=50000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    gens = [("chain", gen_chain, 0.60), ("ring", gen_ring, 0.25),
            ("bridged", gen_bridged, 0.15)]
    with open(args.out, "w") as fh:
        for i in range(args.count):
            r = rng.random()
            acc = 0.0
            for family, fn, w in gens:
                acc += w
                if r <= acc:
                    break
            m = fn(rng)
            fh.write(json.dumps({"id": f"gen-{i}", "family": family,
                                 "atoms": m.atoms}) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
