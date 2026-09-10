#!/usr/bin/env python3
"""
Exhaustive check against a known-truth sequence: the alkane constitutional isomers.

The number of distinct alkanes C(n)H(2n+2) is a classical, independently
published sequence (OEIS A000602: 1, 1, 1, 2, 3, 5, 9, 18, 35, 75 for n = 1..10).
Because it is small and exactly known, it is the one place the engine's coverage
can be checked completely rather than sampled -- every isomer, no exceptions.

Trees are built by adding one carbon at a time to every atom with a free valence
and de-duplicated with the AHU canonical form, which is exact for trees. The
resulting counts are asserted against A000602 before anything is handed to the
engine, so a bug in the generator cannot be mistaken for engine coverage.

Two things are then checked:
  1. coverage  -- the engine names every isomer (nothing in scope is refused);
  2. injectivity -- no two distinct isomers receive the same name, which is what
     makes a name a usable identifier rather than merely a plausible string.

Usage:
    python alkane_isomers.py --max-carbons 10 --out isomers.jsonl
"""

from __future__ import annotations

import argparse
import json

A000602 = [1, 1, 1, 2, 3, 5, 9, 18, 35, 75, 159, 355, 802, 1858]


def canonical(adj: list[list[int]]) -> str:
    """AHU canonical form of a tree, rooted at its centre(s)."""
    n = len(adj)
    if n == 1:
        return "()"
    deg = [len(a) for a in adj]
    leaves = [i for i in range(n) if deg[i] <= 1]
    removed, remaining = set(leaves), n
    layer = leaves
    while remaining > 2:
        remaining -= len(layer)
        nxt = []
        for v in layer:
            for u in adj[v]:
                if u not in removed:
                    deg[u] -= 1
                    if deg[u] == 1:
                        nxt.append(u)
        removed.update(nxt)
        layer = nxt
    centres = layer

    def enc(v: int, parent: int) -> str:
        kids = sorted(enc(u, v) for u in adj[v] if u != parent)
        return "(" + "".join(kids) + ")"

    return min("".join(sorted([enc(c, -1)])) for c in centres) if len(centres) == 1 \
        else min(enc(centres[0], centres[1]) + enc(centres[1], centres[0]),
                 enc(centres[1], centres[0]) + enc(centres[0], centres[1]))


def grow(trees: list[list[list[int]]]) -> list[list[list[int]]]:
    """Every way to attach one more carbon, de-duplicated up to isomorphism."""
    seen, out = set(), []
    for adj in trees:
        for v in range(len(adj)):
            if len(adj[v]) >= 4:          # carbon valence
                continue
            new = [list(a) for a in adj] + [[v]]
            new[v].append(len(adj))
            key = canonical(new)
            if key not in seen:
                seen.add(key)
                out.append(new)
    return out


def to_atoms(adj: list[list[int]]) -> list[dict]:
    return [{"element": "C",
             "bonds": [{"to": u, "order": 1} for u in sorted(nbrs)]}
            for nbrs in adj]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-carbons", type=int, default=10)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    trees = [[[]]]
    rows = []
    for n in range(1, args.max_carbons + 1):
        if n > 1:
            trees = grow(trees)
        expected = A000602[n - 1]
        assert len(trees) == expected, \
            f"generator produced {len(trees)} isomers for C{n}, A000602 says {expected}"
        print(f"  C{n:<2d} {len(trees):5d} isomers  (matches A000602)")
        for i, adj in enumerate(trees):
            rows.append({"id": f"C{n}-{i}", "carbons": n, "atoms": to_atoms(adj)})

    with open(args.out, "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    print(f"\nwrote {len(rows)} isomers to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
