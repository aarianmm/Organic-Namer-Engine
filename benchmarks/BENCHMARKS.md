# Organic-Namer-Engine — benchmarks

Measurements of the IUPAC naming engine in `src/OrganicNamer.Core`: how large a
molecule space its grammar covers, how often the names it produces are correct,
and how fast it produces them.

Everything here is reproducible from this directory — see [`README.md`](README.md)
for the commands. Raw outputs are in `results/`.

**Environment.** macOS 15.7.4, Apple Silicon (Arm64), 8 cores, .NET 10.0.9,
Python 3.13 + RDKit 2026.03.6, OPSIN 2.9.0, Java 23. Engine at commit `62553da`
plus the carbon-subsuming prefix fix and the epoxide/oxa-bridge ring-routing
fix, both described in §3.4.

---

## 1. Summary

| Measure | Result |
|---|---|
| Structural round-trip accuracy, generated corpus | **100.00%** (49,988 / 49,988) |
| Structural round-trip accuracy, real PubChem compounds | **100.00%** (332 / 332) |
| Exhaustive alkane coverage, C₁–C₁₀ | **99.33%** (149 / 150 isomers) |
| Name collisions across 635 distinct isomers | **0** |
| Name agreement with PubChem — strict / convention-normalised | **34.94% / 94.58%** |
| Throughput, single-threaded | **12,375 molecules/s** (p50 24 µs, p99 929 µs) |
| Distinct molecules nameable, ≤ 11 heavy atoms | **≥ 86,023,812** |
| Distinct molecules nameable, acyclic family | **≥ 1.12 × 10²⁵** |

Every name the engine emits — across both the 49,988 generated molecules it
accepts and the 332 real PubChem compounds it accepts — denotes exactly the
molecule it was given, verified by an independent parser, not by string
comparison (§3.1). There are no known false accepts on either corpus (§3.4).

---

## 2. Size of the nameable molecule space

`tools/combinatorics.py` counts the space analytically — no molecules are
generated and the engine is never run. The rules were read out of `IUPAC.cs`,
`ElementGraph.cs` and `SpecificationData.cs`, then each one was confirmed against
the live engine (§2.3).

Molecules are partitioned by the engine's own three-way routing: **chain**
(acyclic, no bridging heteroatom), **cyclic** (one carbon ring), and **bridged**
(one chain-breaking `C–X–C` heteroatom). The chain family is walked with a
transfer matrix along the main chain; chain reversal is quotiented out with
Burnside's lemma, so a molecule and its mirror numbering count once. The
engine's "at most 9 identical substituent names" cap (`numericalPrefixes` holds
`di`…`nona`) is a whole-molecule constraint and is applied by inclusion–exclusion
over the substituent buckets that can physically reach ten occurrences.

### 2.1 By heavy-atom budget

The useful figure. Cumulative count of distinct constitutional isomers the
grammar admits, across all three families:

| Heavy atoms ≤ | Distinct molecules |
|---:|---:|
| 6 | 13,475 |
| 8 | 439,088 |
| 10 | 14,747,508 |
| **11** | **86,023,812** |
| 12 | 504,174,744 |
| 14 | 17,391,816,737 |

For scale, GDB-11 — the standard enumeration of small organic molecules up to 11
heavy atoms of C, N, O and F — contains about 26.4 million structures. This
engine's grammar is broader at that size mainly because it also admits S, Cl, Br
and I, and applies no chemical-stability filter.

### 2.2 Whole-grammar totals

| Family | Distinct molecules |
|---|---:|
| Chain (acyclic) | ≥ 1.12 × 10²⁵ |
| Bridged (`C–X–C`) | 1,759 |
| Cyclic | ~10⁸¹ *(overestimate — see below)* |

The cyclic figure is reported but **should not be quoted**. The engine imposes
no coupling between ring size and substituent size, so a cyclodecane may legally
carry twenty branched decyl substituents — verified: `decylcyclohexane` and
`1,2,3,4,5,6-hexahexyl-1,2,3,4,5,6-hexapentylcyclohexane` both name successfully.
The family is therefore astronomically large and dominated by structures of no
chemical meaning, and the Pólya enumeration used for it cannot express the ≤9 cap,
making the number an upper rather than a lower bound. The heavy-atom-stratified
table above counts rings exactly, because within a small budget the cap provably
cannot bind (a ring of size *r* with *s* substituents needs *r + s ≤ H* and
*s ≤ 2r*, capping *s* at 9 for *H ≤ 14*).

### 2.3 Why the counts are trustworthy

Every figure is a **lower bound**, and the size of the gap is measured rather
than assumed.

The model requires the main chain to be a longest path in the molecule. The
engine is slightly more permissive: it discards branch-invalid candidate chains
*first*, then takes the longest survivor — so a molecule whose longest path is
unnameable but which has a shorter nameable chain is counted here as zero and
named by the engine anyway. The bound therefore always errs low.

`tools/alkane_isomers.py` measures that gap exactly. It builds every alkane
constitutional isomer up to C₁₂ by canonical-form de-duplication, asserts the
counts against the independently published sequence OEIS A000602
(1, 1, 1, 2, 3, 5, 9, 18, 35, 75, 159, 355), and runs all 664 through the engine:

| | C₉ | C₁₀ | C₁₁ | C₁₂ |
|---|---:|---:|---:|---:|
| Analytic model | 35 | 72 | 146 | 307 |
| Engine actually names | 35 | 74 | 154 | 332 |
| Bound tightness | 100% | 97.3% | 94.8% | 92.5% |

Two further checks:

- **The scope rules are the engine's own.** A separate structural criterion —
  *"the longest branch-valid leaf-to-leaf path is ≤ 10 carbons"* — reproduces the
  engine's accept/reject decision on **all 664 alkane isomers exactly**, C₁ through
  C₁₂, including the single decane and the eight undecanes it refuses.
- **The grammar model round-trips.** `tools/generate_corpus.py` builds molecules
  from the same rules `combinatorics.py` counts. The engine names **100.00%** of
  50,000 sampled molecules (49,995 named; the 5 rejections are the ≤9 cap firing
  on legitimately over-substituted samples). An earlier draft scored 99.35%, and
  the failures pinpointed two real scope errors in the model, both since fixed:
  a branch tip may not carry a nitro group (its oxygens are not adjacent to the
  spine tip), and an aromatic ring's branch tips are restricted to halogens.

---

## 3. Accuracy

### 3.1 Structural round-trip (primary measure)

String-comparing against a reference name conflates two very different things:
a name that is *wrong*, and a name that is *right but written under a different
convention*. Composing the engine with a parser separates them:

```
graph --engine--> name --OPSIN--> graph'
```

OPSIN (Open Parser for Systematic IUPAC Nomenclature, Lowe et al., Cambridge)
converts a name back to a structure. If `canonical(graph) == canonical(graph')`
the name unambiguously denotes the molecule it was generated from, whatever
convention it used. Comparison is constitutional only — stereochemistry is
stripped and nitro groups are normalised to one tautomer, since OPSIN emits the
charge-separated form and the engine's input model uses the neutral hypervalent
one. Where the engine emits several equally-valid names, the round-trip passes if
any one of them resolves correctly.

| Corpus | Molecules | Round-trip match | OPSIN parse failures |
|---|---:|---:|---:|
| Generated from the grammar | 49,988 | **100.00%** (49,988) | 0 |
| Real compounds (PubChem CIDs 1–6000) | 332 | **100.00%** (332) | 0 |

Zero parse failures on 50,382 generated names is itself a result: every name the
engine produced was well-formed enough for an independent parser to read.

The generated corpus is exhaustive over the engine's own grammar, so 100% here
means something specific: for every molecule the engine accepts, the name it
emits resolves back to that exact molecule. The PubChem corpus now round-trips
at 100% too — the epoxide false accept described in §3.4 was a ring-routing
defect, since fixed, and there are no known false accepts left on either
corpus.

### 3.2 Exhaustive alkane coverage

The one place coverage can be checked *completely* rather than sampled, since the
ground-truth count is known exactly:

| | C₁–C₉ | C₁₀ | C₁–C₁₀ | C₁₁ | C₁₂ |
|---|---:|---:|---:|---:|---:|
| Isomers | 75 | 75 | 150 | 159 | 355 |
| Named | 75 | 74 | **149 (99.33%)** | 154 | 332 |

The single decane the engine cannot name is the one isomer that is not a
*caterpillar* tree — its branch points do not lie on a common path, so every
candidate main chain leaves a branch hanging off a branch. That is the documented
"branches on branches: not modelled" limitation, and C₁₀ is exactly where it
first bites. C₁₁ and C₁₂ losses are the same limitation plus the 10-carbon root
name table.

**No two distinct isomers received the same name** across all 635 named
structures — the engine's names are usable as identifiers, not merely as
plausible strings.

### 3.3 Agreement with PubChem's reference names

| Comparison | Matches |
|---|---:|
| Strict exact string | 116 / 332 = **34.94%** |
| After convention normalisation | 314 / 332 = **94.58%** |

The gap is almost entirely one habit: the engine never elides a locant that IUPAC
treats as redundant, writing `butan-1-al` for `butanal` and `hexane-1,6-dioic
acid` for `hexanedioic acid`. That single rule accounts for 174 of the 198
normalised matches; a retained-name alias table (`ethanoic acid` ↔ `acetic acid`,
`phenylamine` ↔ `aniline`, …) accounts for 29 more.

Normalisation is purely lexical and table-driven — never structural, never
fuzzy-matched — and is verified not to mask any genuine error that §3.1
identified independently. `score_names.py` asserts this three ways and fails if
any of them breaks: the set of known genuine mismatches must stay empty, the
twelve names fixed by the `oxo` rename must stay correct, and the six molecules
the engine now refuses must stay refused. The residual 5.71% is dominated by
retained trivial names the alias table does not carry — `oxaldehydic acid` for
the engine's (structurally correct) `2-oxoethan-1-oic acid`, for instance. Those
are left unaliased on purpose so this figure stays a measurement rather than
something tuned.

### 3.4 What the failures were, and what was fixed

This section analyses the 97 round-trip failures measured at commit `62553da`,
before the fix. **96 of the 97 shared a single root cause**: carbon-subsuming
prefixes placed on the group's own chain carbon.

`-CHO` (`formyl`), `-C≡N` (`cyano`) and `-CONH₂` (`carbamoyl`) each include their
own carbon in the prefix name. When such a group sits on the main chain and is
not the principal suffix, `FindPrefixCarbons` cites it at that chain carbon's
locant, so the name denotes a molecule with one carbon too many:

| | |
|---|---|
| Input | `O=CCC(=O)O` — 3 carbons (3-oxopropanoic acid) |
| Engine | `3-formylpropan-1-oic acid` |
| That name means | `O=CCCC(=O)O` — 4 carbons |

The engine already guards against exactly this for *branch* substituents — the
`Q3` comment in `NameCarbonSubstituent` explains that carbon-subsuming prefixes
are structurally wrong there and rejects them — but the main-chain path has no
equivalent check. It should emit `oxo` for a chain-terminal `C=O`. Breakdown of
the 79 generated-corpus failures: `formyl` 72, `cyano` 5, `chlorocarbonyl` 2.

Across both corpora, **every** name the engine ever emitted carrying one of these
prefixes failed the round-trip — there was no partial-correctness case to
preserve:

| Prefix | Names emitted | Round-trip matches |
|---|---:|---:|
| `formyl` | 84 | **0** |
| `cyano` | 5 | **0** |
| `carbamoyl` | 5 | **0** |
| `chlorocarbonyl` | 2 | **0** |
| `carboxy` | 0 | — |
| **Total** | **96** | **0** |

The remaining failure is a **false accept**, which matters more than its count
suggests because it violates the engine's stated "name it exactly or reject it"
contract. PubChem CID 2859 (`7-oxabicyclo[4.1.0]heptane-2,3,4,5-tetrol`, an
epoxide-bridged bicyclic) is named `1,2,3,4,5,6-hexahydroxycyclohexane` — the
epoxide bridge is silently dropped and two hydroxyls invented. Fused/bridged ring
systems are supposed to be refused.

The carbon-subsuming prefix defect has since been fixed: an on-chain carbonyl is
now always named `oxo` rather than `formyl`, and a `COOH`/`COCl`/`CONH₂`/`C≡N`
group on the main chain that is not the principal characteristic group is now
refused rather than misnamed. All figures elsewhere in this document are
post-fix measurements; the 97 failures analysed above no longer occur, and the
five molecules that can no longer be named are reflected in the reduced corpus
sizes (49,988 generated, 333 PubChem).

**The epoxide false accept has since been fixed.** Its root cause was a routing
test, `etherOnRing` in `IUPAC.cs`, added so molecules like `methoxycyclohexane`
reach the cyclic naming path: it checked that a carbon ring existed and that
the bridge was an oxygen, but never that the bridge's own two attachment
carbons sat on that ring. For a genuine ether they don't; for an epoxide, or
any oxa-bridge across a ring, they do — and `etherOnRing` routed those around a
guard that already existed for exactly this shape. That guard lives in
`SplitAtBridgingAtom` (`ElementGraph.cs`), whose comment calls itself "the ONLY
place they are caught," because `IsCyclic()` only walks the carbon subgraph and
cannot see a ring closed through a heteroatom. The fix adds one conjunct,
`!atoms.BridgeClosesARing(...)`, to `etherOnRing`, so these molecules fall
through to the bridged path and hit that guard instead. No new error message;
the diff is 15 lines added, 1 changed.

Three shapes were affected, with two visibly different corruptions depending on
input atom order:

| SMILES | Was named |
|---|---|
| `C1CCC2OC2C1` (cyclohexene oxide) | `cycloheptane` — the bridge oxygen absorbed into the ring and counted as a carbon |
| `C1CC2CCC1O2` (1,4-epoxycyclohexane) | `1,4-dihydroxycyclohexane` — the single bridging oxygen read as two independent hydroxyls |
| `C1(C(C(C2C(C1O)O2)O)O)O` (PubChem CID 2859) | `1,2,3,4,5,6-hexahydroxycyclohexane` |

All three now hit the same clean rejection as any other ring closed through a
heteroatom: `ChemistryError: Rings containing a heteroatom are not supported`.
CID 2859 is the molecule that drops the PubChem-named count elsewhere in this
document from the 333 quoted earlier in this section to 332. **There are no
known false accepts on either corpus.**

### 3.5 Scope discrimination

Of 4,747 PubChem compounds converted successfully, the engine named 332 (7.0%)
and refused 4,415. A low acceptance rate is the intended behaviour, not a
weakness — the supported vocabulary is deliberately A-level-sized, and the design
rule is to refuse anything outside it rather than guess. The count that matters
is that refusals are *clean*: every one of the 4,415 came back as a caught
`ChemistryError`, with no crashes, hangs, or silently wrong names.

Five of those refusals are new, and are a deliberate trade: a non-principal
`COOH`/`COCl`/`CONH₂`/`C≡N` on the main chain has no correct substituent name at
this level of the grammar, so the engine refuses instead of guessing. Recovering
them requires excluding those carbons from chain selection, which is a separate
piece of work.

A caveat worth stating: this benchmark can measure false accepts — §3.4 found
one, and it is now fixed, so the count stands at zero — but cannot cheaply
measure false *rejects*, since that would need a ground-truth label for whether
each of the 4,415 is genuinely out of scope. That fix is also why the refused
count here is one higher than in earlier measurements: CID 2859 (§3.4) now
refuses cleanly instead of round-tripping to a wrong name.

---

## 4. Performance

Measured in-process against `OrganicNamer.Core` — no HTTP, no serialisation in
the timed region. All 50,000 corpus molecules are parsed into memory first; the
timed region covers only `FillImplicitHydrogens` → `FromJsonAtoms` →
`new IUPAC(...)`. Three warm-up passes are discarded, then ten timed passes.
Single-threaded by design, so the latency figures are honest per-molecule costs.

| | |
|---|---|
| Throughput | **12,375 molecules/s** (σ = 89 across 10 reps) |
| Latency mean | 80.8 µs |
| p50 | 23.5 µs |
| p99 | 929 µs |
| Max | 13.9 ms |

The heavy tail is chain selection: `FindEveryLongestPath` enumerates every
leaf-to-leaf path and each is filtered for branch validity, so cost grows sharply
with the number of chain ends. The p50-to-max spread of ~590× is that effect.

---

## 5. Reproducing

See [`README.md`](README.md). In short:

```bash
# scope
python tools/combinatorics.py --strata-max 14 --out results/combinatorics.json
python tools/alkane_isomers.py --max-carbons 12 --out data/alkane_isomers.jsonl

# accuracy
python tools/generate_corpus.py --count 50000 --seed 7 --out data/generated_corpus.jsonl
dotnet run -c Release --project runner/OrganicNamer.Benchmarks -- \
    name --in data/generated_corpus.jsonl --out results/generated_named.jsonl
python tools/opsin_roundtrip.py --in results/generated_named.jsonl \
    --out results/generated_roundtrip.jsonl
python tools/score_names.py --in results/pubchem_named.jsonl \
    --out results/name_scores.jsonl

# performance
dotnet run -c Release --project runner/OrganicNamer.Benchmarks -- \
    bench --in data/generated_corpus.jsonl --out results/throughput.json
```

The PubChem fetch (`tools/fetch_pubchem.py`) caches to `data/pubchem_raw.jsonl`
and supports `--resume`, so the network step runs once.
