# Supported Molecules

Scope reference for the naming engine (`src/OrganicNamer.Core`), derived from the code, not
from the plan documents. Rule: **name it exactly or reject it** — an unsupported shape
throws a `ChemistryError` (HTTP 422), it is never guessed at.

## 1. Input model

- Molecules arrive as a JSON atom list: `{ "element": "C", "bonds": [{ "to": i, "order": n }] }`.
- Elements known to the periodic table: **C, H, O, N, F, Cl, Br, I, S**. Anything else → reject.
- Hydrogens may be omitted; they are filled from valency deficit. Explicit H is also fine.
- Bond orders 1/2/3. Atoms exceeding their valency are rejected (only the nitro N is exempt).
- **Aromaticity must be drawn Kekulé-style**: alternating 1,2,1,2,1,2 bond orders around a
  6-carbon ring. A 6-ring drawn with all single bonds is cyclohexane, not benzene.
- **Nitro** is accepted *only* as the hypervalent neutral form `C–N(=O)(=O)`. The
  charge-separated form N⁺(–O⁻)(=O) is unrepresentable (no formal-charge field) → reject.

### Specification sets

| Set | Vocabulary |
|---|---|
| `AllGroups` (default) | everything below |
| `Hydrocarbons` | alkane/alkene/alkyne only — any heteroatom group rejects |
| `Alkanes` | alkane only — C=C and C≡C also reject |

## 2. Routing

The molecule takes exactly one of three paths:

1. **Bridged** — exactly one chain-breaking heteroatom (C–X–C, X = O/N/S) → ester, ether,
   amine, amide, anhydride. Exception: an O bridge on a ring with no adjacent carbonyl
   (aryl ether) routes to the cyclic path instead.
2. **Cyclic** — a carbon ring is present.
3. **Chain** — everything else.

More than one bridging heteroatom → reject outright.

## 3. Chain path (acyclic, no bridge)

### Supported

- Straight and branched alkanes/alkenes/alkynes, **1–10 carbons in the main chain**.
- Correct chain selection: longest → most suffix groups → lowest suffix locants → lowest
  ene/yne locants → lowest prefix locants. Ties emit multiple names.
- Principal groups, in priority order (suffix when highest, prefix otherwise):

  | Group | Suffix | Prefix |
  |---|---|---|
  | `COOH` | -oic acid | carboxy |
  | `COCl` | -oyl chloride | chlorocarbonyl |
  | `CONH₂` | -amide | carbamoyl |
  | `C≡N` | -nitrile | cyano |
  | `CHO` (chain end) | -al | formyl |
  | `C=O` (mid-chain) | -one | oxo |
  | `–OH` | -ol | hydroxy |
  | `–SH` | -thiol | sulfanyl |
  | `–NH₂` | -amine | amino |
  | `=NH` | -imine | imino |

- Prefix-only groups: **fluoro, chloro, bromo, iodo, nitro**.
- Multiplying prefixes di–nona (up to 9 identical groups/substituents).

### Substituent grammar (one level only)

A branch is nameable if it is a carbon spine (attachment = locant 1) plus **either**
one functional group on the spine tip **or** methyl-only branches — never both:

- plain alkyl: `methyl`, `ethyl`, `propyl`, …
- one tip group: `(chloromethyl)`, `(hydroxymethyl)`, `(2-aminoethyl)`, …
- methyl branches: `(1-methylethyl)`, `(1,1-dimethylethyl)`, …

### Rejected

| Shape | Why |
|---|---|
| Main chain > 10 carbons | no root name beyond `dec-` |
| Branches on branches | not modelled |
| Groups on non-tip branch carbons | not modelled |
| Two or more groups on one branch tip (`–CF₃`, `–CHCl₂`) | not modelled |
| A branch tip carrying `COOH`/`COCl`/`CONH₂`/`C≡N` | prefix would double-count the tip carbon |
| Unsaturated substituents | not modelled |
| Cyclic substituents (cyclohexyl-, biphenyl-) | not modelled |
| ≥10 identical substituents | no multiplying prefix |

If *every* candidate chain carries an unnameable branch, the whole molecule is rejected.

## 4. Cyclic path

Exactly **one** all-carbon ring. Ring size follows the root names (up to 10).

### Non-aromatic rings

- `cyclopropane` … `cyclodecane`, with any G1-nameable substituents:
  alkyl, branched alkyl, halo, hydroxy, amino, nitro, alkoxy, `(hydroxymethyl)`, etc.
- Lowest-locant ring numbering with alphabetical tiebreak.
- **Rejected:** C=C or C≡C in the ring (cyclohexene), ring ketones (cyclohexanone —
  the C=O attaches by a double bond), fused/bridged/spiro systems, rings containing a
  heteroatom (THF, epoxides, pyrrolidine, lactones), two separate rings.

### Aromatic rings (benzene only)

- `benzene`, and mono/poly-substituted rings with prefix naming (`methylbenzene`,
  `1-chloro-2-methylbenzene`, `1,2,4,5-tetramethylbenzene`, `1,3-dinitrobenzene`).
- **Retained parents** (highest-priority substituent becomes the parent, locant 1):

  | Substituent | Parent |
  |---|---|
  | `–OH` | phenol |
  | `–NH₂` | phenylamine (never "aniline") |
  | `–COOH` | benzoic acid |
  | `–COCl` | benzoyl chloride |
  | `–CONH₂` | benzamide |
  | `–C≡N` | benzonitrile |
  | `–CHO` | benzaldehyde |

  e.g. `2-hydroxybenzoic acid`, `2,4,6-trinitrophenol`, `3-methylbenzaldehyde`,
  `4-nitrophenylamine`. A lower-priority parent group demotes to its prefix.
- **Chain parents — monosubstituted only:** `phenylmethanol` (–CH₂OH),
  `phenylethanone` (–COCH₃), `ethenylbenzene` (–CH=CH₂).
- **Aryl ethers** (routed here): `methoxybenzene`, `ethoxybenzene`, `2-methoxyphenol`,
  `4-methoxybenzaldehyde`. Alkoxy arms must be straight, saturated, all-carbon.

### Rejected on rings

| Shape | Example |
|---|---|
| Two substituents of the same principal class | benzene-1,2-diol, benzene dicarboxylic acid |
| Chain-parent substituent plus any other substituent | chloro-benzyl alcohol, 2-phenylethanol |
| Suffix-capable group on a ring substituent's tip | aminomethylbenzene, cyanomethylbenzene |
| Two groups on one substituent tip | dichloromethylbenzene |
| Unsaturated or cyclic ring substituents | ethynylbenzene, cyclohexylbenzene, vinylcyclohexane |
| Fused rings | naphthalene, 2,3-dihydrobenzofuran |
| Branched or substituted alkoxy arms | isopropoxybenzene, 2-chloroethoxybenzene |
| Two alkoxy groups (two bridges) | 1,2-dimethoxybenzene |
| Diaryl ether | diphenyl ether |

## 5. Bridged path (one C–X–C heteroatom)

Every side must be a **straight, saturated, unbranched, all-carbon chain attached at its
end**, or an unsubstituted **phenyl** ring, or the acid side. No other functional groups
anywhere on a side.

| Class | Bridge | Supported | Examples |
|---|---|---|---|
| Ether | O, 0 carbonyls | two alkyl sides | methoxymethane, methoxyethane, ethoxypropane |
| Ester | O, 1 carbonyl | alkyl or phenyl alcohol side; aliphatic or benzoate acid side | methyl ethanoate, ethyl propanoate, methyl benzoate, phenyl ethanoate, phenyl benzoate |
| Anhydride | O, 2 carbonyls | aliphatic only, symmetric or mixed | ethanoic anhydride, ethanoic propanoic anhydride |
| 2° amine | N, 0 carbonyls, 2 C | two alkyls, or one alkyl + phenyl | N-methylethanamine, N-methylphenylamine |
| 3° amine | N, 0 carbonyls, 3 C | **exactly one** phenyl + alkyls | N,N-dimethylphenylamine, N-ethyl-N-methylphenylamine |
| 2° amide | N, 1 carbonyl, 2 C | alkyl or phenyl on either side | N-methylethanamide, N-phenylethanamide, N-methylbenzamide |
| 3° amide | N, 1 carbonyl, 3 C | two N-substituents (alkyl/phenyl) | N,N-dimethylethanamide, N,N-diphenylethanamide, N,N-dimethylbenzamide |

### Substituted benzoate rings — esters only

An ester's benzoate ring may carry substituents nameable without group merging —
**nitro, halo, alkyl, hydroxy, amino, alkoxy**: `methyl 3-nitrobenzoate`,
`methyl 2-hydroxybenzoate`, `methyl 4-aminobenzoate`, `ethyl 4-methylbenzoate`.
Ring substituents needing merged groups (`–COOH`, `–CHO`, `–CONH₂`, `–COCl`) are rejected.
A substituted **benzamide** acid ring is rejected outright (`3-bromo-N-methylbenzamide`),
as is a substituted aromatic amine ring (`N,N-dimethyl-4-methylphenylamine`).

### Rejected on the bridged path

| Shape | Example |
|---|---|
| Branched or mid-chain-attached side | 2-methoxypropane |
| Unsaturated side | methyl vinyl ether |
| Extra functional group on a side | 2-methoxyethanol, aspirin, methyl nitroethanoate |
| Non-phenyl ring side | cyclohexyl ethanoate, benzyl ethanoate |
| Substituted phenyl side (non-acid) | 4-nitrophenyl ethanoate, paracetamol shape |
| S bridge | dimethyl sulfide |
| Aliphatic tertiary amine | trimethylamine |
| Diaryl / triaryl amine | diphenylamine, triphenylamine |
| Imides (N with 2 carbonyls) | diacetimide |
| Aromatic anhydride | benzoic anhydride |
| Two bridges | dimethoxymethane, 1,2-dimethoxybenzene |
| Ring closed through the bridge | THF, ethylene oxide, pyrrolidine, lactones |

## 6. Global rejections

- Two heteroatoms bonded to each other — peroxides, hydrazines, nitroso, nitrite and
  nitrate esters, N₂O₂ shapes. The nitro pattern is the only exemption.
- A heteroatom making a multiple bond to more than one carbon.
- More than one bridging heteroatom.
- Any atom over its valency; malformed bond indices.
- Charged species, isotopes, radicals, stereochemistry (E/Z, R/S), organometallics —
  not representable in the input model at all.
- Nitro under the `Hydrocarbons`/`Alkanes` specs (vocabulary gate).

## 7. Not modelled at all

Ring assemblies and fused polycyclics, heterocycles, spiro and bridged bicyclics,
multiplicative/di-suffix nomenclature (dioic acids, diols, diamines on rings), sulfonic
acids and other S/P oxyacids, ureas, carbamates, imines beyond `=NH`, azo/diazo,
alkyl halide chains longer than the substituent grammar allows, and any name requiring
enclosing marks beyond a single parenthesis level.
