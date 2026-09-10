#!/usr/bin/env python3
"""Name-comparison scorer for the pubchem_named.jsonl benchmark.

Strict exact-string matching against PubChem's reference IUPACName
undercounts the engine badly, because most of the gap is pure naming
*convention* -- both names denote the same molecule, they just follow
different (both legitimate) IUPAC styling choices. This script applies a
conservative, table-driven, purely lexical normalisation layer on top of
strict matching and reports both scores.

A molecule counts as a match if ANY of the engine's emitted `names` equals
the reference (strict), or equals the reference after normalisation
(normalised). Normalisation is applied symmetrically to both the reference
name and every candidate engine name -- never structurally, never fuzzily,
never by edit distance. Every rule below must be (and is) justifiable to a
chemist in one sentence.

IMPORTANT: this script hard-codes a regression check. 18 CIDs in the
benchmark are genuine structural mismatches (the engine names a different,
larger molecule than PubChem's reference -- see "formyl/carbamoyl" rules
below) and must NEVER be turned into normalised matches. The script asserts
this after scoring; if a future rule change causes any of them to start
matching, that rule is too aggressive and the assertion will fail the run.

===========================================================================
NORMALISATION RULES
===========================================================================

Rule A -- terminal_group_locant_elision
    A locant is dropped when it is the ONLY chemically possible position.
    Carboxylic acid, aldehyde, amide, nitrile and acyl-halide carbons are
    always the terminal (or, for the "di" forms, both terminal) carbon(s)
    of their chain by definition -- there is no other place they could be
    -- so IUPAC omits the "-1-" (or "-1,N-") locant regardless of chain
    length: butan-1-al -> butanal, pentan-1-oic acid -> pentanoic acid,
    hexane-1,6-dioic acid -> hexanedioic acid, pentane-1,5-dial ->
    pentanedial, ethan-1-amide -> ethanamide, propane-1-nitrile ->
    propanenitrile, propan-1-oyl chloride -> propanoyl chloride.
    Verified against the data: this single rule accounts for the large
    majority (163/187) of the newly-recovered matches, essentially every
    "-an-1-oic acid" / "-ane-1,N-dioic acid" / "-an-1-al" pair in the set.

Rule B -- short_chain_suffix_locant_elision
    Unlike the groups in Rule A, -ol and -thiol are NOT inherently
    terminal (propan-2-ol is a real, different compound from propan-1-ol),
    so their locant is normally required. The exception is methane and
    ethane: with only 1-2 carbons the "1" conveys no information the
    unlocanted name doesn't already convey (verified against the data --
    PubChem keeps the locant for propan-1-ol, octan-1-ol, etc., but drops
    it only for ethan-1-ol -> ethanol, 2-chloroethan-1-ol ->
    2-chloroethanol, methan-1-ol -> methanol), so it is dropped only for
    those two stems.

Rule C -- methane_locant_elision
    Methane has exactly one carbon, so EVERY locant attached to it (on
    substituents or on the suffix) is redundant, not just a suffix "-1-":
    1,1,1,1-tetrachloromethane -> tetrachloromethane, 1,1-dibromomethane
    -> dibromomethane, methane-1,1-dial -> methanedial, 1-aminomethan-1-
    oic acid -> aminomethanoic acid. Applied only when the name contains
    "methan" and no digit other than 1 appears anywhere in it (a name with
    a 2+ digit cannot be a plain methane derivative), so it can never fire
    on a longer, genuinely locant-bearing chain.

Rule D -- multiplying_prefix_vowel_elision
    The final "a" of a multiplying prefix (tetra-, penta-, hexa-, ...) is
    elided before the vowel-initial suffix "-ol", a standard orthographic
    contraction: tetraol -> tetrol, pentaol -> pentol, hexaol -> hexol.
    Verified in the data for tetra/penta/hexa (e.g. hexane-1,2,3,4,5,6-
    hexaol vs PubChem's hexane-1,2,3,4,5,6-hexol); hepta/octa/nona/deca are
    included for completeness of the same orthographic rule even though no
    example of that chain length occurs in this dataset.

Rule "bis_di_sulfanyl_convention"
    "bis(sulfanyl)" and "disulfanyl" name the same pair of separate -SH
    substituents when both locants are given explicitly (6,8-
    bis(sulfanyl)octanoic acid == 6,8-disulfanyloctanoic acid); PubChem's
    generator prefers the bis(...) form to avoid any visual confusion with
    a -S-S- disulfane chain. Purely orthographic, symmetric substitution.

Rule "branched_alkyl_convention"
    "1-methylethyl" (pre-2013 IUPAC) and "propan-2-yl" (2013 IUPAC) name
    the identical isopropyl substituent under different recommendation
    editions; "(1,1-dimethylethyl)" / "tert-butyl" are the equivalent pair
    for tert-butyl (included for the same reason though not present in
    this dataset). A locant-prefixed substituent no longer needs enclosing
    parentheses once it no longer reads as ambiguous (4-(1-methylethyl)...
    -> 4-propan-2-yl...), but a multiplying prefix (di/tri/bis) still
    needs them (2,6-di(1-methylethyl)... -> 2,6-di(propan-2-yl)...) --
    both forms are observed verbatim in the data.

Rule "retained_trivial_name_alias" (Class 2 alias table)
    PubChem prefers a set of IUPAC-retained trivial names over the fully
    systematic substitutive name; the engine only produces the systematic
    form. Every pair below is a name PubChem actually returned as
    `reference_name` for at least one row, paired with the systematic name
    the engine actually produced for that same row (verified by manual
    structure check, not just string shape) -- see ALIAS_TABLE. The
    substitution is substring-level (with a guard, see code comment) so
    that substituted derivatives also match: 2-chloroethanoic acid ->
    2-chloroacetic acid, N-phenylethanamide -> N-phenylacetamide.

Rule "amine_prefix_suffix_convention"
    IUPAC 2013 requires the amine group, when it is the principal
    characteristic group, to be cited as the suffix "-amine"/"-diamine"
    rather than the prefix "amino-"; PubChem follows this, the engine uses
    the prefix form. Both denote the identical molecule. Because moving a
    suffix to a prefix reorders the whole name (not just a substring), this
    is implemented as four explicit whole-name pairs, one per case that
    actually occurs in the data (see AMINE_WHOLE_MAP) -- not a generative
    rule, so it cannot misfire on an unseen ring size or substitution
    pattern.

===========================================================================
ALIAS TABLE (Class 2 retained/trivial names) -- with justification
===========================================================================
  ethanoic acid          <-> acetic acid        PubChem's preferred retained name for CH3COOH
  methanoic acid         <-> formic acid        retained name for HCOOH
  ethanal                <-> acetaldehyde       retained name for CH3CHO
  methanal               <-> formaldehyde       retained name for HCHO
  ethanamide             <-> acetamide          retained name for CH3CONH2
  methanamide            <-> formamide          retained name for HCONH2
  phenylamine            <-> aniline            retained name for C6H5NH2
  methylbenzene          <-> toluene            retained name for C6H5CH3
  ethanedioic acid       <-> oxalic acid        retained name for HOOC-COOH
  aminomethanoic acid    <-> carbamic acid      retained name for H2N-COOH
  hydroxymethanoic acid  <-> carbonic acid      retained name for HO-COOH
  methanedial            <-> carbon dioxide     the "dialdehyde of methane" is O=C=O
  methanenitrile         <-> formonitrile       retained name for H-C#N (HCN)
  hydroxymethanenitrile  <-> cyanic acid        HO-C#N
  sulfanylmethanenitrile <-> thiocyanic acid    HS-C#N
  aminomethanamide       <-> urea               H2N-CO-NH2
  iminomethanediamine    <-> guanidine          HN=C(NH2)2
  tribromomethane        <-> bromoform          retained name for CHBr3

Two names PubChem returns that visually resemble candidates for this table
were deliberately NOT added: "oxamic acid" (H2N-CO-COOH) and "oxaldehydic
acid" (HOOC-CHO). For both, the engine's actual output uses a "carbamoyl"/
"formyl" substituent that adds an extra carbon not present in the true
molecule (see the formyl/carbamoyl note below) -- i.e. the engine's name for
those two rows is a genuine structural error, not a styling difference, and
these are 2 of the 18 hard-mismatch CIDs. Aliasing them would have silently
defeated the regression check, so they are intentionally absent.

===========================================================================
THE GENUINE MISMATCHES -- NOT NORMALISED, BY DESIGN
===========================================================================
This section used to list 18 CIDs (29, 30, 48, 207, 236, 296, 523, 738,
742, 760, 829, 844, 868, 974, 1112, 1122, 2859, 5961). Seventeen of them
shared one cause: the engine named an on-chain -CHO or -CONH2 group with a
"formyl"/"carbamoyl" SUBSTITUENT prefix, which brings its own carbon atom,
while PubChem's reference put that carbon inside the main chain via "oxo"
(no extra carbon). E.g. CID 29: reference "2-amino-3-oxopropanoic acid" is
OHC-CH(NH2)-COOH (3 carbons); the engine's "2-amino-3-formylpropan-1-oic
acid" described propanoic acid (3 carbons) PLUS the formyl carbon = 4 --
a genuinely different, larger molecule.

That defect has since been fixed in the engine, and the 17 split in two:

  * FIXED_BY_OXO_RENAME (11 CIDs) + FIXED_BUT_RETAINED_NAME (1 CID, 760,
    glyoxylic acid, where PubChem uses a retained trivial name this file
    does not alias) -- an on-chain carbonyl now takes the
    prefix "oxo", so the engine emits the same structure PubChem does and
    these are now correct. They still need Rule A (terminal-group locant
    elision) to match, since the engine writes "propan-1-oic acid" where
    PubChem writes "propanoic acid" -- so they are NORMALISED matches, not
    strict ones. `_assert_genuine_mismatches` now asserts they DO match:
    if one ever stops matching, the fix has regressed.

  * REFUSED_AFTER_FIX (5 CIDs) -- a non-principal COOH/COCl/CONH2/C=N on
    the main chain cannot be named as a substituent at all (its carbon
    does not belong to the parent chain), so the engine now REFUSES these
    rather than emitting a wrong name. They are absent from the scorable
    rows entirely, and the guard asserts that absence.

HARD_MISMATCH_CIDS is therefore now empty. The 18th CID, 2859, was a
separate defect: a bicyclic epoxide-tetrol whose bridging oxygen was
misread by ring-system routing -- the engine treated the epoxide like an
ordinary ether substituent instead of noticing its two attachment carbons
were already part of the carbocycle, and emitted the name of an unrelated
hexahydroxycyclohexane, silently dropping the epoxide. That routing check
has since been fixed to reject a bridge whose own two carbons already
close a ring, so CID 2859 is now REFUSED rather than misnamed and has
moved into REFUSED_AFTER_FIX (6 CIDs). The engine has zero known genuine
structural mismatches remaining on this corpus; this is asserted in code
below.

===========================================================================
KNOWN GAPS -- LEFT UNNORMALISED ON PURPOSE (not lexically fixable / not a
convention difference)
===========================================================================
* Stereodescriptors: rows like "(2S)-2-aminopropanoic acid" vs the engine's
  "2-aminopropan-1-oic acid" differ by more than convention -- the engine
  does not emit stereodescriptors at all. Stripping "(2S)-" etc. to force a
  match would hide a real capability gap, not a styling choice, so it is
  not done. (CIDs 5780, 5851, 5852, 5862, 5950, 5951, 5960, 5962, 5984;
  5961 is also stereo-affected but already counted in the hard-18 above.)
* Chain/substituent selection: rows like CID 87 ("3-hydroxy-2-
  methylpropanoic acid" vs the engine's "2-(hydroxymethyl)propan-1-oic
  acid") describe the SAME molecule -- the engine just chose the shorter
  branch as the main chain instead of the longer one, which is a violation
  of IUPAC's "longest chain through the principal group" rule, not a
  matter of two equally-valid conventions. Fixing this needs a structural
  chain-reselection algorithm, which is explicitly out of scope for a
  lexical scorer. (CIDs 38, 87, 848, 1531, 3009.)
* CID 892: "cyclohexane-1,2,3,4,5,6-hexol" vs the engine's "1,2,3,4,5,6-
  hexahydroxycyclohexane" -- the engine expresses the principal
  characteristic group (-OH) entirely as prefixes rather than the -ol
  suffix. That is a structural reordering of the whole name, not a
  substring swap, so it is left unfixed (distinct from CID 2859's genuine
  structural error above, though the engine happens to reuse the same
  wrong string for both).
* CID 1254: mixes the chain-selection issue above with the same suffix-vs-
  prefix issue as CID 892.
* CID 3094: "N-octyloctan-1-amine" (reference, locanted) vs the engine's
  "N-octyloctanamine" (no locant). This is the reverse of Rule B: the
  engine drops a locant that IS chemically required (octan-2-amine is a
  different molecule), which looks like an engine inconsistency rather
  than a legitimate alternate convention, so no rule was added to
  "restore" the locant.

Usage:
    python score_names.py --in benchmarks/results/pubchem_named.jsonl \\
        --out benchmarks/results/name_scores.jsonl
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from typing import Callable, Optional

# ---------------------------------------------------------------------------
# CIDs independently confirmed (via OPSIN) to be genuine structural mismatches
# that the normalised scorer must NEVER match. See the module docstring
# section "THE GENUINE MISMATCHES" above for how this list got down to zero.
# Kept as a set (rather than removed) so a future regression that reintroduces
# a genuine mismatch has somewhere to go.
# ---------------------------------------------------------------------------
HARD_MISMATCH_CIDS: set[str] = set()

# Previously hard mismatches, fixed by the formyl -> oxo rename. These must now
# MATCH (after Rule A locant elision). Asserted positively: a regression in the
# engine's prefix handling would show up here first.
FIXED_BY_OXO_RENAME = {
    "29", "30", "207", "296", "523", "742", "829", "844", "868",
    "1112", "1122",
}

# Also fixed by the oxo rename and confirmed structurally correct by the OPSIN
# round-trip, but still not a string match: PubChem cites a retained trivial
# name this file's alias table does not carry. CID 760 is glyoxylic acid --
# engine "2-oxoethan-1-oic acid" vs PubChem "oxaldehydic acid". Asserted only
# to be scorable (i.e. named, not refused); NOT asserted to match, and
# deliberately not aliased, so the headline normalised percentage stays a
# measurement rather than something tuned to satisfy this guard.
FIXED_BUT_RETAINED_NAME = {
    "760",
}

# Previously hard mismatches, now refused outright: a non-principal
# COOH/COCl/CONH2/C=N on the main chain has no valid substituent name. These
# must be ABSENT from the scorable rows (status != "named"). CID 2859 is the
# odd one out here -- refused via the separate ring-routing fix (bicyclic
# epoxide-tetrol), not the formyl/carbamoyl fix that produced the other five.
REFUSED_AFTER_FIX = {
    "48", "236", "738", "974", "2859", "5961",
}


# ---------------------------------------------------------------------------
# Rule A: terminal-group locant elision (any chain length -- the group is
# inherently terminal, so "-1-" / "-1,N-" carries no information).
# ---------------------------------------------------------------------------
_RULE_A_PATTERNS = [
    (re.compile(r'-1,\d+-dioic acid'), 'dioic acid'),
    (re.compile(r'-1-oic acid'), 'oic acid'),
    (re.compile(r'-1,\d+-dial'), 'dial'),
    (re.compile(r'-1-al'), 'al'),
    (re.compile(r'-1-amide'), 'amide'),
    (re.compile(r'-1-nitrile'), 'nitrile'),
    (re.compile(r'-1-oyl (chloride|bromide|fluoride|iodide)'), r'oyl \1'),
]


def _rule_terminal_group_locant_elision(name: str) -> str:
    out = name
    for pat, repl in _RULE_A_PATTERNS:
        out = pat.sub(repl, out)
    return out


# ---------------------------------------------------------------------------
# Rule B: -ol / -thiol locant elision, restricted to methane/ethane stems.
# ---------------------------------------------------------------------------
_RULE_B_PATTERN = re.compile(r'(meth|eth)an-1-(ol|thiol)')


def _rule_short_chain_suffix_locant_elision(name: str) -> str:
    return _RULE_B_PATTERN.sub(r'\1an\2', name)


# ---------------------------------------------------------------------------
# Rule C: methane parent -- every locant on it is redundant (one carbon).
# Guarded so it can only fire on a genuine methane derivative: the name
# must mention "methan" and must not contain any digit other than 1
# (a real multi-carbon locant set would include a 2, 3, ...).
# ---------------------------------------------------------------------------
_RULE_C_OTHER_DIGIT = re.compile(r'[2-9]')
_RULE_C_MID = re.compile(r'-1(,1)*-')
_RULE_C_LEAD = re.compile(r'^1(,1)*-')


def _rule_methane_locant_elision(name: str) -> str:
    if 'methan' not in name or _RULE_C_OTHER_DIGIT.search(name):
        return name
    out = _RULE_C_MID.sub('', name)
    out = _RULE_C_LEAD.sub('', out)
    return out


# ---------------------------------------------------------------------------
# Rule D: elide the final "a" of a multiplying prefix before -ol.
# ---------------------------------------------------------------------------
_RULE_D_MAP = [
    ('tetraol', 'tetrol'), ('pentaol', 'pentol'), ('hexaol', 'hexol'),
    ('heptaol', 'heptol'), ('octaol', 'octol'), ('nonaol', 'nonol'),
    ('decaol', 'decol'),
]


def _rule_multiplying_prefix_vowel_elision(name: str) -> str:
    out = name
    for a, b in _RULE_D_MAP:
        out = out.replace(a, b)
    return out


# ---------------------------------------------------------------------------
# bis(sulfanyl) <-> disulfanyl.
# ---------------------------------------------------------------------------
def _rule_bis_di_sulfanyl_convention(name: str) -> str:
    return name.replace('bis(sulfanyl)', 'disulfanyl')


# ---------------------------------------------------------------------------
# Branched-alkyl substituent naming: pre-2013 vs 2013 IUPAC recommendations.
# ---------------------------------------------------------------------------
_RULE_C3_MAP = [
    ('1-methylethyl', 'propan-2-yl'),
    ('1,1-dimethylethyl', 'tert-butyl'),
]
_RULE_C3_UNPAREN_PROPAN2YL = re.compile(r'(\d)-\(propan-2-yl\)')
_RULE_C3_UNPAREN_TERTBUTYL = re.compile(r'(\d)-\(tert-butyl\)')


def _rule_branched_alkyl_convention(name: str) -> str:
    out = name
    for a, b in _RULE_C3_MAP:
        out = out.replace(a, b)
    # A single locant-prefixed substituent no longer needs parentheses once
    # translated; a multiplying prefix (di/tri/bis) still does, and is left
    # untouched.
    out = _RULE_C3_UNPAREN_PROPAN2YL.sub(r'\1-propan-2-yl', out)
    out = _RULE_C3_UNPAREN_TERTBUTYL.sub(r'\1-tert-butyl', out)
    return out


# ---------------------------------------------------------------------------
# Class 2: retained/trivial name <-> fully systematic name.
# See the module docstring's ALIAS TABLE section for the one-line
# justification of every entry.
# ---------------------------------------------------------------------------
ALIAS_TABLE: list[tuple[str, str]] = [
    ('ethanoic acid', 'acetic acid'),
    ('methanoic acid', 'formic acid'),
    ('ethanal', 'acetaldehyde'),
    ('methanal', 'formaldehyde'),
    ('ethanamide', 'acetamide'),
    ('methanamide', 'formamide'),
    ('phenylamine', 'aniline'),
    ('methylbenzene', 'toluene'),
    ('ethanedioic acid', 'oxalic acid'),
    ('aminomethanoic acid', 'carbamic acid'),
    ('hydroxymethanoic acid', 'carbonic acid'),
    ('methanedial', 'carbon dioxide'),
    ('methanenitrile', 'formonitrile'),
    ('hydroxymethanenitrile', 'cyanic acid'),
    ('sulfanylmethanenitrile', 'thiocyanic acid'),
    ('aminomethanamide', 'urea'),
    ('iminomethanediamine', 'guanidine'),
    ('tribromomethane', 'bromoform'),
]
# Longest pattern first, so e.g. "aminomethanoic acid" -> "carbamic acid"
# is tried before the more generic "methanoic acid" -> "formic acid" gets a
# chance to fire on the same substring.
_ALIAS_TABLE_SORTED = sorted(ALIAS_TABLE, key=lambda pair: -len(pair[0]))


def _rule_retained_trivial_name_alias(name: str) -> str:
    out = name
    for systematic, trivial in _ALIAS_TABLE_SORTED:
        if systematic.startswith('ethan'):
            # "methan..." always contains "ethan..." as a tail substring
            # (methan|e = m + ethane); guard so the ethane-rooted alias
            # can't corrupt a methane-rooted name (methane has its own
            # entries in this table).
            out = re.sub(r'(?<!m)' + re.escape(systematic), trivial, out)
        else:
            out = out.replace(systematic, trivial)
    return out


# ---------------------------------------------------------------------------
# Amine principal-group suffix vs. prefix convention. Whole-name pairs only
# (see module docstring) -- every entry is a full string that actually
# appears in either `reference_name` or `names` in this dataset.
# ---------------------------------------------------------------------------
AMINE_WHOLE_MAP: dict[str, str] = {
    'aminocyclopentane': 'cyclopentanamine',
    'aminocycloheptane': 'cycloheptanamine',
    'aminocyclooctane': 'cyclooctanamine',
    '1,2-diaminocyclohexane': 'cyclohexane-1,2-diamine',
}


def _rule_amine_prefix_suffix_convention(name: str) -> str:
    return AMINE_WHOLE_MAP.get(name, name)


# Rules are applied in this fixed order. Order matters for two reasons:
# (1) Rule A/B/C must run before the alias table, since e.g. "ethan-1-oic
#     acid" must become "ethanoic acid" before it can be recognised as an
#     alias of "acetic acid"; (2) the alias table's own internal ordering
#     (longest-pattern-first) is handled inside that rule.
RULES: list[tuple[str, Callable[[str], str]]] = [
    ('terminal_group_locant_elision', _rule_terminal_group_locant_elision),
    ('short_chain_suffix_locant_elision', _rule_short_chain_suffix_locant_elision),
    ('methane_locant_elision', _rule_methane_locant_elision),
    ('multiplying_prefix_vowel_elision', _rule_multiplying_prefix_vowel_elision),
    ('bis_di_sulfanyl_convention', _rule_bis_di_sulfanyl_convention),
    ('branched_alkyl_convention', _rule_branched_alkyl_convention),
    ('retained_trivial_name_alias', _rule_retained_trivial_name_alias),
    ('amine_prefix_suffix_convention', _rule_amine_prefix_suffix_convention),
]


def normalize_tracked(name: str) -> tuple[str, list[str]]:
    """Apply every rule in sequence, returning (normalised_name, rules_that_fired).

    A rule is recorded as having "fired" on this string only if it actually
    changed the string at the point it ran (so e.g. a rule that would be a
    no-op after an earlier rule already did the same job isn't credited
    twice).
    """
    out = name
    fired: list[str] = []
    for rule_name, fn in RULES:
        new_out = fn(out)
        if new_out != out:
            fired.append(rule_name)
        out = new_out
    return out, fired


def normalize(name: str) -> str:
    return normalize_tracked(name)[0]


def score_row(row: dict) -> dict:
    """Score a single scorable ('status' == 'named') benchmark row."""
    ref = row['reference_name']
    names = row.get('names') or []

    strict_match = ref in names

    norm_ref, ref_fired = normalize_tracked(ref)
    normalised_match = False
    rule_applied: Optional[str] = None

    if strict_match:
        normalised_match = True
    else:
        for candidate in names:
            norm_candidate, candidate_fired = normalize_tracked(candidate)
            if norm_candidate == norm_ref:
                normalised_match = True
                combined = sorted(set(ref_fired) | set(candidate_fired))
                rule_applied = ','.join(combined) if combined else None
                break

    return {
        'id': row['id'],
        'reference_name': ref,
        'names': names,
        'strict_match': strict_match,
        'normalised_match': normalised_match,
        'rule_applied': rule_applied,
    }


def _assert_genuine_mismatches(scored_by_id: dict[str, dict]) -> None:
    """Three-way regression guard over the formerly-18 hard-mismatch CIDs.

    See the module docstring section "THE GENUINE MISMATCHES" for the full
    story. Asserts, in order:

      1. HARD_MISMATCH_CIDS are still scorable and still do NOT match. A
         failure here means a normalisation rule became too aggressive and
         needs narrowing -- do not loosen the assertion.
      2. FIXED_BY_OXO_RENAME are scorable and DO normalised-match, and
         FIXED_BUT_RETAINED_NAME are at least scorable. A failure here means
         the formyl -> oxo fix has regressed in the engine.
      3. REFUSED_AFTER_FIX are absent from the scorable rows entirely. A
         failure here means the engine started naming a molecule it should
         refuse -- i.e. it is emitting a carbon-subsuming prefix again.
    """
    missing = []
    violations = []
    for cid in sorted(HARD_MISMATCH_CIDS, key=int):
        rec = scored_by_id.get(cid)
        if rec is None:
            missing.append(cid)
            continue
        if rec['normalised_match']:
            violations.append(cid)
    if missing:
        raise AssertionError(
            f"Expected hard-mismatch CIDs missing from scorable rows: {missing}"
        )
    if violations:
        raise AssertionError(
            "Normalisation incorrectly matched known-genuine structural "
            f"mismatches (CIDs {violations}). A rule is too aggressive -- "
            "see the module docstring's 'THE GENUINE MISMATCHES' section."
        )

    absent = sorted(
        (cid for cid in FIXED_BUT_RETAINED_NAME if cid not in scored_by_id),
        key=int,
    )
    if absent:
        raise AssertionError(
            f"CIDs {absent} should be named after the formyl -> oxo fix but "
            "are not scorable at all -- the engine is refusing a molecule it "
            "should name."
        )

    unfixed = []
    for cid in sorted(FIXED_BY_OXO_RENAME, key=int):
        rec = scored_by_id.get(cid)
        if rec is None or not rec['normalised_match']:
            unfixed.append(cid)
    if unfixed:
        raise AssertionError(
            f"CIDs {unfixed} should match after the formyl -> oxo fix but do "
            "not. Either the engine regressed to emitting a carbon-subsuming "
            "prefix, or Rule A (terminal-group locant elision) stopped firing."
        )

    resurrected = sorted(
        (cid for cid in REFUSED_AFTER_FIX if cid in scored_by_id), key=int
    )
    if resurrected:
        raise AssertionError(
            f"CIDs {resurrected} are scorable but should be refused: a "
            "non-principal COOH/COCl/CONH2/C=N on the main chain has no "
            "valid substituent name, so the engine must reject it."
        )


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--in', dest='in_path', required=True,
                         help='Path to pubchem_named.jsonl (or equivalent)')
    parser.add_argument('--out', dest='out_path', required=True,
                         help='Path to write per-molecule scoring JSONL')
    args = parser.parse_args(argv)

    scored: list[dict] = []
    with open(args.in_path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            if row.get('status') != 'named':
                continue
            scored.append(score_row(row))

    total = len(scored)
    strict_count = sum(1 for r in scored if r['strict_match'])
    normalised_count = sum(1 for r in scored if r['normalised_match'])

    rule_credit: dict[str, int] = {name: 0 for name, _ in RULES}
    for r in scored:
        if r['rule_applied']:
            for rule_name in r['rule_applied'].split(','):
                rule_credit[rule_name] += 1

    scored_by_id = {r['id']: r for r in scored}
    _assert_genuine_mismatches(scored_by_id)

    with open(args.out_path, 'w', encoding='utf-8') as f:
        for r in scored:
            f.write(json.dumps(r, ensure_ascii=False) + '\n')

    def pct(n: int) -> str:
        return f"{(n / total * 100):.2f}%" if total else "n/a"

    print(f"Total scorable (status == 'named'): {total}")
    print(f"Strict exact matches:     {strict_count} / {total} ({pct(strict_count)})")
    print(f"Normalised matches:       {normalised_count} / {total} ({pct(normalised_count)})")
    print()
    print("Per-rule breakdown (count of newly-matched rows this rule "
          "contributed to; a row needing several rules together counts "
          "toward each of them):")
    for rule_name, _ in RULES:
        print(f"  {rule_name}: {rule_credit[rule_name]}")
    print()
    print(f"Regression check passed: {len(HARD_MISMATCH_CIDS)} known-genuine "
          f"mismatch(es) still unmatched, {len(FIXED_BY_OXO_RENAME)} "
          f"oxo-rename fixes still matching (+{len(FIXED_BUT_RETAINED_NAME)} "
          f"fixed but retained-name), {len(REFUSED_AFTER_FIX)} "
          "formerly-misnamed molecules still refused.")
    print(f"Wrote per-molecule results to {args.out_path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
