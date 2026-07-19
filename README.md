# Organic-Namer-Engine
A tool for computing the IUPAC name from the structure of an organic chemical

## Input conventions

**Nitro groups (–NO₂):** the JSON schema has no formal-charge field, so
nitro is accepted only in the hypervalent neutral form — one N with a single
bond to carbon and two double bonds to terminal oxygens
(`N(=O)(=O)`). The charge-separated form N⁺(–O⁻)(=O) is unrepresentable and
is rejected. Implicit-hydrogen filling never adds H to a nitro N.

**Substituted benzoate rings (esters only):** a benzoate ring side
(`methyl benzoate`, `phenyl benzoate`, ...) may carry ring substituents that
G1 can name pre-merge — nitro, halo, alkyl, hydroxy, amino, alkoxy — e.g.
`methyl 3-nitrobenzoate`, `methyl 2-hydroxybenzoate`. A ring substituent
whose naming needs merged groups (`-COOH`, `-CHO`, `-CONH2`, `-COCl`) is
rejected rather than misnamed. This vocabulary is **ester-only**: a
substituted benzamide acid ring (`3-bromo-N-methylbenzamide`) is rejected,
because the `N-` prefix scaffold can't interleave locant kinds with a
ring-numeric prefix in the correct alphabetical order. Unsubstituted
benzamide acid rings (`N-methylbenzamide`) are unaffected.
