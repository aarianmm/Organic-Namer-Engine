# Organic-Namer-Engine
A tool for computing the IUPAC name from the structure of an organic chemical

## Input conventions

**Nitro groups (–NO₂):** the JSON schema has no formal-charge field, so
nitro is accepted only in the hypervalent neutral form — one N with a single
bond to carbon and two double bonds to terminal oxygens
(`N(=O)(=O)`). The charge-separated form N⁺(–O⁻)(=O) is unrepresentable and
is rejected. Implicit-hydrogen filling never adds H to a nitro N.
