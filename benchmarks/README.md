# benchmarks

Reproducible measurement harness for the naming engine. Results and methodology
are in [`BENCHMARKS.md`](BENCHMARKS.md).

```
benchmarks/
  BENCHMARKS.md                    the report
  runner/OrganicNamer.Benchmarks/  C# harness (references Core directly, no HTTP)
  tools/                           Python analysis + dataset scripts
  data/                            corpora (generated and fetched)
  results/                         measurement outputs
```

## Setup

The C# runner needs only the .NET SDK. The Python tools need a virtualenv on
Python 3.13 (RDKit has no 3.14 wheel yet):

```bash
/opt/homebrew/bin/python3.13 -m venv tools/.venv
tools/.venv/bin/pip install rdkit requests
```

OPSIN is vendored as `tools/opsin-2.9.0.jar` and needs a JRE on `PATH`.

Run the Python tools with `tools/.venv/bin/python`. Paths below are relative to
the repository root.

## The runner

```bash
# name every molecule in a JSONL corpus
dotnet run -c Release --project benchmarks/runner/OrganicNamer.Benchmarks -- \
    name --in <in.jsonl> --out <out.jsonl>

# throughput + latency percentiles
dotnet run -c Release --project benchmarks/runner/OrganicNamer.Benchmarks -- \
    bench --in <in.jsonl> --out <results.json> [--warmup 3] [--reps 10]
```

Input is one molecule per line: `{"id": ..., "atoms": [...]}` in the engine's
JSON format (hydrogens may be omitted). Any other fields are passed through to
the output, so downstream scoring keeps its metadata. Output adds `status`
(`named` / `rejected` / `invalid`), `names`, `micros`, and `error`.

## The tools

| Script | Purpose |
|---|---|
| `combinatorics.py` | Analytic count of the nameable molecule space. Pure enumeration — never runs the engine. |
| `alkane_isomers.py` | Builds every alkane isomer up to C*n*, asserted against OEIS A000602. Ground truth for coverage. |
| `generate_corpus.py` | Samples molecules from the supported grammar. Feeds the round-trip test and validates the combinatorics scope model. |
| `smiles_to_graph.py` | SMILES ↔ engine JSON graph. `--selftest` round-trips 51 reference molecules. |
| `fetch_pubchem.py` | Fetches reference names + SMILES from PubChem PUG-REST. Rate-limited, cached, `--resume`-able. |
| `opsin_roundtrip.py` | Structural accuracy: name → OPSIN → structure, compared to the input graph. |
| `score_names.py` | String agreement with PubChem, strict and convention-normalised. |

## End to end

```bash
cd benchmarks
V=tools/.venv/bin/python

# --- scope ---
$V tools/combinatorics.py --strata-max 14 --out results/combinatorics.json
$V tools/alkane_isomers.py --max-carbons 12 --out data/alkane_isomers.jsonl

# --- corpora ---
$V tools/generate_corpus.py --count 50000 --seed 7 --out data/generated_corpus.jsonl
$V tools/fetch_pubchem.py --cids 1-6000 --out data/pubchem_raw.jsonl
$V tools/smiles_to_graph.py --in data/pubchem_raw.jsonl \
    --out data/pubchem_graphs.jsonl --fail-out data/pubchem_unconvertible.jsonl

# --- name everything (from the repo root) ---
cd ..
for c in generated_corpus alkane_isomers; do
  dotnet run -c Release --project benchmarks/runner/OrganicNamer.Benchmarks -- \
      name --in benchmarks/data/$c.jsonl --out benchmarks/results/${c%_corpus}_named.jsonl
done

# --- accuracy + performance ---
cd benchmarks
$V tools/opsin_roundtrip.py --in results/generated_named.jsonl \
    --out results/generated_roundtrip.jsonl
$V tools/score_names.py --in results/pubchem_named.jsonl \
    --out results/name_scores.jsonl
cd .. && dotnet run -c Release --project benchmarks/runner/OrganicNamer.Benchmarks -- \
    bench --in benchmarks/data/generated_corpus.jsonl --out benchmarks/results/throughput.json
```

Seeds are fixed, so corpora regenerate identically. The PubChem step is the only
one needing network access.
