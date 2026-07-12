# HTTP smoke tests

Ad-hoc smoke-test harness used during development of the cyclic/aromatic naming
work. Not a formal test suite (no CI integration, no xUnit) - just the
generate-payloads-and-POST-to-the-live-API setup, kept here instead of a
throwaway scratch directory so it survives between sessions.

## Usage

Start the API, then run the harness against it:

```bash
dotnet run --project src/OrganicNamer.Api --urls "http://localhost:5211" &
python3 tests/http_smoke/run.py
```

Run a single case:

```bash
python3 tests/http_smoke/run.py --case chlorobenzene
```

## Files

- `builders.py` - small composable helpers for constructing atom/bond JSON
  payloads (rings, chains, substituents) in the project's explicit-hydrogen
  convention.
- `cases.py` - named cases (`{"payload": ..., "expect": ...}`), covering the
  chain-based, cyclic, branch-tip-group, and aromatic naming paths, plus
  cases that must be rejected rather than silently misnamed, plus a few
  adversarial/malformed inputs that must not crash the server.
- `run.py` - POSTs each case to the running API and reports pass/fail.

Cases track whatever's actually implemented on the current branch - update
`cases.py` as new molecule classes land (see `Algorithm-Extension-Plan.md`
at the repo's parent directory for the phase-by-phase design).
