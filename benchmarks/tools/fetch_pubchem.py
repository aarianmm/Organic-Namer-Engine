#!/usr/bin/env python3
"""Fetch a reference dataset of compounds from PubChem PUG-REST.

Fetches, per CID: IUPACName, a canonical SMILES (PubChem renamed the
CanonicalSMILES property to SMILES in 2025 -- this script tries SMILES
first and falls back to CanonicalSMILES for older API behaviour),
MolecularFormula, and Title.

Endpoint pattern:
    https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/<cids>/property/<props>/JSON

Rate limit: PubChem allows at most 5 requests/second; this script sleeps
between requests to stay under that.

Raw API responses are cached to a JSONL file (one JSON record per line, one
line per successfully-fetched CID) so re-runs with --resume don't re-hit
the API for CIDs already cached.

CLI:
    python fetch_pubchem.py --cids 1-6000 --out benchmarks/data/pubchem_raw.jsonl
    python fetch_pubchem.py --cid-file cids.txt --out benchmarks/data/pubchem_raw.jsonl --resume
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from typing import Iterable, Optional

import requests

BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids}/property/{props}/JSON"
PROPERTIES = ["IUPACName", "SMILES", "CanonicalSMILES", "MolecularFormula", "Title"]
USER_AGENT = "Organic-Namer-Engine-Benchmark/1.0 (research/benchmarking script; https://github.com/)"

BATCH_SIZE = 100
REQUESTS_PER_SECOND = 5
MIN_INTERVAL = 1.0 / REQUESTS_PER_SECOND


def parse_cid_range(spec: str) -> list[int]:
    """Parse a spec like '1-6000' or '1-100,250,300-310' into a sorted list
    of unique CIDs."""
    cids: set[int] = set()
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start_str, end_str = part.split("-", 1)
            start, end = int(start_str), int(end_str)
            cids.update(range(start, end + 1))
        else:
            cids.add(int(part))
    return sorted(cids)


def load_cid_file(path: str) -> list[int]:
    cids: list[int] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cids.append(int(line))
    return cids


def load_cached_cids(out_path: str) -> set[int]:
    """Return the set of CIDs already present in the raw cache file, so
    --resume can skip re-fetching them."""
    cached: set[int] = set()
    try:
        with open(out_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    record = json.loads(line)
                except json.JSONDecodeError:
                    continue
                cid = record.get("CID")
                if cid is not None:
                    cached.add(int(cid))
    except FileNotFoundError:
        pass
    return cached


def batched(items: list[int], size: int) -> Iterable[list[int]]:
    for i in range(0, len(items), size):
        yield items[i:i + size]


def fetch_batch(session: requests.Session, cids: list[int]) -> Optional[list[dict]]:
    """Fetch one batch of CIDs. Returns the list of property dicts from
    PubChem's response ("PropertyTable"."Properties"), or None on failure.
    """
    props = ",".join(PROPERTIES)
    url = BASE_URL.format(cids=",".join(str(c) for c in cids), props=props)
    try:
        resp = session.get(url, timeout=30)
    except requests.RequestException as exc:
        sys.stderr.write(f"request error for CIDs {cids[0]}-{cids[-1]}: {exc}\n")
        return None

    if resp.status_code == 404:
        # None of the CIDs in this batch exist / have properties -- not a
        # hard error, just an empty result.
        return []
    if not resp.ok:
        sys.stderr.write(
            f"HTTP {resp.status_code} for CIDs {cids[0]}-{cids[-1]}: {resp.text[:300]}\n"
        )
        return None

    try:
        data = resp.json()
    except ValueError as exc:
        sys.stderr.write(f"invalid JSON for CIDs {cids[0]}-{cids[-1]}: {exc}\n")
        return None

    return data.get("PropertyTable", {}).get("Properties", [])


def _get_smiles(record: dict) -> Optional[str]:
    """Return the SMILES for a PubChem property record, trying the 2025+
    'SMILES' property name first and falling back to the legacy
    'CanonicalSMILES' name."""
    return record.get("SMILES") or record.get("CanonicalSMILES")


def _is_usable(record: dict) -> bool:
    """A record is usable only if it has a non-empty IUPAC name and a
    non-empty SMILES (under either property name)."""
    return bool(record.get("IUPACName")) and bool(_get_smiles(record))


def run_fetch(cids: list[int], out_path: str, resume: bool) -> None:
    if resume:
        already = load_cached_cids(out_path)
        cids = [c for c in cids if c not in already]
        print(f"Resuming: {len(already)} CIDs already cached, {len(cids)} left to fetch")
        file_mode = "a"
    else:
        file_mode = "w"

    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})

    n_returned = 0
    n_written = 0
    n_skipped_unusable = 0
    n_batches_failed = 0
    last_request_time = 0.0

    with open(out_path, file_mode, encoding="utf-8") as fout:
        batches = list(batched(cids, BATCH_SIZE))
        for batch_no, batch in enumerate(batches, 1):
            elapsed = time.monotonic() - last_request_time
            if elapsed < MIN_INTERVAL:
                time.sleep(MIN_INTERVAL - elapsed)
            last_request_time = time.monotonic()

            records = fetch_batch(session, batch)
            if records is None:
                n_batches_failed += 1
                continue

            for record in records:
                n_returned += 1
                if not _is_usable(record):
                    n_skipped_unusable += 1
                    continue
                fout.write(json.dumps(record) + "\n")
                n_written += 1
            fout.flush()

            if batch_no % 10 == 0 or batch_no == len(batches):
                print(f"  batch {batch_no}/{len(batches)}: {n_written} usable records written so far "
                      f"({n_returned} returned, {n_skipped_unusable} skipped as unusable)")

    print(f"Done. PubChem returned {n_returned} property records; "
          f"wrote {n_written} usable records (skipped {n_skipped_unusable} with no "
          f"IUPACName/SMILES); {n_batches_failed} batches failed -> {out_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cids", help="CID range spec, e.g. '1-6000' or '1-100,250,300-310'")
    parser.add_argument("--cid-file", help="File with one CID per line")
    parser.add_argument("--out", required=True, help="Output JSONL cache file for raw records")
    parser.add_argument("--resume", action="store_true", help="Skip CIDs already present in --out")
    args = parser.parse_args()

    if not args.cids and not args.cid_file:
        parser.error("one of --cids or --cid-file is required")

    if args.cids:
        cids = parse_cid_range(args.cids)
    else:
        cids = load_cid_file(args.cid_file)

    print(f"Fetching {len(cids)} CIDs in batches of {BATCH_SIZE} "
          f"(rate limit {REQUESTS_PER_SECOND} req/s)")
    run_fetch(cids, args.out, args.resume)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
