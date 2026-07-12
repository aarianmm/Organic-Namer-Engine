"""
HTTP smoke harness for OrganicNamer.Api.

Not a formal test suite - just the ad-hoc POST-to-the-live-API setup used
during development, kept here so it isn't lost between sessions. Requires
the API to already be running (dotnet run --project src/OrganicNamer.Api).

Usage:
    dotnet run --project src/OrganicNamer.Api --urls "http://localhost:5211" &
    python3 tests/http_smoke/run.py [--url http://localhost:5211]
"""

import argparse
import json
import sys
import urllib.error
import urllib.request

from cases import CASES


def post(url, body):
    data = json.dumps(body).encode()
    req = urllib.request.Request(
        url, data=data, headers={"Content-Type": "application/json"}, method="POST"
    )
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            return resp.status, resp.read().decode()
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode()
    except Exception as e:
        return None, str(e)


def check(expect, status, body):
    if expect["type"] == "exact":
        if status != 200:
            return False, f"expected HTTP 200, got {status}: {body}"
        try:
            names = json.loads(body).get("names", [])
        except json.JSONDecodeError:
            return False, f"non-JSON response: {body}"
        if expect["name"] not in names:
            return False, f"expected name '{expect['name']}', got {names}"
        return True, ""

    if expect["type"] == "reject":
        if status is None:
            return False, f"request failed outright: {body}"
        if not (400 <= status < 500):
            return False, f"expected a 4xx rejection, got {status}: {body}"
        return True, ""

    if expect["type"] == "no_crash":
        if status is None:
            return False, f"request failed outright (server likely crashed): {body}"
        try:
            json.loads(body)
        except json.JSONDecodeError:
            return False, f"non-JSON response (possible unhandled exception): {body}"
        return True, ""

    return False, f"unknown expect type {expect['type']}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", default="http://localhost:5211/api/name")
    parser.add_argument("--case", default=None, help="run only the named case")
    args = parser.parse_args()

    names = [args.case] if args.case else list(CASES.keys())
    failures = []

    for name in names:
        if name not in CASES:
            print(f"unknown case: {name}")
            sys.exit(1)
        case = CASES[name]
        status, body = post(args.url, case["payload"])
        ok, reason = check(case["expect"], status, body)
        mark = "PASS" if ok else "FAIL"
        print(f"[{mark}] {name}")
        if not ok:
            print(f"       {reason}")
            failures.append(name)

    print()
    print(f"{len(names) - len(failures)}/{len(names)} passed")
    if failures:
        print("Failed:", ", ".join(failures))
        sys.exit(1)


if __name__ == "__main__":
    main()
