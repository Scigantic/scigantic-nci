"""Fetch the test fixture mirror from the public buckets.

Reads tests/fixtures/manifest.txt (bucket, key, bytes) and downloads each object over
plain HTTPS into tests/fixtures/mirror/<gdc_mirror|idc_mirror>/<key>, the same layout as
the full local mirrors. Files already present with the recorded size are skipped, so the
directory can be cached between CI runs (keyed on the manifest hash).

An object whose live size differs from the manifest (the bucket owners rebuilt it) is
fetched anyway and listed at the end; the tests then run against the live object. Pass
--update-manifest to write the live sizes back into manifest.txt. A truncated download
is still an error.

    python tests/fixtures/fetch_mirror.py [--dest DIR] [--workers N] [--update-manifest]
"""
from __future__ import annotations

import argparse
import os
import sys
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.join(HERE, "manifest.txt")
DEFAULT_DEST = os.path.join(HERE, "mirror")
LAYOUT = {"scigantic-gdc-open": "gdc_mirror", "scigantic-idc-open": "idc_mirror"}
CHANGED: list[tuple[str, str, int, int]] = []
GONE: list[tuple[str, str]] = []


def entries() -> list[tuple[str, str, int]]:
    out = []
    with open(MANIFEST) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            bucket, key, size = line.rsplit(" ", 2)
            out.append((bucket, key, int(size)))
    return out


def fetch(bucket: str, key: str, size: int, dest: str) -> tuple[str, bool]:
    path = os.path.join(dest, LAYOUT[bucket], key)
    if os.path.isfile(path) and os.path.getsize(path) == size:
        return path, False
    os.makedirs(os.path.dirname(path), exist_ok=True)
    url = f"https://{bucket}.s3.amazonaws.com/{urllib.request.quote(key)}"
    tmp = path + ".part"
    try:
        resp = urllib.request.urlopen(url, timeout=120)
    except urllib.error.HTTPError as e:
        if e.code in (403, 404):  # removed from the bucket since the manifest was recorded
            GONE.append((bucket, key))
            return path, False
        raise
    with resp, open(tmp, "wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    got = os.path.getsize(tmp)
    declared = resp.headers.get("Content-Length")
    if declared is not None and got != int(declared):
        os.remove(tmp)
        raise RuntimeError(f"{url}: truncated download, {got} of {declared} bytes")
    os.replace(tmp, path)
    if got != size:
        # the bucket object was legitimately rebuilt since the manifest was recorded; the tests
        # run against what is live, and a stale expectation shows up as a test failure instead
        CHANGED.append((bucket, key, size, got))
    return path, True


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dest", default=DEFAULT_DEST)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--update-manifest", action="store_true", help="write live sizes of changed objects into manifest.txt")
    a = ap.parse_args(argv)
    items = entries()
    total = sum(s for _, _, s in items)
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        results = list(ex.map(lambda e: fetch(*e, a.dest), items))
    fetched = sum(1 for _, new in results if new)
    print(f"{len(items)} objects, {total / 1e6:.1f} MB in manifest; fetched {fetched}, already present {len(items) - fetched - len(GONE)}; mirror at {a.dest}")
    if GONE:
        print(f"{len(GONE)} objects are no longer in the buckets; re-record the manifest (see tests/README.md):")
        for bucket, key in sorted(GONE):
            print(f"  {bucket}/{key}")
        return 1
    if CHANGED:
        print(f"{len(CHANGED)} objects changed size since the manifest was recorded (tests use the live copies):")
        for bucket, key, want, got in sorted(CHANGED):
            print(f"  {bucket}/{key}: {want} -> {got} bytes")
        if a.update_manifest:
            live = {(b, k): g for b, k, _, g in CHANGED}
            with open(MANIFEST) as fh:
                lines = fh.read().splitlines()
            out = []
            for line in lines:
                if line.strip() and not line.startswith("#"):
                    bucket, key, _ = line.rsplit(" ", 2)
                    if (bucket, key) in live:
                        line = f"{bucket} {key} {live[(bucket, key)]}"
                out.append(line)
            with open(MANIFEST, "w") as fh:
                fh.write("\n".join(out) + "\n")
            print(f"updated {len(live)} sizes in {MANIFEST}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
