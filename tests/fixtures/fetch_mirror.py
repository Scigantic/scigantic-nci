"""Fetch the test fixture mirror from the public buckets.

Reads tests/fixtures/manifest.txt (bucket, key, bytes) and downloads each object over
plain HTTPS into tests/fixtures/mirror/<gdc_mirror|idc_mirror>/<key>, the same layout as
the full local mirrors. Files already present with the right size are skipped, so the
directory can be cached between CI runs (keyed on the manifest hash).

    python tests/fixtures/fetch_mirror.py [--dest DIR] [--workers N]
"""
from __future__ import annotations

import argparse
import os
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.join(HERE, "manifest.txt")
DEFAULT_DEST = os.path.join(HERE, "mirror")
LAYOUT = {"scigantic-gdc-open": "gdc_mirror", "scigantic-idc-open": "idc_mirror"}


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
    with urllib.request.urlopen(url, timeout=120) as resp, open(tmp, "wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    got = os.path.getsize(tmp)
    if got != size:
        os.remove(tmp)
        raise RuntimeError(f"{url}: expected {size} bytes, got {got}")
    os.replace(tmp, path)
    return path, True


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dest", default=DEFAULT_DEST)
    ap.add_argument("--workers", type=int, default=16)
    a = ap.parse_args(argv)
    items = entries()
    total = sum(s for _, _, s in items)
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        results = list(ex.map(lambda e: fetch(*e, a.dest), items))
    fetched = sum(1 for _, new in results if new)
    print(f"{len(items)} objects, {total / 1e6:.1f} MB in manifest; fetched {fetched}, already present {len(items) - fetched}; mirror at {a.dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
