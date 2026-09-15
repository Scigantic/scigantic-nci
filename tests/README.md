# Tests

The suite reads a local copy of the two public mirror buckets and never needs credentials.
`tests/conftest.py` picks the copy in this order: `$SCIGANTIC_NCI_TEST_ROOT`, the full
build-time mirrors when they exist on the machine that built the buckets, then
`tests/fixtures/mirror`.

`tests/fixtures/manifest.txt` lists every object the suite reads (bucket, key, bytes;
about 270 MB across 718 objects). Fetch them once and run the tests:

```
python tests/fixtures/fetch_mirror.py
pytest -q
```

CI does exactly that, with `tests/fixtures/mirror` cached on the manifest hash. Set
`SCIGANTIC_NCI_TEST_FIXTURES_ONLY=1` to force the fixture copy locally and reproduce a CI run.

When a test starts reading a new table, regenerate the manifest against a full mirror:

```
python -m pytest -q -p tests.fixtures.record_manifest tests
```

The three methylation tests need a TCGA-CHOL build that includes `methylation/`, which is
not in the public bucket; they skip everywhere else. Tests marked `network` read the live
buckets and the raw NCI buckets directly (`-m "not network"` deselects them).
