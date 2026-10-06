# Contributing

Run the documented environment and tests before submitting a change:

```bash
python -m pip install -r requirements.txt
MPLBACKEND=Agg python -m pytest -q
```

Describe the observed failure, expected behavior, dependency versions and the
smallest input that reproduces it. Include a regression test for corrected
behavior; prefer known synthetic signals, independent numerical expectations
and real public interfaces over tests that duplicate implementation details.

Explain changed units, boundary conditions, random seeds or scientific assumptions.
Do not change numerical expectations only to silence a failing test. Keep examples
small and runnable without private files, preserve input data, and avoid downloads
or processing during import. Never commit credentials or restricted datasets.

Report skipped integration checks explicitly. A passing synthetic benchmark does
not prove scientific validity for every dataset; document untested conditions.
