# Checking the results

Run the commands from the repository root.

## Aggregate verification

```bash
python3 -B scripts/run_public_checks.py
```

The checks cover all 14 aggregate tables, the 41-threshold CDF grid, the 50/60/70 threshold results, pairwise CDF distances, mean-score arithmetic, working interval endpoints, domestic-source sensitivity, overlap summaries, source-program fingerprints, file allowlisting, basic direct-identifier screening and file checksums. Tests of malformed inputs run in temporary copies.

## Reconstruct figures from aggregate tables

```bash
python3 -m pip install -r requirements-figures.txt
python3 -B scripts/regenerate_aggregate_figures.py --out generated_figures
python3 -B scripts/reconstruct_figure_s3_distances.py --out generated_figures
```

These figures illustrate the released aggregate results and are distinct from the original manuscript artwork. Pairwise CDF distances are integrated numerically from the rounded 41-threshold grid.

## Synthetic-data reference tests

```bash
python3 -m pip install -r requirements-reference.txt
python3 -B tests/test_reference_method.py
python3 -B tests/test_reference_method_nontrivial.py
```

The tests exercise the separate AIPW reference implementation with synthetic data.

The E20 coverage summary is also checked against its machine-readable [scope metadata](../data/metadata/e20_coverage_scope.json) in the main verification command.

## E20 synthetic simulation test

```bash
python3 -B tests/test_original_e20_smoke.py
```

The E20 test runs a reduced synthetic case. The simulation summary is in [`data/aggregates/synthetic_simulation_coverage_summary.csv`](../data/aggregates/synthetic_simulation_coverage_summary.csv).

## Source-code verification

```bash
python3 -B scripts/verify_original_code_archive.py
```

This command verifies the source files' SHA-256 values and syntax. See [Analysis code](../code/README.md) and [Software environment](SOFTWARE_ENVIRONMENT.md) for execution requirements.
