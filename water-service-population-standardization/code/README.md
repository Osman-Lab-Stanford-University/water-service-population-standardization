# Analysis code

The `original/` directory contains the nine original E04–E23 analysis programs and the four E27–E28 sensitivity source files identified in the [code guide](../docs/ANALYSIS_CODE.md). Their SHA-256 values are listed in [`ORIGINAL_SOURCE_SHA256SUMS.txt`](ORIGINAL_SOURCE_SHA256SUMS.txt).

| Component | Purpose |
|---|---|
| `original/E04/` and `original/E04_engine/` | Cross-fitted AIPW mean standardization |
| `original/E09/` | Four-contrast working uncertainty |
| `original/E15/` | Standardized CDF estimation |
| `original/E20/` | Synthetic validation |
| `original/E23/` | Numerical reconciliation |
| `original/S1_generator_1/` | Supplementary-figure generation |
| `original/dependency_*` | Additional support modules |
| `original/E27/` | E27 outcome learner comparison (HGBR, Random Forest, Ridge) |
| `original/E28/` | Exact E27 copy, in-memory OOF interface, E28 nested HGBR/paired RMSE/working IF diagnostics |

`runner_support/` provides a Python import layout for the three common modules. The empirical programs require restricted survey inputs, scoring and decoding components, and intermediate results. E04 includes outcome-model fitting, and E09 extends the analysis to lower-tail estimation and four-contrast working inference.

From the repository root:

```bash
python3 -B scripts/verify_original_code_archive.py
```

The E27 and E28 scripts require authorized restricted inputs, and the E28 run additionally verifies the original successful E27 aggregate stage. The repository does not contain those inputs. See [outcome-model sensitivity](../docs/OUTCOME_MODEL_SENSITIVITY.md).
