# Analysis code

The `original/` directory contains the nine analysis programs identified in the [code guide](../docs/ANALYSIS_CODE.md). Their SHA-256 values are listed in [`ORIGINAL_SOURCE_SHA256SUMS.txt`](ORIGINAL_SOURCE_SHA256SUMS.txt).

| Component | Purpose |
|---|---|
| `original/E04/` and `original/E04_engine/` | Cross-fitted AIPW mean standardization |
| `original/E09/` | Four-contrast working uncertainty |
| `original/E15/` | Standardized CDF estimation |
| `original/E20/` | Synthetic validation |
| `original/E23/` | Numerical reconciliation |
| `original/S1_generator_1/` | Supplementary-figure generation |
| `original/dependency_*` | Additional support modules |

`runner_support/` provides a Python import layout for the three common modules. The empirical programs require restricted survey inputs, scoring and decoding components, and intermediate results. E04 includes outcome-model fitting, and E09 extends the analysis to lower-tail estimation and four-contrast working inference.

From the repository root:

```bash
python3 -B scripts/verify_original_code_archive.py
```
