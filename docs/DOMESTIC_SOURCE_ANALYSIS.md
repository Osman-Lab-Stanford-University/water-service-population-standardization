# Primary domestic-water source sensitivity (Q1.1)

The primary analysis classifies the 3,066 household reports by **primary drinking-water source (Q1.3)**. The supplementary sensitivity analysis instead uses **primary domestic-water source (Q1.1)**, which includes water used for household purposes beyond drinking.

Within the drinking-water analytical cohort, 2,012 of 3,066 reports (65.62%) have matching domestic and drinking source classifications. The five modeled domestic-water groups contain 3,044 reports (municipal piped 1,335; groundwater 1,449; surface 160; tanker 50; commercial 50); the remaining 22 have a domestic-source classification of *Other*.

| Source category | Q1.3 standardized WSQS | Q1.1 standardized WSQS | Q1.1 ESS / observed group |
|---|---:|---:|---:|
| Municipal piped | 74.847 | 74.908 | 96.7% |
| Groundwater | 75.449 | 75.148 | 97.3% |
| Surface water | 76.365 | 73.363 | 89.6% |
| Tanker delivery | 71.530 | 77.510 | 41.2% |
| Commercial water | 72.299 | 67.759 | 65.3% |

The standardized piped-minus-groundwater difference is approximately −0.60 WSQS points using Q1.3 and −0.24 using Q1.1. Domestic-source tanker and commercial groups each contain 50 observations. Their effective sample fractions are 41.2% and 65.3%, respectively; maximum normalized weights are 12.757% and 6.410%, and maximum balance gaps are 17.51 and 8.11 percentage points. The estimated domestic-minus-drinking standardized mean differences are +5.981 WSQS points for tanker delivery and −4.539 for commercial water. The two specifications also differ in group classification and analytical-sample composition, with 22 domestic-source *Other* reports excluded from the Q1.1 analysis. The Q1.1 sensitivity reports standardized means and support diagnostics, not a complete distributional re-estimation.

The aggregate values are provided in [`si_table_s3_drinking_domestic_source_comparison.csv`](../data/aggregates/si_table_s3_drinking_domestic_source_comparison.csv) and [`si_table_s4_domestic_source_sensitivity.csv`](../data/aggregates/si_table_s4_domestic_source_sensitivity.csv). Additional methodological details are available in [Numerical provenance](PROVENANCE.md).
