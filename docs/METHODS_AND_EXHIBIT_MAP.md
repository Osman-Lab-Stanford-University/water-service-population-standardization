# Manuscript and repository exhibit map

The table identifies the released data associated with each manuscript exhibit and the type of verification available in this repository.

| Manuscript item | Released result or source | Verification available |
|---|---|---|
| Table 1: analytic counts and standardized means | `main_table_1_observed_standardized_means.csv` | Source totals and mean differences |
| Figure 4(a–b): observed and standardized means | Main Table 1 CSV | Aggregate arithmetic |
| Figure 4(c): selected mean uncertainty | Tables S5/S6 CSVs | Reported working interval endpoints |
| Figure 5(a–b): standardized CDFs | `full_41_threshold_standardized_cdf_percent.csv` | All 41 thresholds and saved E23 output agreement at reported precision |
| Figure 5(c) and Table 2: selected lower-tail probabilities | Main Table 2 CSV | WSQS thresholds 50, 60 and 70 |
| Figure 1(b): drinking–domestic source agreement | SI Table S3 CSV | Source counts and agreement summary |
| Figure 3(c): propensity/weighting diagnostics | Overlap diagnostics CSV | Source-level diagnostics and archived E04 comparisons |
| Figure S1: category balance | Figure 3 source-level SMD maxima CSV | Five source maxima and separately documented 165-combination source verification |
| Figure S2: propensity support | Overlap diagnostics CSV | Reported group-level support summaries |
| Figure S3: pairwise CDF distances | Full standardized CDF CSV | Ten aggregate-derived integrated distances |
| SI Table S1: covariate descriptions | Manuscript Supporting Information | Covariate definitions |
| SI Table S2: questionnaire items | Questionnaire mapping CSV; separate survey repository | Question text and response categories |
| SI Tables S3–S4: domestic-source sensitivity | SI Tables S3/S4 CSVs | Classification counts and standardized mean comparisons |
| SI Tables S5–S6: working intervals | SI Tables S5/S6 CSVs | Reported intervals and contrast arithmetic |
| SI Table S7: alternative domains | SI Table S7 CSV | Released group-level results |
| Eight-variable adjustment sensitivity | `eight_covariate_sensitivity.csv` | Standardized mean comparisons |
| Synthetic simulation coverage | `synthetic_simulation_coverage_summary.csv` | Piped–tanker 39-interior-threshold candidate-band coverage; see simulation note |
| Weight capping at 95th/99th percentiles | E04 analysis code and manuscript | Sensitivity calculations |
| Observed source-specific WSQS dispersion | `observed_wsqs_sd_by_source.csv` | Five observed sample SD values |

After weighting, the largest absolute category-level standardized mean differences remain 0.156 for surface water and 0.188 for tanker delivery. The overlap support screens address a different diagnostic and do not imply uniformly small postweighting imbalance. Working intervals in Tables S5–S6 use independent-record calculations rather than sampling-design or cluster-robust variance estimates.

The [analysis-code guide](ANALYSIS_CODE.md) describes the E04–E23 workflow, and the [data dictionary](DATA_DICTIONARY.md) defines the aggregate tables. The [synthetic simulation summary](SIMULATION_RESULTS.md) describes the E20 candidate-band results.
