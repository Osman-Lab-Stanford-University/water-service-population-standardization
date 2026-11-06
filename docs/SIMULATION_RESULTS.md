# E20 synthetic simultaneous-band experiment

E20 studies two candidate methods for simultaneous confidence bands using a synthetic data-generating process with five source groups, correlated binary characteristics and discretized WSQS outcomes. The nuisance probabilities are known from the generating process; the simulation does not use the household survey.

The four archived aggregate values summarize the **municipal piped–tanker contrast** over the **39 non-endpoint WSQS thresholds** in the `FULL_INTERIOR_39` family. The original E20 program also evaluates a piped–groundwater contrast and two additional threshold families. The full 41-point grid includes the two endpoints, which are excluded from this interior family. The compact four-row table retains the term `full_grid` in its original column name; here it refers to the full interior family.

| Scenario | Replications | Studentized max-statistic coverage | Absolute maximum-error coverage | Nominal |
|---|---:|---:|---:|---:|
| Synthetic null | 240 | 57.9% | 95.0% | 95% |
| Synthetic alternative | 120 | 42.5% | 95.8% | 95% |

The studentized procedure shows substantial undercoverage in this synthetic setting. The unstudentized absolute-error procedure is near nominal coverage in both scenarios. These findings provide a numerical assessment of the candidate procedures under the specified oracle-nuisance data-generating process; they are not an empirical coverage guarantee for full-CDF bands estimated from survey data. Neither procedure is selected for the manuscript's full-grid empirical inference on the basis of this experiment.

The aggregate values are in [`synthetic_simulation_coverage_summary.csv`](../data/aggregates/synthetic_simulation_coverage_summary.csv). The original [E20 program](../code/original/E20/simulate_holdout.py) sets a sample size of 3,066 per replication, seed 20261030 and 1,000 multiplier draws per replication. Additional contrast and threshold-family definitions are specified in that program.
