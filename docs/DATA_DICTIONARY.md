# Data dictionary — aggregate CSVs

The 14 CSV files in `data/aggregates/` contain the manuscript and Supporting Information aggregates. The tables below identify their fields and units.

| File | Rows | Description |
|---|---:|---|
| `eight_covariate_sensitivity.csv` | 5 | Six-X and eight-X standardized means with difference in score points. |
| `figure3_source_level_max_absolute_smd.csv` | 5 | Five source-group maxima of absolute binary-indicator SMD before/after weighting; fixed common-reference Bernoulli denominator. |
| `full_41_threshold_standardized_cdf_percent.csv` | 41 | WSQS threshold (points), group CDF probabilities (percent), and piped-minus-groundwater gap (percentage points); 41 rows. |
| `main_table_1_observed_standardized_means.csv` | 5 | Source sample count and observed/standardized WSQS means (0–100 points), differences in points. |
| `main_table_2_standardized_threshold_probabilities.csv` | 5 | Selected standardized CDF probabilities for WSQS at or below 50, 60 and 70 (percent). |
| `observed_wsqs_sd_by_source.csv` | 5 | Unweighted observed household WSQS sample standard deviation in score points. |
| `overlap_diagnostics_original_aggregate.csv` | 5 | Source-specific sample size, effective sample size, propensity tail, weight share and maximum balance gaps; see columns. |
| `si_table_s2_questionnaire_item_mapping.csv` | 8 | Public English survey question text and response categories for source/adjustment variables. |
| `si_table_s3_drinking_domestic_source_comparison.csv` | 5 | Source-group sizes and observed WSQS under two reported primary-source classifications. |
| `si_table_s4_domestic_source_sensitivity.csv` | 5 | Domestic-source sensitivity standardized means and source-specific support diagnostics. |
| `si_table_s5_four_contrast_working_intervals.csv` | 4 | Four selected piped–groundwater point estimates and working marginal/four-contrast intervals. |
| `si_table_s6_ten_pairwise_mean_intervals.csv` | 10 | Ten pairwise standardized mean contrasts and conservative simultaneous working intervals. |
| `si_table_s7_core_domain_sensitivity.csv` | 5 | Standardized group means for CORE-8 and two alternative seven-domain composites. |
| `synthetic_simulation_coverage_summary.csv` | 4 | Piped–tanker coverage for two simultaneous-band candidates in the synthetic E20 experiment; [design and interpretation](SIMULATION_RESULTS.md). |

## Field names

The actual headers below are the machine-readable variable definitions. Percentages are on 0–100 scale unless a column explicitly says proportion or fraction (0–1). WSQS points are on its 0–100 scoring scale. Differences between CDFs are in **percentage points**, not percentages.


### `eight_covariate_sensitivity.csv`

`primary_drinking_source`, `six_X_standardized_WSQS`, `eight_X_standardized_WSQS`, `eight_X_minus_six_X_WSQS`

### `figure3_source_level_max_absolute_smd.csv`

`source`, `n`, `max_abs_SMD_before`, `max_abs_SMD_after`

### `full_41_threshold_standardized_cdf_percent.csv`

`WSQS_threshold`, `CDF_PIPED_percent`, `CDF_GROUNDWATER_percent`, `CDF_SURFACE_percent`, `CDF_TANKER_percent`, `CDF_COMMERCIAL_percent`, `PIPED_minus_GROUNDWATER_pp`

### `main_table_1_observed_standardized_means.csv`

`Primary drinking-water source`, `n`, `Observed mean WSQS`, `Standardized mean WSQS`, `Observed − standardized`

### `main_table_2_standardized_threshold_probabilities.csv`

`Primary drinking-water source`, `WSQS ≤50 (%)`, `WSQS ≤60 (%)`, `WSQS ≤70 (%)`

### `observed_wsqs_sd_by_source.csv`

`source`, `n`, `observed_WSQS_sample_SD`

### `overlap_diagnostics_original_aggregate.csv`

`source`, `n`, `ess`, `ess_fraction`, `largest_normalized_weight_share`, `max_weighted_category_prop_gap_pp`, `p01_target`, `E04_screen_pass`

### `si_table_s2_questionnaire_item_mapping.csv`

`QID`, `Measure`, `English questionnaire wording`, `Response options`

### `si_table_s3_drinking_domestic_source_comparison.csv`

`Source`, `Drinking n`, `Domestic n`, `Same source n (%)`, `Drinking observed WSQS`, `Domestic observed WSQS`

### `si_table_s4_domestic_source_sensitivity.csv`

`Source`, `Drinking standardized WSQS`, `Domestic standardized WSQS`, `Domestic − Drinking`, `Domestic ESS/n (%)`, `Maximum weight share (%)`, `Maximum balance gap (pp)`

### `si_table_s5_four_contrast_working_intervals.csv`

`Contrast`, `Estimate`, `Marginal 95% interval`, `Four-contrast simultaneous 95% interval`

Both interval fields are *working* 95% intervals computed from independent-record influence-score approximations, without probability-sampling or recruitment-cluster design variance.

### `si_table_s6_ten_pairwise_mean_intervals.csv`

`Contrast (first source minus second)`, `Difference (WSQS points)`, `Simultaneous 95% interval`

The interval is a conservative working simultaneous interval derived from the five-mean Bonferroni rectangle under an independent-record approximation; it is not a clustered or sampling-design interval.

### `si_table_s7_core_domain_sensitivity.csv`

`Source`, `CORE-8`, `CORE7A`, `CORE7B`

### `synthetic_simulation_coverage_summary.csv`

`simulation_source`, `scenario`, `method`, `full_grid_coverage_percent`, `nominal_percent`

These rows describe the piped–tanker contrast over E20's 39 interior thresholds (240 null and 120 alternative replications). The source column name `full_grid_coverage_percent` refers to this interior family, not all 41 points. See [simulation results](SIMULATION_RESULTS.md) and the [machine-readable comparison scope](../data/metadata/e20_coverage_scope.json).

## CORE-8 WSQS score

Eight questionnaire dimensions contribute equally when applicable: affordability, access, reliability, efficiency, acceptability, perceived water quality, quantity and customer service. The empirical analysis includes households with **seven or eight scored responses**. For substantive item-severity codes `s = 1, 2, 3, 4, 5` (1 best; 5 worst), the corresponding dimension score is `25 × (5 − s)`, ranging from 100 to 0. WSQS is the arithmetic mean of the seven or eight applicable dimension scores. The original response-label-to-severity coding crosswalk is maintained separately and is not included in this archive.

## Important conventions

- The 41 evaluated thresholds are 0, 2.5, ..., 100 WSQS points.
- `PIPED_minus_GROUNDWATER_pp` means the piped cumulative proportion **minus** the groundwater cumulative proportion in percentage points.
- Summary statistics are descriptive of the retained Karur-recruitment analytical records, not estimates for all households in Karur or Tamil Nadu.
- The [English and Tamil questionnaire repository](https://github.com/Osman-Lab-Stanford-University/water-service-survey-2025) maintains the survey reference instruments separately.
- The Figure S1 category-level source table is not distributed.
