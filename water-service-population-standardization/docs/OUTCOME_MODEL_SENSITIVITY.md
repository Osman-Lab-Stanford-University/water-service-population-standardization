# Outcome-learner sensitivity (Empirical 27 and Empirical 28)

## Scope

This extension accompanies Section 4.6, Section 5.5 and Supporting Information Section S5/Table S8 of the revised manuscript *Comparing Household Water Service Quality Across Primary Drinking-Water Sources Using Population Standardization*. It checks whether standardized **mean** WSQS depends on the conditional-outcome learner. The original E04 histogram gradient boosting (HGBR) primary estimator and the distributional CDF analysis remain unchanged.

All outcome models use the same 3,066 observations, five drinking-water groups, six background variables (33 category indicators), source-stratified five-fold splits (seed 20261022), original source-membership logistic model and uncapped AIPW estimating equation. Results describe a convenience-recruited analytical cohort, not a representative population or causal effect.

## Learners and tuning

- Primary HGBR: squared-error loss; 200 iterations, learning rate 0.05, up to 7 leaves, depth 3, minimum 10 observations per leaf, L2 = 1.0, no early stopping.
- Random Forest: 500 bootstrap trees; squared-error criterion; all predictors eligible at each split; unrestricted tree depth; minimum 10 observations per leaf; seed 20261022; one thread.
- Ridge: alpha = 10, SVD solver, intercept; predictions constrained to the 0–100 score scale.
- Nested cross-validation HGBR: four fixed configurations were compared within each source-specific outer training fold, using three inner folds and held-out prediction RMSE. The best configuration was fitted to that full outer training subset and predicted its held-out fold. Hyperparameters are in [`outcome_model_sensitivity_specification.json`](../data/metadata/outcome_model_sensitivity_specification.json). The conservative configuration (learning rate 0.05, 150 iterations, 5 leaves, depth 3, minimum leaf 15, L2 = 5) was chosen in all 25 source-by-fold selections. This grid is a sensitivity check, not evidence of global hyperparameter optimality.

The original HGBR was selected before these model-comparison exercises; the alternatives are sensitivity specifications rather than a retrospective claim of predictive superiority.

## Results and interpretation

The corresponding reported values are in [`si_table_s8_outcome_model_sensitivity.csv`](../data/aggregates/si_table_s8_outcome_model_sensitivity.csv). The group ordering remains surface water, groundwater, piped, commercial, tanker under all learners. Across five sources, the largest absolute change relative to primary HGBR is 0.406 WSQS points (tanker under Ridge). The piped-minus-groundwater contrast is −0.602 for original HGBR, −0.656 for Random Forest, −0.674 for Ridge and −0.666 for CV-selected HGBR.

Pooled out-of-fold RMSE is 16.095, 15.850, 15.763 and 15.913, respectively. The three alternatives improve this predictive diagnostic modestly, but lower predictive error need not imply a meaningfully different or more precise AIPW estimate. The paired RMSE comparisons use 2,000 source-stratified bootstrap resamples of the held-out prediction errors, without refitting learners. Their percentile intervals are conditional on the fitted models and are not design-based confidence intervals.

Working AIPW standard errors are computed from the variation of centered empirical score contributions, with paired per-observation differences used for contrasts. They treat observations as independent and do not account fully for nuisance fitting or potential dependence from recruitment. They are labeled **working** diagnostics, not publication-certified inference. The primary mean/CDF uncertainty results remain those reported separately in the manuscript.

## Verification and empirical replay

The saved E28 rerun verified the original E02 logistic out-of-fold log loss and 30 source-group support summaries, reproduced 20 original E04 HGBR aggregates, and reconciled 60 values with the earlier E27 stage. The aggregate figures here were transcribed from the **successful author-supplied Empirical 28 report**, preserving reported numerical precision; no participant-level files are in this repository. Both E27 and E28 source files are retained with SHA-256 hashes under [`code/original`](../code/README.md).

From the repository root, anyone can run the standard-library aggregate checks without access to survey records:

```bash
python3 -B scripts/run_public_checks.py
```

With the dependencies in `requirements-reference.txt`, the additional synthetic-only code test is:

```bash
python3 -B tests/test_outcome_learner_synthetic.py
```

Full empirical refitting requires authorized, restricted original source data, the original decoder and scoring crosswalk, the original E04 stage reports, and the previously completed E27 aggregate stage. The archived programs require those files as explicit arguments. They must be run in a separate private work directory, not with respondent-level outputs committed to GitHub. This package is therefore a verifiable aggregate-plus-source release, not an independently runnable original-data reproduction.

## Files

| Path under `data/aggregates` | Rows | Content |
|---|---:|---|
| `si_table_s8_outcome_model_sensitivity.csv` | 8 | Five mean WSQS values, selected mean contrast, pooled prediction RMSE and working contrast SE across four learners |
| `outcome_model_oof_rmse.csv` | 6 | Source-specific and pooled out-of-fold RMSE |
| `outcome_model_paired_rmse_bootstrap.csv` | 18 | Three learner-versus-HGBR paired RMSE differences and conditional bootstrap percentile endpoints across five sources and pooled |
| `outcome_model_hgbr_tuning_counts.csv` | 5 | Number of outer folds selecting each of the four HGBR configurations |
| `outcome_model_working_if_means.csv` | 20 | Source-specific mean and working influence-score SE for each of four learners |
| `outcome_model_working_if_piped_groundwater.csv` | 4 | Contrast estimates, working SEs, learner shifts and working SE of paired shifts |

## Model-method references

Breiman, L. (2001). Random forests. *Machine Learning*, 45, 5–32. https://doi.org/10.1023/A:1010933404324

Hoerl, A. E., & Kennard, R. W. (1970). Ridge regression: Biased estimation for nonorthogonal problems. *Technometrics*, 12(1), 55–67. https://doi.org/10.1080/00401706.1970.10488634

Varma, S., & Simon, R. (2006). Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics*, 7, 91. https://doi.org/10.1186/1471-2105-7-91

Friedman, J. H. (2001). Greedy function approximation: A gradient boosting machine. *The Annals of Statistics*, 29(5), 1189–1232. https://doi.org/10.1214/aos/1013203451
