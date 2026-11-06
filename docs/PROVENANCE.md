# Numerical provenance

The repository includes 14 aggregate tables corresponding to the manuscript and Supporting Information. Their published-precision values were compared with saved computational results.

## Primary analysis

The analytical cohort contains 3,066 households in five groups defined by primary drinking-water source (Q1.3). Five-fold cross-fitted AIPW uses a common reference distribution of six respondent and household characteristics. The Karur-centered sample was recruited by convenience sampling with support from local partners.

## Results and supporting evidence

| Result | Evidence |
|---|---|
| Standardized WSQS means | E04 group means, sample sizes and overlap diagnostics agree with the source calculations. |
| Standardized CDFs | The E23 register matches 205 source-specific cumulative proportions and 41 piped–groundwater differences at reported precision; the separate 25H comparison agrees. |
| Domestic-source sensitivity | The Q1.1 five-group results agree with saved 26E/26F analyses. |
| Figure S1 | Numerical comparisons cover 165 indicator-by-source combinations, the associated figure fingerprint and the display-screening flags. |
| Analysis programs | SHA-256 values for the nine source programs appear in [`code/ORIGINAL_SOURCE_SHA256SUMS.txt`](../code/ORIGINAL_SOURCE_SHA256SUMS.txt). |

Detailed values and tolerances are in [CDF verification](CDF_VERIFICATION.md) and [Figure S1 verification](FIGURE_S1_VERIFICATION.md). The [code guide](ANALYSIS_CODE.md) describes the analytical workflow.

## Reproduction scope

The included tests reproduce aggregate arithmetic and selected figures. Empirical model re-estimation requires restricted respondent records, survey-scoring and decoding programs, and the corresponding stage inputs.
