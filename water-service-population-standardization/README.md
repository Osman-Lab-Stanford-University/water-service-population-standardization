# Household Water Service Quality Standardization

Research code and aggregate results accompanying *Comparing Household Water Service Quality Across Primary Drinking-Water Sources Using Population Standardization*, a manuscript prepared for *Water Resources Research*.

## Study overview

The study compares household-reported water-service quality across five groups defined by primary drinking-water source (Q1.3). The Water Service Quality Score (WSQS, 0–100) summarizes eight dimensions of household water-service experience.

The analysis draws on 3,066 household reports from a larger survey conducted in Tamil Nadu, India, in September 2025. Recruitment used convenience sampling with support from local partners. The empirical application focuses on a Karur-centered analytical subset.

Five-fold cross-fitted augmented inverse-probability weighting (AIPW) standardizes the source-specific WSQS means and distributions to the empirical joint distribution of six measured respondent and household characteristics. The results describe the analytical sample and do not establish causal effects or regional population estimates.

A supplementary sensitivity analysis uses **primary domestic-water source (Q1.1)** in place of primary drinking-water source (Q1.3). The classifications agree for 2,012 of 3,066 households (65.6%). The five domestic-source groups contain 3,044 households, excluding 22 classified as *Other*. See [Domestic-source sensitivity](docs/DOMESTIC_SOURCE_ANALYSIS.md).

The E27–E28 sensitivity analysis compares Random Forest, Ridge and nested-CV-selected histogram gradient boosting with the original conditional-outcome learner, retaining the original target and AIPW score. See [outcome-model sensitivity](docs/OUTCOME_MODEL_SENSITIVITY.md) and aggregate Supporting Information Table S8.

The separate [Household Water Service Survey (2025)](https://github.com/Osman-Lab-Stanford-University/water-service-survey-2025) repository provides English and Tamil reference questionnaires and supporting documentation.

## Contents

- [`data/aggregates/`](data/README.md) — 20 tables containing observed and standardized WSQS results, cumulative distributions, overlap diagnostics, working uncertainty intervals and sensitivity analyses.
- [`code/original/`](code/README.md) — research programs for AIPW mean and CDF estimation, uncertainty, simulation, figure generation, and the E27–E28 outcome-model sensitivity analyses.
- [`scripts/`](docs/REPRODUCING_RESULTS.md) and [`tests/`](docs/REPRODUCING_RESULTS.md) — aggregate-result and source-integrity checks. [`src/aipw_reference.py`](src/aipw_reference.py) is a separate synthetic-data reference implementation.
- [`examples/`](examples/README.md) — figures generated from the included aggregate tables, distinct from the manuscript artwork.

## Verification

From the repository root, with Python 3.10 or later:

```bash
python3 -B scripts/run_public_checks.py
```

The checks use the Python standard library. See [Reproducing Results](docs/REPRODUCING_RESULTS.md) for figure-generation and synthetic-test instructions and the [Data Dictionary](docs/DATA_DICTIONARY.md) for variable definitions.

## Data and code availability

The repository includes aggregate results and the associated analysis programs. [Numerical provenance](docs/PROVENANCE.md) describes comparisons with saved computational results; the [analysis-code guide](docs/ANALYSIS_CODE.md) identifies the programs and dependencies.

The participant-level survey data are restricted. Complete empirical model re-estimation requires the protected records, original scorer and decoder, and intermediate stage inputs. The released tables and synthetic tests can be checked without those materials.

## Citation and licensing

The manuscript citation and a persistent archive identifier will be added following publication and archival registration. Repository citation metadata are supplied in [`CITATION.cff`](CITATION.cff).

The code is offered under the [MIT License](LICENSE-CODE), and aggregate data and documentation under [CC BY 4.0](LICENSE-DATA). See [license coverage](LICENSES.md) for file-level scope.
