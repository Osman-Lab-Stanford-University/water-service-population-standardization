# Analysis programs

The `code/original/` directory contains the nine original E04–E23 programs and the four E27–E28 outcome-model sensitivity source files. File hashes are recorded for reproducibility.

## Source inventory

| Program | Role | SHA-256 |
|---|---|---|
| `E04/run_aipw.py` | Primary six-characteristic AIPW means | `439fc907ff72fdc93071bb61e4ea37cd0ecff68b2948a2fd6942208bcf9d9bf1` |
| `E04_engine/core.py` | Frozen AIPW scoring engine | `10fafe121b64b9ff53bcb992c076d45a7b1097fd0ba76e7b1702dbb2e0cd4645` |
| `E09/logistic_inference.py` | Four-contrast working uncertainty | `f12f0240cf929003dabc6eab188f52ca28aa0a78a7add186852d562cc6014222` |
| `E15/run_cdf.py` | Full 41-threshold standardized CDF estimation | `3277b60e82d4d5e4714ddd7d15e4b6478797284391dd8a6b320db8404346708d` |
| `E20/simulate_holdout.py` | Synthetic oracle-nuisance validation | `42149a134cfba1c3dc63aa594aa74f0731bcf0d5d9de7a40a082920a2737effe` |
| `E23/reconcile.py` | Reconciliation of saved empirical outputs | `df85dbcd7f6d6e8d0371590ddac1c9cfe8bc3cc8f136e776b306f9479056d3d2` |
| `S1_generator_1/build_si_figures.py` | Supporting-figure generation | `6985a2ea05c7782cfb80bd1881f66f2e9bb529d3bc0e87adaa42848a28895e11` |
| `dependency_support_policy/support_policy.py` | Overlap/support diagnostics | `d2ae5f0484e58d216dcc2deead8eb27a2452fa2b09cac562a8f34f7f3182beca` |
| `dependency_demo/demo.py` | Original synthetic demonstration | `eec9c3e94b33085b6784ce8d65622ab9e33284d7f1c7fdc72ffb5abb8bf8d39b` |

Paths in the table are relative to `code/original/`.

## Analytical workflow

**E04** constructs the five-source analytical cohort (N=3,066) and estimates cross-fitted propensity and outcome models for the six-characteristic standardization. The `E04_engine` module supplies AIPW scoring functions.

**E09** calculates four paired piped–groundwater contrasts and working uncertainty intervals from influence scores, including multiplier-based joint intervals. The intervals use independent-record calculations rather than survey-design or cluster-robust variance estimation.

**E15** estimates standardized cumulative distributions at 41 WSQS thresholds (0–100, spaced 2.5 points apart) using cross-fitted source-specific logistic outcome models.

**E20** compares two candidate simultaneous-band procedures in oracle-nuisance simulations. On the piped–tanker contrast across 39 interior thresholds, the studentized method covered the true curve in 57.9% of null and 42.5% of alternative replications, compared with 95.0% and 95.8% for the absolute-error method. These findings concern the synthetic candidate bands, not the fitted empirical WSQS point estimates; see [simulation results](SIMULATION_RESULTS.md). **E23** reconciles numerical results from earlier estimation stages.

The **Figure S1 generator** computes category-level balance diagnostics from 33 binary indicators representing six adjustment variables.

## Code requirements

The complete empirical workflow uses the original survey records, the questionnaire decoder, CORE-8 scoring code and response crosswalk, and intermediate stage inputs. Those restricted and upstream materials are not included here.

`code/runner_support/src/water_standardize/` supplies the core, support-policy and demonstration modules needed by the archived Python imports. The additional `src/aipw_reference.py` is a standalone synthetic-data example; the verification scripts operate on the released aggregate tables.

To verify the code files and Python syntax from the repository root:

```bash
python3 -B scripts/verify_original_code_archive.py
```

## Outcome-model learner sensitivity (E27–E28)

Archived programs: `code/original/E27/run_outcome_model_sensitivities.py`, `code/original/E28/E27_REFERENCE_UNMODIFIED.py`, `code/original/E28/E27_extended_in_memory.py` and `code/original/E28/run_e28_diagnostics.py`. E28 imports the adjacent in-memory E27 extension; the unmodified archived reference lets the E28 run verify the prior E27 stage. See [methodology, public aggregate files and restricted prerequisites](OUTCOME_MODEL_SENSITIVITY.md).
