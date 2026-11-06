# Software environment and model settings

The aggregate verification commands require Python 3.10 or later and use the standard library. Optional reference tests use NumPy, SciPy, scikit-learn and pandas; aggregate-based figures use NumPy and Matplotlib. See `requirements-reference.txt` and `requirements-figures.txt` at the repository root.

E04 uses five source-stratified folds (seed 20261022). Propensity estimation uses L2 logistic regression (C=1.0, `lbfgs`, maximum 2,000 iterations, tolerance 1e-6). Source-specific outcome models use histogram gradient boosting with 200 iterations, learning rate 0.05, a seven-leaf maximum, depth three, minimum leaf count ten, L2 regularization 1.0, no early stopping and seed 20261022. E15 fits logistic threshold models with up to 3,000 iterations and tolerance 1e-6 at 41 thresholds from 0 to 100 in increments of 2.5. E09 uses 4,000 multiplier replications for four-contrast working intervals.

The original computation environment was not saved as a dependency lockfile. Empirical execution requires the restricted survey records, decoder, scorer, coding crosswalk and intermediate results.
