"""Statistical core. Does not read survey files."""

import numpy as np
import pandas as pd

from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

REGIMES = (
    "PIPED", "GROUNDWATER", "SURFACE",
    "TANKER", "COMMERCIAL"
)

NUMERIC = (
    "income", "vehicles", "household_size",
    "children_u18", "children_u5"
)

CATEGORICAL = (
    "dwelling", "tenure", "upbringing"
)

FEATURES = NUMERIC + CATEGORICAL

INDICATORS = tuple(
    name + "_missing" for name in NUMERIC
)

LEVELS = {
    "dwelling": (
        "HOUSE", "APARTMENT", "TEMPORARY",
        "OTHER", "MISSING"
    ),
    "tenure": (
        "RENT", "OWN_CLEAR", "OWN_LOAN",
        "OTHER", "MISSING"
    ),
    "upbringing": (
        "URBAN", "SUBURBAN", "RURAL", "MISSING"
    )
}

SEED = 20261022


def validate(frame):
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("Expected pandas DataFrame")

    needed = set(FEATURES) | {"regime", "wsqs"}
    absent = needed - set(frame.columns)

    if absent:
        raise ValueError(
            "Missing required columns: " +
            ",".join(sorted(absent))
        )

    if frame.empty or frame.index.has_duplicates:
        raise ValueError("Empty data or duplicate indices")

    if (
        frame["regime"].isna().any()
        or not frame["regime"].isin(REGIMES).all()
        or set(frame["regime"]) != set(REGIMES)
    ):
        raise ValueError("Five populated regimes required")

    y = pd.to_numeric(
        frame["wsqs"], errors="coerce"
    ).to_numpy(dtype=float)

    if not (
        np.isfinite(y).all()
        and ((y >= 0) & (y <= 100)).all()
    ):
        raise ValueError("WSQS outside valid range")

    for field in NUMERIC:
        v = pd.to_numeric(
            frame[field], errors="coerce"
        )

        bad = frame[field].notna() & v.isna()

        if (
            bad.any()
            or np.isinf(v.to_numpy(float)).any()
            or v.dropna().empty
        ):
            raise ValueError(
                f"Invalid numeric variable: {field}"
            )

    for field in CATEGORICAL:
        values = frame[field].dropna()

        if not values.isin(LEVELS[field]).all():
            raise ValueError(
                f"Invalid categorical variable: {field}"
            )

    both = frame[
        ["children_u18", "children_u5"]
    ].notna().all(axis=1)

    if (
        frame.loc[both, "children_u5"].astype(float)
        >
        frame.loc[both, "children_u18"].astype(float)
    ).any():
        raise ValueError("Impossible child-count ordering")

    return {
        "validated": True,
        "n": int(len(frame)),
        "source_counts": {
            z: int((frame.regime == z).sum())
            for z in REGIMES
        }
    }


def preprocessor():
    numeric = Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler())
    ])

    categorical = Pipeline([
        (
            "impute",
            SimpleImputer(
                strategy="constant",
                fill_value="MISSING"
            )
        ),
        (
            "onehot",
            OneHotEncoder(
                categories=[
                    list(LEVELS[field])
                    for field in CATEGORICAL
                ],
                handle_unknown="error",
                sparse_output=False
            )
        )
    ])

    return ColumnTransformer(
        [
            (
                "numeric",
                numeric,
                list(NUMERIC + INDICATORS)
            ),
            (
                "categorical",
                categorical,
                list(CATEGORICAL)
            )
        ],
        sparse_threshold=0
    )


def aipw_scores(observed, y, e, m):
    z = np.asarray(observed)
    y = np.asarray(y, dtype=float)
    e = np.asarray(e, dtype=float)
    m = np.asarray(m, dtype=float)

    n = len(y)
    k = len(REGIMES)

    if (
        z.shape != (n,)
        or e.shape != (n, k)
        or m.shape != (n, k)
    ):
        raise ValueError("Invalid estimator matrix shapes")

    if not (
        np.isfinite(y).all()
        and np.isfinite(e).all()
        and np.isfinite(m).all()
    ):
        raise ValueError("Nonfinite model inputs")

    if (e <= 1e-6).any() or (e > 1).any():
        raise ValueError(
            "Propensity outside allowed numerical range"
        )

    if not np.allclose(
        e.sum(axis=1),
        1.0,
        atol=1e-8,
        rtol=0
    ):
        raise ValueError(
            "Propensity probabilities do not sum to one"
        )

    if ((m < -1e-8) | (m > 100 + 1e-8)).any():
        raise ValueError(
            "Outcome prediction outside 0–100"
        )

    indicator = np.column_stack([
        z == label for label in REGIMES
    ]).astype(float)

    if not np.all(indicator.sum(axis=1) == 1):
        raise ValueError("Unknown regime")

    return (
        m
        + indicator * (y[:, None] - m) / e
    )


def capped_scores(observed, y, e, m, percentile):
    if percentile not in (95, 99):
        raise ValueError("Only approved caps supported")

    # Perform all uncapped input validations.
    aipw_scores(observed, y, e, m)

    z = np.asarray(observed)
    y = np.asarray(y, dtype=float)
    e = np.asarray(e, dtype=float)
    m = np.asarray(m, dtype=float)

    result = m.copy()
    caps = {}

    for j, label in enumerate(REGIMES):
        idx = np.flatnonzero(z == label)
        weights = 1.0 / e[idx, j]

        limit = float(
            np.percentile(weights, percentile)
        )

        result[idx, j] += (
            np.minimum(weights, limit)
            * (y[idx] - m[idx, j])
        )

        caps[label] = limit

    return result, caps


def summarize(scores):
    scores = np.asarray(scores, dtype=float)

    if (
        scores.ndim != 2
        or scores.shape[1] != len(REGIMES)
        or scores.shape[0] < 2
        or not np.isfinite(scores).all()
    ):
        raise ValueError("Invalid score matrix")

    n = len(scores)
    mean = scores.mean(axis=0)

    centered = scores - mean

    covariance = (
        centered.T @ centered
        / (n * (n - 1))
    )

    difference = mean[0] - mean[1]

    difference_variance = (
        covariance[0, 0]
        + covariance[1, 1]
        - 2 * covariance[0, 1]
    )

    return {
        "mean": {
            label: float(mean[j])
            for j, label in enumerate(REGIMES)
        },
        "se": {
            label: float(
                np.sqrt(max(0, covariance[j, j]))
            )
            for j, label in enumerate(REGIMES)
        },
        "covariance": covariance,
        "contrast_piped_minus_groundwater":
            float(difference),
        "contrast_se": float(
            np.sqrt(max(0, difference_variance))
        )
    }


def fit_crossfitted(frame, *, seed=SEED, folds=5, return_internal_nuisance=False):
    validate(frame)

    if folds != 5:
        raise ValueError("G1 specifies five folds")

    if frame.regime.value_counts().min() < folds:
        raise ValueError(
            "Insufficient observations per regime"
        )

    X = frame[list(FEATURES)].copy()

    # Fixed missingness indicators.
    # Only each record's missingness is used here.
    for field in NUMERIC:
        X[field + "_missing"] = (
            X[field].isna().astype(float)
        )

    y = frame.wsqs.to_numpy(dtype=float)
    z = frame.regime.to_numpy(dtype=str)

    n = len(frame)
    k = len(REGIMES)

    e = np.full((n, k), np.nan)
    m = np.full((n, k), np.nan)
    fold_ids = np.full(n, -1, dtype=int)

    splitter = StratifiedKFold(
        n_splits=folds,
        shuffle=True,
        random_state=seed
    )

    overlaps = []

    for fold, (train, test) in enumerate(
        splitter.split(X, z)
    ):
        overlap = np.intersect1d(train, test).size
        overlaps.append(int(overlap))

        if overlap:
            raise AssertionError(
                "Training/validation leakage"
            )

        fold_ids[test] = fold

        pre = preprocessor()
        Xt = pre.fit_transform(X.iloc[train])
        Xh = pre.transform(X.iloc[test])

        propensity = LogisticRegression(
            C=1.0,
            solver="lbfgs",
            max_iter=3000,
            tol=1e-6
        )

        propensity.fit(Xt, z[train])

        probabilities = propensity.predict_proba(Xh)

        for j, label in enumerate(REGIMES):
            matched = np.flatnonzero(
                propensity.classes_ == label
            )

            if len(matched) != 1:
                raise ValueError(
                    "Training fold missing regime"
                )

            e[test, j] = probabilities[
                :, int(matched[0])
            ]

            subset = z[train] == label

            if not np.any(subset):
                raise ValueError(
                    "No outcome observations for regime"
                )

            outcome = HistGradientBoostingRegressor(
                loss="squared_error",
                learning_rate=0.05,
                max_iter=200,
                max_leaf_nodes=7,
                max_depth=3,
                min_samples_leaf=10,
                l2_regularization=1.0,
                early_stopping=False,
                random_state=seed
            )

            outcome.fit(
                Xt[subset],
                y[train][subset]
            )

            m[test, j] = outcome.predict(Xh)

    if (
        not np.all(fold_ids >= 0)
        or not np.isfinite(e).all()
        or not np.isfinite(m).all()
    ):
        raise ValueError("Incomplete cross-fitting")

    raw = aipw_scores(z, y, e, m)

    cap99, limits99 = capped_scores(
        z, y, e, m, 99
    )

    cap95, limits95 = capped_scores(
        z, y, e, m, 95
    )

    diagnostics = {}

    for j, label in enumerate(REGIMES):
        observed = z == label
        weight = 1.0 / e[observed, j]

        count = int(observed.sum())

        ess = float(
            weight.sum() ** 2
            / np.sum(weight ** 2)
        )

        concentration = float(
            weight.max() / weight.sum()
        )

        p01 = float(
            np.percentile(e[:, j], 1)
        )

        severe = (
            ess < 30
            or ess / count < 0.25
            or concentration > 0.10
            or p01 < 0.001
        )

        diagnostics[label] = {
            "n": count,
            "ess": ess,
            "ess_fraction": ess / count,
            "largest_weight_fraction":
                concentration,
            "p01_target": p01,
            "severe_support_flag": bool(severe)
        }

    result = {
        "estimate": summarize(raw),
        "capped99": summarize(cap99)["mean"],
        "capped95": summarize(cap95)["mean"],
        "caps99": limits99,
        "caps95": limits95,
        "observed": {
            label: float(y[z == label].mean())
            for label in REGIMES
        },
        "diagnostics": diagnostics,
        "source_counts": {
            label: int((z == label).sum())
            for label in REGIMES
        },
        "folds": folds,
        "fold_train_test_overlap": overlaps,
        "status":
            "SYNTHETIC_DEMONSTRATION_ONLY_NOT_EMPIRICALLY_CERTIFIED"
    }
    if return_internal_nuisance:
        result["_internal_nuisance"] = {
            "e": e,
            "m": m,
            "fold_ids": fold_ids
        }
    return result
