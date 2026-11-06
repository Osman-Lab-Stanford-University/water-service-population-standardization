"""Independent synthetic generator. No empirical data access."""

import numpy as np
import pandas as pd

from .core import REGIMES, SEED


def generate(n_per_group=10):
    if n_per_group != 10:
        raise ValueError(
            "Public example requires exactly 50 records"
        )

    rng = np.random.default_rng(SEED)
    n = 50

    regimes = np.repeat(REGIMES, 10)
    rng.shuffle(regimes)

    income = rng.integers(
        1, 6, n
    ).astype(float)

    vehicles = rng.integers(
        0, 3, n
    ).astype(float)

    household_size = rng.integers(
        1, 9, n
    ).astype(float)

    children_u18 = np.array([
        rng.integers(0, int(h) + 1)
        for h in household_size
    ], dtype=float)

    children_u5 = np.array([
        rng.integers(0, int(c) + 1)
        for c in children_u18
    ], dtype=float)

    for arr in (children_u18, children_u5):
        arr[[1, 13, 24, 38]] = np.nan

    dwelling = rng.choice(
        ["HOUSE", "APARTMENT", "TEMPORARY", "OTHER"],
        n
    )

    tenure = rng.choice(
        ["RENT", "OWN_CLEAR", "OWN_LOAN", "OTHER"],
        n
    )

    upbringing = rng.choice(
        ["URBAN", "SUBURBAN", "RURAL"],
        n
    )

    # Arbitrary fictional effects, not study estimates.
    effects = dict(zip(
        REGIMES,
        [2.0, 0.0, -3.0, -5.0, 1.0]
    ))

    quality = (
        68
        + np.array([effects[z] for z in regimes])
        + 1.4 * (income - 3)
        + rng.normal(0, 8, n)
    )

    return pd.DataFrame({
        "regime": regimes,
        "wsqs": np.clip(
            np.round(quality, 3), 0, 100
        ),
        "income": income,
        "vehicles": vehicles,
        "household_size": household_size,
        "children_u18": children_u18,
        "children_u5": children_u5,
        "dwelling": dwelling,
        "tenure": tenure,
        "upbringing": upbringing
    })
