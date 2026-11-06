"""Restricted empirical propensity diagnostic.
No household-level output, models saved or outcome estimation.
"""

import ast
import collections
import csv
import hashlib
import json
import math
import posixpath
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import log_loss

source, prior_file, raw, ref, package, output = map(
    Path, sys.argv[1:7]
)

sys.path.insert(0, str(package / "src"))

from water_standardize.support_policy import (
    assess_regime, LIMITS
)

prior = json.loads(prior_file.read_text())

# Reuse the previously verified original-survey reader
# and decoders, without executing its entire audit.
tree = ast.parse(source.read_text())

constants = {
    "SOURCE", "DOM", "CODES", "WORDNUM", "ZEROS"
}

functions = {
    "normalize", "empty", "integer", "sha",
    "qid", "col", "workbook_rows", "decode"
}

selected_nodes = []

for node in tree.body:
    if isinstance(node, ast.FunctionDef):
        if node.name in functions:
            selected_nodes.append(node)
    elif isinstance(node, ast.Assign):
        if any(
            isinstance(t, ast.Name) and t.id in constants
            for t in node.targets
        ):
            selected_nodes.append(node)

found_functions = {
    n.name for n in selected_nodes
    if isinstance(n, ast.FunctionDef)
}

if found_functions != functions:
    raise RuntimeError("Verified decoder interface changed")

W = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
R = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}'
P = '{http://schemas.openxmlformats.org/package/2006/relationships}'

ns = {
    "re": re, "math": math, "hashlib": hashlib,
    "zipfile": zipfile, "posixpath": posixpath,
    "ET": ET, "W": W, "R": R, "P": P
}

code = compile(
    ast.Module(body=selected_nodes, type_ignores=[]),
    str(source),
    "exec"
)
exec(code, ns)

# Reconstruct the frozen lookup exactly as in G7.
CORE = (
    "Q2.2", "Q3.1", "Q4.1", "Q5.1",
    "Q6.1", "Q7.1", "Q8.1", "Q9.1"
)

QIDS = (
    "Q15.17", "Q15.13", "Q15.15",
    "Q15.14", "Q15.10", "Q15.5"
)

crosswalk = ref / "questionbook/core8_response_coding.csv"

with crosswalk.open(
    encoding="utf-8-sig", newline=""
) as f:
    entries = list(csv.DictReader(f))

if len(entries) != 43:
    raise RuntimeError("Frozen crosswalk changed")

severities = {q: {} for q in CORE}

for r in entries:
    q = r["Question_ID"]
    label = ns["normalize"](r["Original_response_label"])

    if q not in severities or label in severities[q]:
        raise RuntimeError("Crosswalk ambiguity")

    severities[q][label] = (
        int(r["Final_severity_1_best_5_worst"])
        if r["Response_classification"]
        == "Scored substantive response"
        else None
    )

ns["severities"] = severities
ns["CORE"] = CORE
ns["labelcov"] = {
    q: {
        ns["normalize"](label): code
        for code, label in values.items()
    }
    for q, values in ns["CODES"].items()
}

decode = ns["decode"]
normalize = ns["normalize"]
empty = ns["empty"]
qid = ns["qid"]
rows = ns["workbook_rows"]
sha = ns["sha"]

if sha(raw) != prior["original_sha256_private"]:
    raise RuntimeError("Original workbook hash differs from G7")

assert decode("Q1.3", "Groundwater") == ("label", 2)
assert decode("Q1.1", "13") == ("code", 2)
assert decode("Q15.13", "None") == ("label", 11)
assert decode("Q15.10", "8 or more") == ("label", 8)

# Select the same primary question columns without
# assuming that headings contain only a bare QID.
first = []

for item in rows(raw):
    first.append(item)
    if len(first) >= 12:
        break

candidates = [
    (n, r) for n, r in first
    if (
        "Q1.1" in {qid(v) for v in r.values()}
        and "Q1.3" in {qid(v) for v in r.values()}
        and sum(
            q in {qid(v) for v in r.values()}
            for q in CORE
        ) >= 6
    )
]

if len(candidates) != 1:
    raise RuntimeError("Survey header ambiguous")

header_number, header = candidates[0]

required = ("Q1.1", "Q1.3") + CORE + QIDS

possible = {
    q: [
        j for j, value in header.items()
        if qid(value) == q
    ]
    for q in required
}

if any(not v for v in possible.values()):
    raise RuntimeError("Required question column missing")

valid_kinds = {
    "code", "label", "other_with_text",
    "count", "explicit_zero", "count_word",
    "count_with_unit", "numeric_display_label",
    "missing_status", "structural_na"
}

profiles = {
    q: {j: [0, 0] for j in matches}
    for q, matches in possible.items()
}

for n, row in rows(raw):
    if n <= header_number:
        continue
    for q, matches in possible.items():
        for j in matches:
            value = row.get(j)
            if not empty(value):
                profiles[q][j][1] += 1
                if decode(q, value)[0] in valid_kinds:
                    profiles[q][j][0] += 1

selected = {}

for q, matches in possible.items():
    if len(matches) == 1:
        selected[q] = matches[0]
        continue

    ranked = sorted(
        (
            (counts[0], counts[1], j)
            for j, counts in profiles[q].items()
        ),
        reverse=True
    )

    best, second = ranked[:2]

    if not (
        best[0] >= 20
        and best[0] >= 4 * max(1, second[0])
        and best[0] >= 0.5 * max(1, best[1])
    ):
        raise RuntimeError("Primary column ambiguity for " + q)

    selected[q] = best[2]

if len(selected) != 16:
    raise RuntimeError("Six-covariate column resolution failed")

id_columns = [
    j for j, value in header.items()
    if normalize(value) in (
        "responseid", "responseidentifier"
    )
]

if len(id_columns) != 1:
    raise RuntimeError("Response ID column ambiguous")

id_column = id_columns[0]

# Predeclared categorical codes. These are fixed
# questionnaire categories, not fitted encoders.
categories = {
    "Q15.17": (1, 2, 3, 4, 5),
    "Q15.13": (9, 10, 11),
    "Q15.15": (9, 10, 11, 5),
    "Q15.14": (1, 2, 3, 4),
    "Q15.10": (20, 21, 22, 4, 5, 6, 7, 8),
    "Q15.5": (2, 14, 15)
}

regime_names = {
    1: "PIPED",
    2: "GROUNDWATER",
    3: "SURFACE",
    4: "TANKER",
    5: "COMMERCIAL"
}

records = []
regimes = []
seen = set()
counts = collections.Counter()
missing_counts = collections.defaultdict(collections.Counter)

for n, row in rows(raw):
    if n <= header_number:
        continue

    rid = str(row.get(id_column, "") or "").strip()

    if not re.fullmatch(r"R_[A-Za-z0-9]+", rid):
        continue

    token = hashlib.sha256(rid.encode()).digest()

    if token in seen:
        raise RuntimeError("Duplicate Qualtrics response ID")

    seen.add(token)
    counts["substantive"] += 1

    source = decode(
        "Q1.3", row.get(selected["Q1.3"])
    )[1]

    if source is None or source not in range(1, 7):
        raise RuntimeError("Invalid drinking-source code")

    if source == 6:
        counts["other_source"] += 1
        continue

    scores = [
        decode(q, row.get(selected[q]))
        for q in CORE
    ]

    if any(kind == "invalid" for kind, _ in scores):
        raise RuntimeError("Invalid CORE-8 response")

    if sum(value is not None for _, value in scores) < 7:
        raise RuntimeError("WSQS eligibility mismatch")

    features = []

    for q in QIDS:
        kind, value = decode(
            q, row.get(selected[q])
        )

        if value is None and kind != "missing":
            raise RuntimeError(
                "Unrecognized primary covariate " + q
            )

        if value is not None and value not in categories[q]:
            raise RuntimeError(
                "Unexpected covariate category " + q
            )

        features.append(value)

        if value is None:
            missing_counts[regime_names[source]][q] += 1

    records.append(features)
    regimes.append(source)

if counts["substantive"] != 3165:
    raise RuntimeError("Substantive cohort mismatch")

if counts["other_source"] != 99:
    raise RuntimeError("Sixth source count mismatch")

y = np.asarray(regimes, dtype=int)

expected = {1:1275, 2:1212, 3:234, 4:120, 5:225}

if len(records) != 3066:
    raise RuntimeError("Five-regime cohort mismatch")

if any(int(np.sum(y == z)) != n for z,n in expected.items()):
    raise RuntimeError("Five-regime sample counts mismatch")

# One-hot columns are prespecified across all
# possible questionnaire categories, including missing.
# This does not use empirical outcome information.
Xparts = []

for j, q in enumerate(QIDS):
    options = list(categories[q]) + [None]

    block = np.zeros(
        (len(records), len(options)),
        dtype=float
    )

    location = {v:k for k,v in enumerate(options)}

    for i, row in enumerate(records):
        block[i, location[row[j]]] = 1.0

    Xparts.append(block)

X = np.column_stack(Xparts)

assert np.isfinite(X).all()
assert len(X) == len(y) == 3066


import importlib.util
import inspect
from water_standardize import core as engine
from water_standardize import demo

# Independently invoke the frozen measurement scorer.
scorer = ref / "src/wsqs/scoring.py"

spec = importlib.util.spec_from_file_location(
    "frozen_core8_empirical03", scorer
)

scoring_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scoring_module)

assert tuple(q for q, _ in scoring_module.CORE8) == CORE

scores = []
sources = []
domains = collections.Counter()

# Outcome reconstruction occurs in memory only.
for number, row in rows(raw):

    if number <= header_number:
        continue

    identifier = str(
        row.get(id_column, "") or ""
    ).strip()

    if not re.fullmatch(
        r"R_[A-Za-z0-9]+", identifier
    ):
        continue

    source = decode(
        "Q1.3", row.get(selected["Q1.3"])
    )[1]

    if source == 6:
        continue

    if source not in range(1, 6):
        raise RuntimeError("Invalid empirical source")

    answers = {
        q: decode(q, row.get(selected[q]))
        for q in CORE
    }

    if any(
        k == "invalid"
        for k, v in answers.values()
    ):
        raise RuntimeError("Invalid CORE-8 answer")

    sevs = {
        q: value
        for q, (kind, value) in answers.items()
    }

    count = sum(
        v is not None
        for v in sevs.values()
    )

    if count not in (7, 8):
        raise RuntimeError("Ineligible CORE-8 score")

    score = float(scoring_module.core8(sevs))

    independent = sum(
        25 * (5-v)
        for v in sevs.values()
        if v is not None
    ) / count

    if not math.isclose(
        score, independent,
        rel_tol=0,
        abs_tol=1e-12
    ):
        raise RuntimeError("CORE-8 scorer parity failed")

    if not 0 <= score <= 100:
        raise RuntimeError("WSQS outside 0–100")

    sources.append(source)
    scores.append(score)
    domains[str(count)] += 1

sources = np.asarray(sources, dtype=int)
scores = np.asarray(scores, dtype=float)

# Must match the independently reconstructed
# six-covariate cohort and order exactly.
if (
    len(scores) != 3066
    or not np.array_equal(sources, y)
):
    raise RuntimeError(
        "Outcome and propensity cohort mismatch"
    )

observed = {}

for z in range(1, 6):
    name = regime_names[z]
    values = scores[sources == z]

    if len(values) != expected[z]:
        raise RuntimeError("Regime sample mismatch")

    observed[name] = {
        "n": int(len(values)),
        "observed_mean_WSQS": float(values.mean()),
        "observed_median_WSQS": float(np.median(values)),
        "observed_sd_WSQS": float(
            np.std(values, ddof=1)
        ),
        "p10_WSQS": float(np.quantile(values, .10)),
        "p90_WSQS": float(np.quantile(values, .90)),
        "observed_risk_WSQS_le_50": float(
            np.mean(values <= 50)
        ),
        "observed_risk_WSQS_le_60": float(
            np.mean(values <= 60)
        ),
        "observed_risk_WSQS_le_70": float(
            np.mean(values <= 70)
        )
    }

SOURCE_WORKBOOK = raw

import warnings
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import log_loss

from water_standardize import core as engine
from water_standardize.support_policy import assess_regime

# Empirical 03 already reconstructed and reconciled
# the outcome, the six covariates and all 3066 records.

reference02 = json.loads(
    (output / "_empirical02_reference.json").read_text()
)

reference03 = json.loads(
    (output / "_empirical03_reference.json").read_text()
)

if reference02["n"] != 3066 or reference03["sample_n"] != 3066:
    raise RuntimeError("Prior sample counts disagree")

if not reference02["previous_HGBR_exactly_reproduced"]:
    raise RuntimeError("Empirical 02 provenance failed")

engine_file = package / "src/water_standardize/core.py"

if sha(engine_file) != reference03[
    "engine_information"
]["engine_core_sha256"]:
    raise RuntimeError("Frozen engine source hash changed")

labels = list(engine.REGIMES)

if labels != [
    "PIPED", "GROUNDWATER", "SURFACE",
    "TANKER", "COMMERCIAL"
]:
    raise RuntimeError("Unexpected frozen regime order")

seed = 20261022
n = len(scores)
k = len(labels)

if (
    n != 3066 or len(X) != n
    or not np.array_equal(sources, y)
    or not np.isfinite(scores).all()
):
    raise RuntimeError("Outcome-covariate cohort mismatch")

z = np.asarray(
    [regime_names[int(v)] for v in y],
    dtype=str
)

e = np.full((n,k), np.nan)
m = np.full((n,k), np.nan)
fold_ids = np.full(n, -1, dtype=int)
fold_counts = []

splits = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=seed
)

for fold, (train, test) in enumerate(
    splits.split(X, y)
):
    if np.intersect1d(train, test).size:
        raise RuntimeError("Cross-fitting leakage")

    fold_ids[test] = fold

    propensity = LogisticRegression(
        penalty="l2",
        C=1.0,
        solver="lbfgs",
        max_iter=2000,
        tol=1e-6,
        fit_intercept=True
    )

    with warnings.catch_warnings():
        warnings.simplefilter(
            "error", ConvergenceWarning
        )
        propensity.fit(X[train], y[train])

    if not np.array_equal(
        propensity.classes_, np.arange(1,6)
    ):
        raise RuntimeError("Propensity class-order mismatch")

    e[test] = propensity.predict_proba(X[test])

    fold_counts.append({
        labels[j]: int(np.sum(y[test] == j+1))
        for j in range(k)
    })

    for j, label in enumerate(labels):

        subset = y[train] == j+1

        if subset.sum() < 5:
            raise RuntimeError(
                "Insufficient outcome training data"
            )

        # Exact learner hyperparameters from the
        # verified frozen estimator implementation.
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
            X[train][subset],
            scores[train][subset]
        )

        m[test,j] = outcome.predict(X[test])

if (
    np.any(fold_ids < 0)
    or not np.isfinite(e).all()
    or not np.isfinite(m).all()
    or np.any(e <= 0)
    or np.any(e >= 1)
    or not np.allclose(e.sum(axis=1),1,atol=1e-9)
):
    raise RuntimeError("Cross-fitting validity failure")

# Exact reproduction of the propensity specification
# independently executed in Empirical 02.
ll = float(log_loss(
    y,e,labels=[1,2,3,4,5]
))

prior_logloss = reference02[
    "model_results"
]["LOGISTIC_L2_C1"]["out_of_fold_logloss"]

if not math.isclose(
    ll,prior_logloss,
    rel_tol=0,abs_tol=1e-9
):
    raise RuntimeError(
        "Empirical 02 propensity log loss did not reproduce"
    )

diagnostics = {}

for j,label in enumerate(labels):
    obs = y == j+1
    w = 1.0 / e[obs,j]

    ess = float(w.sum()**2 / np.square(w).sum())

    d = {
        "n":int(obs.sum()),
        "ess":ess,
        "ess_fraction":float(ess/obs.sum()),
        "largest_weight_fraction":
            float(w.max()/w.sum()),
        "p01_target":float(np.quantile(e[:,j],.01))
    }

    screen = assess_regime(d)

    prior_diag = reference02[
        "model_results"
    ]["LOGISTIC_L2_C1"]["regimes"][label]

    for field in (
        "n","ess","ess_fraction",
        "largest_weight_fraction","p01_target"
    ):
        if not math.isclose(
            d[field],prior_diag[field],
            rel_tol=1e-9,
            abs_tol=1e-9
        ):
            raise RuntimeError(
                "Empirical 02 diagnostic parity failed: "
                +label+" "+field
            )

    weighted = np.average(
        X[obs],axis=0,weights=w
    )

    balance = float(np.max(np.abs(
        weighted-X.mean(axis=0)
    )))

    if not math.isclose(
        balance,
        prior_diag["max_absolute_category_balance_gap"],
        rel_tol=0,
        abs_tol=1e-9
    ):
        raise RuntimeError(
            "Empirical 02 balance parity failed"
        )

    diagnostics[label] = {
        **d,
        "balance_max_absolute_category_gap":balance,
        "screen_pass":bool(screen["screen_pass"]),
        "screen_reasons":list(screen["reasons"])
    }

if not all(
    d["screen_pass"] for d in diagnostics.values()
):
    raise RuntimeError(
        "G5.5 support screen failed; "
        "AIPW estimates not released"
    )

# USE THE EXISTING FROZEN AIPW SCORE FUNCTION.
raw = np.asarray(
    engine.aipw_scores(z,scores,e,m),
    dtype=float
)

if raw.shape != (n,k) or not np.isfinite(raw).all():
    raise RuntimeError("Frozen AIPW score output invalid")

# Independent numerical verification of its estimating
# equation. This is an audit, not a replacement estimator.
independent = np.empty_like(raw)

for j,label in enumerate(labels):
    independent[:,j] = (
        m[:,j]
        + (z == label).astype(float)
        * (scores-m[:,j]) / e[:,j]
    )

max_discrepancy = float(
    np.max(np.abs(raw-independent))
)

if max_discrepancy > 1e-8:
    raise RuntimeError(
        "Frozen AIPW estimating equation did not reconcile"
    )

summary = engine.summarize(raw)

def ordered_values(value):
    if isinstance(value,dict):
        if set(value)!=set(labels):
            raise RuntimeError("Unexpected regime keys")
        a=np.asarray(
            [value[key] for key in labels],dtype=float
        )
    else:
        a=np.asarray(value,dtype=float)
    if a.shape!=(k,) or not np.isfinite(a).all():
        raise RuntimeError("Unexpected estimator summary")
    return a

if not isinstance(summary,dict) or "mean" not in summary:
    raise RuntimeError("Frozen summary interface changed")

means = ordered_values(summary["mean"])

# Reconcile frozen summary with mean household-level
# AIPW scores. No records are retained after this run.
direct = raw.mean(axis=0)

if not np.allclose(means,direct,atol=1e-9,rtol=0):
    raise RuntimeError("Frozen AIPW summary mismatch")

# Apply the same frozen percentile-capping operations
# used in the synthetic estimator.
cap99,limits99 = engine.capped_scores(
    z,scores,e,m,99
)

cap95,limits95 = engine.capped_scores(
    z,scores,e,m,95
)

mean99 = ordered_values(
    engine.summarize(cap99)["mean"]
)

mean95 = ordered_values(
    engine.summarize(cap95)["mean"]
)

# Influence-score SE is recorded as a diagnostic, not
# treated as publication-certified uncertainty.
se = raw.std(axis=0,ddof=1)/np.sqrt(n)

result = {}

for j,label in enumerate(labels):

    observed_mean = float(scores[z==label].mean())

    prior_mean = reference03[
        "observed_by_regime"
    ][label]["observed_mean_WSQS"]

    if not math.isclose(
        observed_mean,prior_mean,
        abs_tol=1e-10,rel_tol=0
    ):
        raise RuntimeError(
            "Observed Empirical 03 means did not reproduce"
        )

    predicted_all = float(m[:,j].mean())

    augmentation = float(
        (raw[:,j]-m[:,j]).mean()
    )

    if not math.isclose(
        predicted_all+augmentation,
        means[j],abs_tol=1e-9,rel_tol=0
    ):
        raise RuntimeError(
            "AIPW decomposition mismatch"
        )

    actual = z==label

    result[label] = {
        "n":int(actual.sum()),
        "observed_mean":observed_mean,
        "standardized_AIPW_mean":float(means[j]),
        "standardized_minus_observed":
            float(means[j]-observed_mean),
        "outcome_model_average":predicted_all,
        "AIPW_augmentation":augmentation,
        "capped99_mean":float(mean99[j]),
        "capped95_mean":float(mean95[j]),
        "diagnostic_influence_SE":float(se[j]),
        "outcome_OOF_RMSE_observed_group":
            float(np.sqrt(np.mean(
                (scores[actual]-m[actual,j])**2
            )))
    }

i=labels.index("PIPED")
j=labels.index("GROUNDWATER")

contrast_scores = raw[:,i]-raw[:,j]

contrast = {
    "direction":"PIPED_MINUS_GROUNDWATER",
    "observed_difference":
        float(result["PIPED"]["observed_mean"]
              -result["GROUNDWATER"]["observed_mean"]),
    "standardized_AIPW_difference":
        float(means[i]-means[j]),
    "capped99_difference":
        float(mean99[i]-mean99[j]),
    "capped95_difference":
        float(mean95[i]-mean95[j]),
    "diagnostic_influence_SE":
        float(np.std(
            contrast_scores,ddof=1
        )/np.sqrt(n))
}

report = {
    "phase":"EMPIRICAL_04_AIPW_STANDARDIZATION",
    "status":"AGGREGATE_EMPIRICAL_AIPW_ESTIMATED",
    "analytical_n":n,
    "target_distribution":
        "Empirical six-covariate distribution, n=3066",
    "provision_variable":
        "Q1.3 primary drinking-water source",
    "six_adjustment_question_ids":list(QIDS),
    "propensity":
        "Five-fold L2 logistic C=1.0, "
        "Empirical 02 preprocessing and exact parity",
    "outcome_learner":
        "Five-fold HGBR, frozen engine hyperparameters",
    "feature_encoding":
        "Predeclared six-covariate one-hot categories, "
        "including explicit missing category",
    "frozen_engine_unchanged":True,
    "frozen_engine_sha256":sha(engine_file),
    "engine_score_function_reused":True,
    "engine_summary_function_reused":True,
    "engine_capping_functions_reused":True,
    "full_frozen_fit_crossfitted_called":False,
    "reason_for_adapter":
        "Original frozen fitting function uses eight "
        "covariates; primary specification uses six",
    "scoring_parity":True,
    "propensity_02_exact_parity":True,
    "observed_03_exact_parity":True,
    "aipw_equation_max_difference":
        max_discrepancy,
    "propensity_out_of_fold_logloss":ll,
    "overlap":diagnostics,
    "fold_regime_counts":fold_counts,
    "by_regime":result,
    "main_contrast":contrast,
    "diagnostic_SE_not_publication_certified":True,
    "full_inference_certified":False,
    "distributional_AIPW_not_yet_estimated":True,
    "household_records_saved":False,
    "propensity_predictions_saved":False,
    "outcome_predictions_saved":False,
    "model_objects_saved":False,
    "public_release_authorized":False
}

(output/"EMPIRICAL_04_REPORT.json").write_text(
    json.dumps(
        report,indent=2,sort_keys=True,allow_nan=False
    )+"\n",
    encoding="utf-8"
)

if sha(SOURCE_WORKBOOK) != prior["original_sha256_private"]:
    raise RuntimeError("Original source changed")

print("ANALYTICAL_SAMPLE:",n)
print("FROZEN_WSQS_SCORING: VERIFIED")
print("EMPIRICAL_02_PROPENSITY_PARITY: PASS")
print("EMPIRICAL_03_OBSERVED_MEANS_PARITY: PASS")
print("FROZEN_AIPW_SCORE_EQUATION_PARITY: PASS")
print("FIVE_G5_5_SUPPORT_SCREENS: PASS")
print("LOGISTIC_PROPENSITY_LOGLOSS:",round(ll,5))
print("")
print("REGIME | N | OBSERVED | STANDARDIZED | CHANGE | CAP99 | CAP95 | DIAGNOSTIC_SE")

for label in labels:
    d=result[label]
    print(
        label,
        d["n"],
        round(d["observed_mean"],3),
        round(d["standardized_AIPW_mean"],3),
        round(d["standardized_minus_observed"],3),
        round(d["capped99_mean"],3),
        round(d["capped95_mean"],3),
        round(d["diagnostic_influence_SE"],3)
    )

print("")
print(
    "PIPED_MINUS_GROUNDWATER_OBSERVED:",
    round(contrast["observed_difference"],4)
)
print(
    "PIPED_MINUS_GROUNDWATER_STANDARDIZED:",
    round(contrast["standardized_AIPW_difference"],4)
)
print(
    "PIPED_MINUS_GROUNDWATER_DIAGNOSTIC_SE:",
    round(contrast["diagnostic_influence_SE"],4)
)
print(
    "PIPED_MINUS_GROUNDWATER_CAP99:",
    round(contrast["capped99_difference"],4)
)
print(
    "PIPED_MINUS_GROUNDWATER_CAP95:",
    round(contrast["capped95_difference"],4)
)
print("")
print("STANDARDIZED_TAIL_RISK: NOT YET ESTIMATED")
print("PUBLICATION_CONFIDENCE_INTERVALS: NOT YET CERTIFIED")
print("RESPONDENT_LEVEL_FILES: NONE")
print("PUBLIC_RELEASE_AUTHORIZED: NO")
print("EMPIRICAL_04_RESULT: COMPLETE")
