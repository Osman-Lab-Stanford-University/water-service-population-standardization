import json
import hashlib
import math
import re
import sys
from pathlib import Path

root, out = map(
    lambda x: Path(x).resolve(),
    sys.argv[1:]
)

assert root.is_dir()
assert out.parent.parent == root

specs = {
    "04": ("AIPW", "EMPIRICAL_04"),
    "08": ("LOGISTIC_TAIL", "EMPIRICAL_08"),
    "13": ("CHILD_RECONCILIATION", "EMPIRICAL_13"),
    "14D": ("LOG_CHILDREN", "EMPIRICAL_14D"),
    "14E": ("EXTREME_COUNT", "EMPIRICAL_14E"),
    "14F": ("ARGFIX", "EMPIRICAL_14F"),
    "15": ("PRIMARY_FULL_CDF", "EMPIRICAL_15"),
    "16": ("SIMULTANEOUS_CDF", "EMPIRICAL_16"),
    "17": ("ORACLE_CDF_COVERAGE", "EMPIRICAL_17"),
    "18": ("FITTED_NUISANCE_COVERAGE", "EMPIRICAL_18"),
    "19": ("TAIL_COVERAGE_FORENSICS", "EMPIRICAL_19"),
    "20": ("INDEPENDENT_CDF_VALIDATION", "EMPIRICAL_20"),
    "21": ("ROW3_CORRECTED", "EMPIRICAL_21"),
    "22": ("TIMESTAMP_CORRECTED", "EMPIRICAL_22")
}

manifest_names = {
    "14F": "EMPIRICAL_14F_ARGFIX_MANIFEST.sha256",
    "21": "EMPIRICAL_21_ROW3_MANIFEST.sha256",
    "22": "EMPIRICAL_22_TIMESTAMP_MANIFEST.sha256"
}

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(
            lambda: f.read(2**20), b""
        ):
            h.update(part)
    return h.hexdigest()

def load(k, suffix, prefix):
    dirs = [
        d for d in root.glob(
            prefix + "_" + suffix + "_*"
        )
        if (
            d / "PRIVATE_VERIFICATION"
            / ("EMPIRICAL_" + k + "_REPORT.json")
        ).is_file()
    ]

    if not dirs:
        raise RuntimeError(
            "Missing frozen stage " + k
        )

    d = max(
        dirs,
        key=lambda p: p.stat().st_mtime_ns
    )

    v = d / "PRIVATE_VERIFICATION"

    mf = v / manifest_names.get(
        k, prefix + "_MANIFEST.sha256"
    )

    if not mf.is_file():
        raise RuntimeError(
            "Missing manifest " + k
        )

    protected = set()

    for line in mf.read_text().splitlines():
        if not line.strip():
            continue

        m = re.fullmatch(
            r"([0-9a-fA-F]{64})\s+\*?(.+)",
            line
        )

        if not m:
            raise RuntimeError(
                "Malformed manifest " + k
            )

        p = (v / m[2]).resolve()

        if (
            not p.is_relative_to(d.resolve())
            or not p.is_file()
            or sha(p) != m[1].lower()
        ):
            raise RuntimeError(
                "Manifest mismatch " + k
            )

        protected.add(p)

    report = v / (prefix + "_REPORT.json")

    if report.resolve() not in protected:
        raise RuntimeError(
            "Report not protected " + k
        )

    data = json.loads(report.read_text())

    assert data["public_release_authorized"] is False

    return data, {
        "stage": d.name,
        "report_sha256": sha(report),
        "manifest_sha256": sha(mf)
    }

reports = {}
provenance = {}

for k, (suffix, prefix) in specs.items():
    reports[k], provenance[k] = load(
        k, suffix, prefix
    )

a, b, c, d, e, f, g, h = [
    reports[k]
    for k in (
        "04", "08", "13", "14D",
        "14E", "14F", "15", "16"
    )
]

def same(x, y, tol=1e-8):
    if (
        not math.isfinite(float(x))
        or not math.isfinite(float(y))
        or abs(float(x)-float(y)) > tol
    ):
        raise RuntimeError(
            "Results discrepancy: "
            + str(x) + " versus " + str(y)
        )

labels = [
    "PIPED",
    "GROUNDWATER",
    "SURFACE",
    "TANKER",
    "COMMERCIAL"
]

assert all(
    reports[k].get("analytical_n") == 3066
    for k in (
        "04", "13", "14D",
        "14E", "14F", "15", "16"
    )
)

assert g["regimes"] == labels
assert len(g["threshold_rows"]) == 41
assert len(h["rows"]) == 41
assert g["grid"]["increment"] == 2.5
assert f["G5_5_all_regimes_pass"] is True

assert c["overall_counts"][
    "both_usable_after_consistency"
] == 2532

assert c["overall_counts"][
    "at_least_one_missing_after_consistency"
] == 534

assert reports["21"][
    "original_export_unique_response_ids"
] == 3165

assert reports["22"][
    "original_unique_response_ids"
] == 3165

assert reports["22"][
    "protocol_verified_session_identifiers"
] is False

assert reports["21"][
    "verified_recruitment_session_identifiers"
] == 0

assert len(
    reports["22"]["design_effect_scenarios"]
) == 12

assert h["e15_exact_estimate_parity_max_pp"] <= 1e-8
assert h["e15_exact_se_parity_max_pp"] <= 1e-8

assert not h["recruitment_clustering_accounted_for"]
assert not h["monte_carlo_coverage_certified"]

# Reconcile standardized mean estimates.

means = []

for label in labels:
    m = float(
        a["by_regime"][label][
            "standardized_AIPW_mean"
        ]
    )

    six = d["mean_results"][label]
    cap = f["capped_sensitivity"]["MEAN"][label]

    same(m, six["six_covariate"])
    same(cap["uncapped"], six["eight_covariate"])

    support = f["G5_5_regime_results"][label]

    assert support["G5_5_screen_pass"]
    assert support["n"] == six["n"]

    means.append({
        "group": label,
        "n": six["n"],
        "primary_mean": m,
        "eight_covariate_mean":
            six["eight_covariate"],
        "eight_ESS": support["ess"],
        "cap99_mean_change":
            cap["cap99"] - cap["uncapped"],
        "cap95_mean_change":
            cap["cap95"] - cap["uncapped"]
    })

assert sum(x["n"] for x in means) == 3066

main = float(
    a["main_contrast"][
        "standardized_AIPW_difference"
    ]
)

same(
    main,
    means[0]["primary_mean"]
    - means[1]["primary_mean"]
)

# Reconcile every raw CDF estimate with the
# independently computed Empirical 16 outputs.

cdf = []
risks = []

for i, (r, s) in enumerate(
    zip(g["threshold_rows"], h["rows"])
):
    threshold = i * 2.5

    same(r["threshold_wsqs"], threshold)
    same(s["threshold_wsqs"], threshold)

    raw = {
        label: float(r["raw_aipw_cdf"][label])
        for label in labels
    }

    if any(
        not 0 <= value <= 1
        for value in raw.values()
    ):
        raise RuntimeError(
            "Out-of-bounds CDF probability"
        )

    diff = 100 * (
        raw["PIPED"] - raw["GROUNDWATER"]
    )

    same(
        diff,
        s["contrast_piped_minus_groundwater_pp"]
    )

    same(
        diff,
        r["raw_piped_minus_groundwater_pp"]
    )

    same(
        s["working_independent_household_se_pp"],
        r["working_independent_household_SE_pp"]
    )

    entry = {
        "threshold": threshold,
        "regime_cdf_percent": {
            label: 100*value
            for label, value in raw.items()
        },
        "contrast_pp": diff,
        "working_simultaneous_95_ci_pp": [
            s["working_simultaneous_95_lower_pp"],
            s["working_simultaneous_95_upper_pp"]
        ],
        "working_familywise_adjusted_p":
            s["familywise_multiplier_adjusted_p"]
    }

    cdf.append(entry)

    if threshold in (50, 60, 70):
        anchor = b["results"][
            str(int(threshold))
        ]

        for label in labels:
            same(
                raw[label],
                anchor["by_regime"][label][
                    "standardized_AIPW_probability"
                ],
                0.0015
            )

        risks.append(entry)

# Confirm monotonicity without modifying values.

for label in labels:
    values = [
        r["regime_cdf_percent"][label]
        for r in cdf
    ]

    assert all(
        values[i+1] >= values[i]-1e-8
        for i in range(40)
    )

    assert g["cdf_diagnostics"][label][
        "raw_steps_decreasing"
    ] == 0

maxrow = max(
    cdf,
    key=lambda x: abs(x["contrast_pp"])
)

same(
    maxrow["threshold"],
    h["maximum_absolute_gap_threshold_wsqs"]
)

same(
    maxrow["contrast_pp"],
    g["max_abs_piped_minus_groundwater_difference_pp"]
)

assert len(risks) == 3

for k in ("17", "18", "19", "20"):
    assert reports[k][
        "publication_inference_certified"
    ] is False

effect = (
    e["mean_piped_minus_groundwater"][
        "eight_covariate"
    ]
    - d["mean_piped_minus_groundwater"][
        "eight_covariate"
    ]
)

assert e["extreme_count_stress_test"][
    "affected_households"
] == 1

register = {
    "stage": "EMPIRICAL_23",
    "status":
        "RECONCILED_NOT_PUBLICATION_CERTIFIED",
    "n": 3066,
    "primary_mean_results": means,
    "primary_piped_minus_groundwater_mean": main,
    "registered_threshold_risks": risks,
    "full_41_threshold_CDF": cdf,
    "working_global_CDF_p":
        h["working_global_no_cdf_difference_p"],
    "working_simultaneous_significant_thresholds":
        len(h["simultaneously_significant_thresholds"]),
    "eight_covariate_piped_minus_groundwater_mean":
        d["mean_piped_minus_groundwater"][
            "eight_covariate"
        ],
    "extreme_count_stress_mean_contrast_change":
        effect,
    "parent_export_unique_response_ids": 3165,
    "verified_recruitment_sessions": 0,
    "synthetic_validation_reports_verified":
        [17, 18, 19, 20],
    "inferential_cautions": [
        "Recruitment dependence unresolved",
        "Full-grid fitted-nuisance coverage not certified",
        "Sparse-regime full-grid studentized bands undercovered",
        "Nonprobability sample; noncausal comparisons"
    ],
    "provenance": provenance,
    "restricted_household_data_read": False,
    "public_release_authorized": False,
    "publication_inference_certified": False
}

(out / "EMPIRICAL_23_REGISTER.json").write_text(
    json.dumps(
        register,
        sort_keys=True,
        indent=2,
        allow_nan=False
    ) + "\n"
)

# Produce manuscript-ready aggregate tables.

lines = [
    "# Empirical 23 | Verified manuscript numbers",
    "",
    "**Status:** Aggregate results verified; "
    "inference not publication-certified.",
    "",
    "## Five source-group standardized means",
    "",
    "| Group | N | Primary WSQS | Eight-covariate WSQS |",
    "|---|---:|---:|---:|"
]

for row in means:
    lines.append(
        f"| {row['group']} | {row['n']:,} | "
        f"{row['primary_mean']:.3f} | "
        f"{row['eight_covariate_mean']:.3f} |"
    )

lines += [
    "",
    f"Primary piped minus groundwater: "
    f"**{main:+.3f} WSQS points**.",
    "",
    "## Registered lower-tail risks",
    "",
    "| Threshold | Piped | Groundwater | Surface | "
    "Tanker | Commercial | Piped − groundwater |",
    "|---|---:|---:|---:|---:|---:|---:|"
]

for r in risks:
    values = r["regime_cdf_percent"]

    lines.append(
        f"| ≤{r['threshold']:g} | "
        + " | ".join(
            f"{values[label]:.3f}%"
            for label in labels
        )
        + f" | {r['contrast_pp']:+.3f} pp |"
    )

lines += [
    "",
    "## Full distribution and inference limitations",
    "",
    "All 41 raw CDF thresholds are bounded and monotone.",
    "",
    f"Maximum absolute piped–groundwater gap: "
    f"{maxrow['contrast_pp']:+.3f} pp at "
    f"WSQS {maxrow['threshold']:g}.",
    "",
    f"**Working** full-grid global p-value: "
    f"{h['working_global_no_cdf_difference_p']:.3f}. "
    f"Simultaneously significant thresholds: "
    f"{len(h['simultaneously_significant_thresholds'])}.",
    "",
    f"Eight-covariate mean contrast: "
    f"{d['mean_piped_minus_groundwater']['eight_covariate']:+.3f}. "
    f"Extreme-count stress change: {effect:+.4f} points.",
    "",
    "No verified recruitment sessions are available. "
    "Hypothetical ICC calculations are not cluster-robust "
    "inference. Sparse-regime full-grid studentized bands "
    "failed synthetic coverage tests.",
    "",
    "The sample is nonprobability-based. "
    "These are descriptive, not causal or "
    "population-representative comparisons.",
    "",
    f"**Provenance:** {len(provenance)} frozen reports "
    "and manifests verified without reading household records.",
    ""
]

(out / "EMPIRICAL_23_MANUSCRIPT_NUMBERS.md").write_text(
    "\n".join(lines)
)

files = sorted(
    p for p in out.iterdir()
    if p.is_file()
)

(out / "EMPIRICAL_23_MANIFEST.sha256").write_text(
    "".join(
        f"{sha(p)}  {p.name}\n"
        for p in files
    )
)

print("EMPIRICAL_23_RESULTS")
print(
    "FROZEN_REPORTS_AND_MANIFESTS_VERIFIED:",
    len(provenance)
)
print("PRIMARY_SAMPLE:", 3066)
print("SOURCE_GROUPS:", len(labels))
print(
    "MEAN_PIPED_MINUS_GROUNDWATER:",
    round(main, 5)
)

for row in means:
    print(
        "MEAN", row["group"],
        round(row["primary_mean"], 4),
        "EIGHT",
        round(row["eight_covariate_mean"], 4)
    )

for r in risks:
    print(
        "TAIL", r["threshold"],
        "PIPED_MINUS_GROUNDWATER_PP",
        round(r["contrast_pp"], 4)
    )

print("FULL_CDF_THRESHOLDS:", len(cdf))
print(
    "MAX_GAP_PP:",
    round(maxrow["contrast_pp"], 4),
    "AT_WSQS:", maxrow["threshold"]
)
print(
    "WORKING_FULL_CDF_GLOBAL_P:",
    round(
        h["working_global_no_cdf_difference_p"],
        6
    )
)
print("RECRUITMENT_SESSION_IDS_VERIFIED: 0")
print("INFERENCE_PUBLICATION_CERTIFIED: NO")
print("HOUSEHOLD_LEVEL_DATA_READ: NO")
print("EMPIRICAL_23_RESULT: COMPLETE")
