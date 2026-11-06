#!/usr/bin/env python3
"""Fail-closed public path/CSV roster and simple identifier heuristic.
This is not an IRB disclosure review or a guarantee against re-identification.
"""
from pathlib import Path
import csv,re,sys
ROOT=Path(__file__).resolve().parents[1]
EXPECTED={
 'data/aggregates/eight_covariate_sensitivity.csv':5,
 'data/aggregates/figure3_source_level_max_absolute_smd.csv':5,
 'data/aggregates/full_41_threshold_standardized_cdf_percent.csv':41,
 'data/aggregates/main_table_1_observed_standardized_means.csv':5,
 'data/aggregates/main_table_2_standardized_threshold_probabilities.csv':5,
 'data/aggregates/observed_wsqs_sd_by_source.csv':5,
 'data/aggregates/overlap_diagnostics_original_aggregate.csv':5,
 'data/aggregates/si_table_s2_questionnaire_item_mapping.csv':8,
 'data/aggregates/si_table_s3_drinking_domestic_source_comparison.csv':5,
 'data/aggregates/si_table_s4_domestic_source_sensitivity.csv':5,
 'data/aggregates/si_table_s5_four_contrast_working_intervals.csv':4,
 'data/aggregates/si_table_s6_ten_pairwise_mean_intervals.csv':10,
 'data/aggregates/si_table_s7_core_domain_sensitivity.csv':5,
 'data/aggregates/synthetic_simulation_coverage_summary.csv':4,
 'examples/reconstructed_figure_s3_pairwise_distances.csv':10,
}
actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*.csv')}
if actual!=set(EXPECTED):
    raise SystemExit('Public CSV allowlist violation: '+str(sorted(actual.symmetric_difference(EXPECTED))))
for p in ROOT.rglob('*'):
    rel=p.relative_to(ROOT).as_posix()
    if p.is_symlink(): raise SystemExit('Symlinks not permitted: '+rel)
    if not p.is_file():continue
    if any(x in rel.split('/') for x in ('.git','__pycache__','.venv','venv','.DS_Store')):raise SystemExit('Generated/private file: '+rel)
    if p.suffix.lower() in ('.xlsx','.xls','.sav','.rds','.parquet','.docx','.pdf','.zip','.pkl','.joblib','.pyc','.ipynb','.sqlite','.db','.npz','.npy','.env'):
        raise SystemExit('Unreviewed data or binary file: '+rel)
    if p.stat().st_size>2_000_000:raise SystemExit('Oversize/unreviewed public file: '+rel)
    if p.suffix.lower() in ('.md','.py','.csv','.json','.txt','.cff','.yml','.yaml') or p.name=='.gitignore':
        s=p.read_text(encoding='utf-8',errors='replace')
        # Specific user home paths and direct identifiers, not ordinary cited URLs.
        if re.search(r'(?i)(?:/Users/[A-Za-z0-9._-]{2,}/|/home/[A-Za-z0-9._-]{2,}/)',s):
            raise SystemExit('Private user path: '+rel)
        if re.search(r'(?i)\b[A-Za-z0-9._%+-]+@(?:[A-Za-z0-9-]+\.)+[A-Za-z]{2,}\b',s):
            raise SystemExit('Potential email address: '+rel)
for name,n in EXPECTED.items():
    with (ROOT/name).open(newline='',encoding='utf-8') as f:
        r=csv.DictReader(f)
        headers=r.fieldnames or []
        if not headers or len(headers)!=len(set(headers)) or any(not h for h in headers):
            raise SystemExit('Malformed CSV header: '+name)
        if any(any(x in h.lower().replace(' ','') for x in ('responseid','respondentid','participantid','email','phone','latitude','longitude','pincode','ipaddress')) for h in headers):
            raise SystemExit('Identifier column: '+name)
        rows=list(r)
        if len(rows)!=n or any(None in row or any(v is None for v in row.values()) for row in rows):
            raise SystemExit('Unexpected CSV shape: '+name)
print('PASS: exact allowlist for 14 empirical aggregate CSVs and 1 derived example CSV; basic direct-identifier screen')
print('NOTICE: small-cell/linkage and institutional disclosure decisions remain manual')
