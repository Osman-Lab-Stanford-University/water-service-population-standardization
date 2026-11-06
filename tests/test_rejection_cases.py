#!/usr/bin/env python3
"""Negative controls: make deliberate copies of public files to ensure checkers reject defects.

All mutations run against temporary copies of repository files.
"""
from pathlib import Path
import shutil, subprocess, sys, tempfile, os

ROOT = Path(__file__).resolve().parents[1]
CASES = [
    ('bad_group_count', 'data/aggregates/main_table_1_observed_standardized_means.csv', '1,275', '1,274', 'tests/test_data_contract.py'),
    ('wrong_CDF_threshold', 'data/aggregates/full_41_threshold_standardized_cdf_percent.csv', '2.500,', '2.600,', 'tests/test_data_contract.py'),
    ('wrong_CDF_endpoint', 'data/aggregates/full_41_threshold_standardized_cdf_percent.csv', '100,100.000,100.000,100.000,100.000,100.000,0', '100,99.000,100.000,100.000,100.000,100.000,0', 'tests/test_data_contract.py'),
    ('wrong_domestic_source', 'data/aggregates/si_table_s3_drinking_domestic_source_comparison.csv', '1,335', '1,334', 'tests/test_data_contract.py'),
    ('wrong_covariate_N', 'data/metadata/six_characteristic_design_schema.json', '3066', '3067', 'tests/test_data_contract.py'),
    ('wrong_original_QID', 'data/metadata/six_characteristic_design_schema.json', 'Q15.17', 'Q15.99', 'scripts/verify_source_design_contract.py'),
    ('fake_email', 'data/aggregates/si_table_s2_questionnaire_item_mapping.csv', 'Q1.1', 'fake'+'@'+'example.edu', 'scripts/release_screen.py'),
    ('wrong_S5_interval', 'data/aggregates/si_table_s5_four_contrast_working_intervals.csv', '[−1.894, 0.691]', '[−1.990, 0.691]', 'scripts/validate_public_package.py'),
    ('unreviewed_csv', 'examples/reconstructed_figure_s3_pairwise_distances.csv', 'distance', 'distance_NEW', 'scripts/verify_file_manifest.py'),
    ('checksum_tamper', 'docs/PROVENANCE.md', '# Numerical provenance', '# EDITED Numerical provenance', 'scripts/verify_file_manifest.py'),
]

def run_case(name, relative, old, new, checker):
    with tempfile.TemporaryDirectory(prefix='wrr_negative_') as tmp:
        path=Path(tmp)/'repository'
        shutil.copytree(ROOT,path,ignore=shutil.ignore_patterns('__pycache__','*.pyc','generated_figures','generated_tables'))
        target=path/relative
        txt=target.read_text(encoding='utf-8')
        if old not in txt: raise AssertionError(f'{name}: mutation target missing: {old}')
        target.write_text(txt.replace(old,new,1),encoding='utf-8')
        env=os.environ.copy(); env['PYTHONDONTWRITEBYTECODE']='1'
        proc=subprocess.run([sys.executable,str(path/checker)],cwd=path,env=env,capture_output=True,text=True)
        if proc.returncode==0: raise AssertionError(f'{name}: BAD RECORD WAS ACCEPTED')
        print('REJECTED:',name)

if __name__=='__main__':
    for case in CASES: run_case(*case)
    with tempfile.TemporaryDirectory(prefix='wrr_extra_csv_') as tmp:
        candidate=Path(tmp)/'repository'
        shutil.copytree(ROOT,candidate,ignore=shutil.ignore_patterns('__pycache__','*.pyc','generated_figures','generated_tables'))
        (candidate/'examples/unreviewed_extra.csv').write_text('x\n1\n',encoding='utf-8')
        p=subprocess.run([sys.executable,str(candidate/'scripts/release_screen.py')],cwd=candidate,env=os.environ,capture_output=True,text=True)
        if p.returncode==0: raise AssertionError('Unapproved CSV outside data/ was accepted')
        print('REJECTED: extra_CSV_outside_allowlist')
    print(f'PASS: {len(CASES)+1} deliberate corruption cases rejected')
