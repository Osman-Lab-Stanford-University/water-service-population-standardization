#!/usr/bin/env python3
"""Run a TINY, isolated no-survey-data E20 smoke test; not original coverage replication.
Requires NumPy/SciPy; intentionally excluded from dependency-free public checks.
"""
from pathlib import Path
import json,os,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/'code/original/E20/simulate_holdout.py'
with tempfile.TemporaryDirectory(prefix='wrr_e20_synthetic_smoke_') as tmp:
    env=os.environ.copy()
    env.update({'E20_NULL_REPS':'2','E20_ALT_REPS':'2','E20_MULTIPLIER_DRAWS':'100','PYTHONDONTWRITEBYTECODE':'1'})
    result=subprocess.run([sys.executable,'-B',str(SCRIPT),str(Path(tmp)/'out')],
                          cwd=tmp,env=env,text=True,capture_output=True,timeout=150)
    if result.returncode:
        raise AssertionError('E20 isolated smoke failed: '+result.stderr[-1500:])
    out=Path(tmp)/'out'
    assert {p.name for p in out.iterdir()}=={'EMPIRICAL_20_REPORT.json'}, 'Unexpected E20 files'
    data=json.loads((out/'EMPIRICAL_20_REPORT.json').read_text(encoding='utf-8'))
    assert data['phase']=='EMPIRICAL_20_INDEPENDENT_ORACLE_FULL_GRID_HOLDOUT'
    assert data['restricted_household_data_read'] is False
    assert data['case_level_data_saved'] is False
    assert data['oracle_nuisance_only'] is True
    assert data['publication_inference_certified'] is False
    assert data['n_per_replication']==3066
    assert data['multiplier_draws']==100
    assert {k:v['replications'] for k,v in data['scenarios'].items()}=={'null':2,'alternative':2}
    assert set(data['three_families_predeclared'])=={'FULL_INTERIOR_39','CENTRAL_20_TO_80','REGISTERED_50_60_70'}
    for case in data['scenarios'].values():
        assert set(case['pairs'])=={'PIPED_GROUNDWATER','PIPED_TANKER'}
        for pair in case['pairs'].values():
            assert set(pair['families'])==set(data['three_families_predeclared'])
print('PASS: reduced E20 original-code simulation smoke; 2+2 synthetic replications, no respondent data')
print('NOTE: this is not reproduction or validation of the manuscript full-size E20 simulation coverage')
