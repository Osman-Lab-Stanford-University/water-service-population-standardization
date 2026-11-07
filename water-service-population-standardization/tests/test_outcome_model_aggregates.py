#!/usr/bin/env python3
"""Fail-closed E27–E28 aggregate consistency checks (no private data)."""
import csv, json, math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
AGG=ROOT/'data/aggregates'
NAMES=['Municipal piped','Groundwater','Surface water','Tanker delivery','Commercial water']
LABELS=['PIPED','GROUNDWATER','SURFACE','TANKER','COMMERCIAL']
COUNTS=[1275,1212,234,120,225]
METHODS=['HGBR','RF','RIDGE','HGBR_TUNED']
COLS=['HGBR_original','Random_Forest','Ridge','CV_selected_HGBR']

def read(name,rows):
    with (AGG/name).open(newline='',encoding='utf-8') as f:
        r=csv.DictReader(f); out=list(r)
    assert len(out)==rows and len(set(r.fieldnames))==len(r.fieldnames),name
    assert all(None not in row and all(v is not None for v in row.values()) for row in out),name
    return out

def near(x,y,eps=.0011):
    assert abs(float(x)-float(y))<=eps,(x,y,eps)

def ordered(data,key,names):
    assert [x[key] for x in data]==names,([x[key] for x in data],names)

s8=read('si_table_s8_outcome_model_sensitivity.csv',8)
ordered(s8,'Metric',NAMES+['Municipal piped minus groundwater','Pooled out-of-fold RMSE','Piped-minus-groundwater working SE'])
assert [int(x['n']) for x in s8[:5]]==COUNTS
assert all(x['n']=='' for x in s8[5:])
primary=read('main_table_1_observed_standardized_means.csv',5)
for i in range(5):
    near(s8[i]['HGBR_original'],primary[i]['Standardized mean WSQS'],.0011)
    for c in COLS:assert 0<=float(s8[i][c])<=100
for c in COLS:
    near(float(s8[0][c])-float(s8[1][c]),s8[5][c],.0021)
for c in COLS:
    order=sorted(range(5),key=lambda i:float(s8[i][c]),reverse=True)
    assert order==[2,1,0,4,3],(c,order)
rmse=read('outcome_model_oof_rmse.csv',6)
ordered(rmse,'Source',NAMES+['Pooled'])
assert [int(r['n']) for r in rmse]==COUNTS+[3066]
R_COLS=['HGBR_original_RMSE','RF_RMSE','Ridge_RMSE','CV_selected_HGBR_RMSE']
for j,c in enumerate(R_COLS):
    for r in rmse:assert float(r[c])>0
    pooled=math.sqrt(sum(COUNTS[i]*float(rmse[i][c])**2 for i in range(5))/3066)
    near(pooled,rmse[5][c],.00013)
    near(rmse[5][c],s8[6][COLS[j]],.001)
paired=read('outcome_model_paired_rmse_bootstrap.csv',18)
for i,src in enumerate(LABELS+['POOLED']):
    for j,m in enumerate(['RF','RIDGE','HGBR_TUNED']):
        r=paired[3*i+j]
        assert r['Source']==src and r['Alternative_learner']==m and int(r['Bootstrap_replicates'])==2000
        assert int(r['n'])==(COUNTS[i] if i<5 else 3066)
        low,high,delta=map(float,[r['Bootstrap_lower_2p5'],r['Bootstrap_upper_97p5'],r['RMSE_difference_vs_HGBR']])
        assert low<=high and low<=delta<=high
        near(delta,float(rmse[i][R_COLS[j+1]])-float(rmse[i][R_COLS[0]]),.00013)
tune=read('outcome_model_hgbr_tuning_counts.csv',5)
ordered(tune,'Source',LABELS)
for row in tune:
    assert int(row['Outer_fold_count'])==5
    assert [int(row[c]) for c in ['Original','Conservative','Flexible','Slow_growth']]==[0,5,0,0]
means=read('outcome_model_working_if_means.csv',20)
for i,src in enumerate(LABELS):
    for j,m in enumerate(METHODS):
        r=means[4*i+j]
        assert r['Source']==src and r['Learner']==m and int(r['n'])==COUNTS[i]
        assert float(r['Working_IF_SE'])>0
        near(r['AIPW_standardized_mean'],s8[i][COLS[j]],.00051)
contrasts=read('outcome_model_working_if_piped_groundwater.csv',4)
for j,m in enumerate(METHODS):
    r=contrasts[j]
    assert r['Learner']==m and float(r['Working_IF_SE'])>0
    near(r['AIPW_contrast'],s8[5][COLS[j]],.00051)
    near(r['Working_IF_SE'],s8[7][COLS[j]],.00051)
    near(r['AIPW_contrast'],float(means[j]['AIPW_standardized_mean'])-float(means[4+j]['AIPW_standardized_mean']),.00001)
    near(r['Difference_from_original_HGBR'],float(r['AIPW_contrast'])-float(contrasts[0]['AIPW_contrast']),.00001)
    if j==0:assert r['Working_IF_SE_of_difference']==''
    else: assert 0<float(r['Working_IF_SE_of_difference'])<1
meta=json.loads((ROOT/'data/metadata/outcome_model_sensitivity_specification.json').read_text())
assert meta['sample_size']==3066 and list(meta['group_counts'].values())==COUNTS
assert meta['outer_seed']==20261022 and meta['outer_cross_folds']==5 and meta['encoded_indicator_count']==33
assert len(meta['covariate_questions'])==6 and len(meta['hgbr_tuning']['grid'])==4
assert meta['hgbr_tuning']['selected_in_25_of_25']=='conservative'
assert meta['paired_rmse_bootstrap']['resamples']==2000
print('PASS: 6 E27–E28 aggregate tables, manuscript Table S8, source contrasts, pooled RMSE, bootstrap, tuning, working AIPW SE and 4-model parameter schema')
