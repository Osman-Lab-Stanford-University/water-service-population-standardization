#!/usr/bin/env python3
"""Synthetic-only tests for the archived outcome-learner diagnostic code."""
import importlib.util
import sys
import tempfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'code/original/E28'
sys.path.insert(0,str(SRC))
spec=importlib.util.spec_from_file_location('wrr_e28',SRC/'run_e28_diagnostics.py')
d=importlib.util.module_from_spec(spec)
spec.loader.exec_module(d)
rng=np.random.default_rng(2800)
g=80;n=5*g
X=rng.normal(size=(n,33))
y=np.repeat(np.arange(1,6),g)
s=np.clip(72+2*X[:,0]+np.repeat([0,1,2,3,4],g)+rng.normal(0,9,n),0,100)
e=np.full((n,5),.2)
base=np.tile(np.clip(72+2*X[:,0],0,100)[:,None],(1,5))
preds={'HGBR':base,'RF':base.copy(),'RIDGE':base.copy(),'HGBR_TUNED':base.copy()}
paired=d.boot_rmse_pairs(y,s,preds,bootstrap=220)
assert len(paired)==18
assert all(abs(row['RMSE_difference'])<1e-12 for row in paired)
class Engine:
    @staticmethod
    def aipw_scores(z,score,proba,pr):
        return pr+(np.asarray(z)[:,None]==np.asarray(d.LABELS)[None,:])*(score[:,None]-pr)/proba
means,contrasts,changes=d.aipw_working_summary(y,s,e,preds,Engine)
assert (len(means),len(contrasts),len(changes))==(20,40,45)
assert abs(next(c['contrast_AIPW'] for c in contrasts if c['learner']=='HGBR' and c['contrast']=='PIPED-minus-GROUNDWATER')-
       (next(v['mean_AIPW'] for v in means if v['learner']=='HGBR' and v['source']=='PIPED')-
        next(v['mean_AIPW'] for v in means if v['learner']=='HGBR' and v['source']=='GROUNDWATER')))<1e-10
# Ensure the error-difference bootstrap reacts to genuine predictive improvement.
perfect=dict(preds);perfect['RF']=np.broadcast_to(s[:,None],(n,5)).copy()
p=d.boot_rmse_pairs(y,s,perfect,bootstrap=220)
assert next(x for x in p if x['source']=='POOLED' and x['comparison']=='RF minus HGBR')['paired_bootstrap_upper']<0
# This is an intentionally smaller two-candidate grid, synthetic data only.
cfg=[dict(d.GRID[0],max_iter=12),dict(d.GRID[1],max_iter=15)]
tuned,picks,inner=d.do_tuned_hgbr(X,y,s,verbose=False,cfgs=cfg)
assert tuned.shape==(n,5) and np.isfinite(tuned).all()
assert all(sum(picks[l].values())==5 for l in d.LABELS)
print('PASS: E28 synthetic-only paired bootstrap, AIPW score contrasts, influence-score comparison and nested HGBR tuning')
