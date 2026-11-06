"""Synthetic smoke test for the educational reference implementation; no study data."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from aipw_reference import fit_reference

rng=np.random.default_rng(3200)
n=150
X=rng.binomial(1,.4,size=(n,9)).astype(float)
source=np.repeat(np.arange(1,6),30)
score=np.full(n,50.0)
result=fit_reference(X,source,score,thresholds=[0,50,100])
assert result.source_counts.tolist()==[30]*5
assert np.allclose(result.observed_means,50,atol=1e-9)
assert np.allclose(result.standardized_means,50,atol=1e-8)
assert np.allclose(result.standardized_cdf_percent[:,0],[0,100,100],atol=1e-8)
assert np.allclose(result.standardized_cdf_percent,np.repeat([[0],[100],[100]],5,axis=1),atol=1e-8)
assert np.all(result.ess>0) and np.all(result.ess<=30+1e-8)
print('PASS — source-stratified fivefold AIPW reference estimator; constant outcome/CDF event fallbacks')
