"""Nontrivial synthetic regression/propensity test; no empirical study records."""
from pathlib import Path
import numpy as np,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from aipw_reference import fit_reference
rng=np.random.default_rng(4021);n=2000
X=rng.binomial(1,[.15,.35,.55,.45,.20,.70],size=(n,6)).astype(float)
# Source propensity depends nontrivially on covariate indicators.
logits=np.column_stack([.1+1.0*X[:,0]-.4*X[:,1],
 .2-.5*X[:,0]+.3*X[:,2],
 -.5+.8*X[:,1]-.2*X[:,3],
 -.7+.6*X[:,0]+.6*X[:,4],
 -.2-.3*X[:,2]+.9*X[:,5]])
prob=np.exp(logits-logits.max(axis=1,keepdims=True));prob/=prob.sum(axis=1,keepdims=True)
source=1+(rng.random(n)[:,None] > np.cumsum(prob,axis=1)).sum(axis=1)
source=np.minimum(source,5)
base=39.+5.*(source-1)
Y=base+7*X[:,0]+5*X[:,1]+2*X[:,2]
result=fit_reference(X,source,Y,thresholds=[0,45,55,75,100])
target=39+5*np.arange(5)+7*X[:,0].mean()+5*X[:,1].mean()+2*X[:,2].mean()
print('group sizes:',result.source_counts.tolist())
print('known standardized target:',np.round(target,3).tolist())
print('fitted standardized means:',np.round(result.standardized_means,3).tolist())
print('max deviation:',max(abs(result.standardized_means-target)))
print('cdf endpoints:',result.standardized_cdf_percent[0].tolist(),result.standardized_cdf_percent[-1].tolist())
if not max(abs(result.standardized_means-target))<1.5: raise AssertionError('Known-covariate mean target mismatch')
if not np.max(abs(result.standardized_cdf_percent[0]))<1e-6: raise AssertionError('CDF at WSQS 0 mismatch')
if not np.max(abs(result.standardized_cdf_percent[-1]-100))<1e-6: raise AssertionError('CDF at WSQS 100 mismatch')
print('PASS — nontrivial covariate-dependent synthetic target; no study data')
