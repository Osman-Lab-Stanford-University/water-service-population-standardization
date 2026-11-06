"""Reference AIPW estimator for five water-service groups, tested on synthetic data.

This module implements reference source-specific mean and distributional
standardization on prepared arrays, with synthetic-data test cases.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence
import warnings
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import log_loss

SOURCE_NAMES = ("Municipal piped", "Groundwater", "Surface water", "Tanker delivery", "Commercial water")

@dataclass
class AggregateResult:
    names: tuple[str,...]
    source_counts: np.ndarray
    observed_means: np.ndarray
    observed_sample_sds: np.ndarray
    standardized_means: np.ndarray
    thresholds: np.ndarray
    standardized_cdf_percent: np.ndarray  # T x 5
    unweighted_max_abs_smd: np.ndarray
    weighted_max_abs_smd: np.ndarray
    weighted_max_category_gap_pp: np.ndarray
    ess: np.ndarray
    max_weight_share: np.ndarray
    first_percentile_propensity: np.ndarray
    propensity_out_of_fold_logloss: float
    constant_indicator_count: int


def fit_reference(X: np.ndarray, source: np.ndarray, score: np.ndarray, *,
                  thresholds: Sequence[float] | None = None, seed: int = 20261022) -> AggregateResult:
    """Independently fit cross-fitted nuisance models and return group-level outputs.

    X must be the **original, privately constructed** numeric one-hot design;
    source is 1,2,3,4,5; score is 0–100. Exactly five source classes expected.
    The caller must authenticate the data loader against the original study.
    No microdata or fitted individual predictions are written to disk.
    """
    X = np.asarray(X, dtype=float)
    raw_source=np.asarray(source)
    if raw_source.ndim!=1 or not np.issubdtype(raw_source.dtype, np.number):
        raise ValueError('Source codes must be a numeric one-dimensional array')
    if not np.all(np.isfinite(raw_source)) or not np.all(raw_source==np.round(raw_source)):
        raise ValueError('Source labels must be finite integer codes (1 through 5)')
    source = raw_source.astype(int)
    score = np.asarray(score, dtype=float)
    n = len(score)
    if X.ndim!=2 or len(X)!=n or len(source)!=n or n<100:raise ValueError('Invalid input dimensions')
    if np.any(~np.isfinite(X)) or np.any(~np.isfinite(score)) or not np.all((score>=0)&(score<=100)):
        raise ValueError('Non-finite features or WSQS outside [0,100]')
    if not np.all((X==0)|(X==1)):
        raise ValueError('The input X must be the study-style binary indicator matrix')
    if not np.array_equal(np.unique(source),np.arange(1,6)):raise ValueError('Expected five source classes labeled 1..5')
    counts=np.array([np.sum(source==j+1) for j in range(5)],dtype=int)
    if np.any(counts<10):raise ValueError('A source group has fewer than 10 records')
    if thresholds is None: thresholds=np.linspace(0,100,41)
    ts=np.asarray(thresholds,dtype=float)
    if ts.ndim!=1 or len(ts)==0 or np.any(~np.isfinite(ts)) or np.any(np.diff(ts)<=0) or np.any((ts<0)|(ts>100)):
        raise ValueError('Supply nonempty, strictly increasing thresholds on the WSQS 0–100 scale')
    k=len(ts)
    e=np.full((n,5),np.nan,dtype=float)
    m=np.full((n,5),np.nan,dtype=float)
    q=np.full((n,5,k),np.nan,dtype=float)
    fold=StratifiedKFold(n_splits=5,shuffle=True,random_state=seed)
    for train,test in fold.split(X,source):
        p=LogisticRegression(C=1.0,solver='lbfgs',max_iter=2000,
                             tol=1e-6,fit_intercept=True)
        with warnings.catch_warnings():
            warnings.simplefilter('error',ConvergenceWarning)
            p.fit(X[train],source[train])
        if not np.array_equal(p.classes_,np.arange(1,6)):raise ValueError('Propensity class ordering changed')
        e[test]=p.predict_proba(X[test])
        for j in range(5):
            selected=train[source[train]==j+1]
            reg=HistGradientBoostingRegressor(loss='squared_error',learning_rate=0.05,
                max_iter=200,max_leaf_nodes=7,max_depth=3,min_samples_leaf=10,
                l2_regularization=1.0,early_stopping=False,random_state=seed)
            reg.fit(X[selected],score[selected])
            m[test,j]=reg.predict(X[test])
            for it,t in enumerate(ts):
                binary=(score[selected]<=t).astype(int)
                if np.all(binary==binary[0]):
                    q[test,j,it]=float(binary[0])
                else:
                    classifier=LogisticRegression(C=1.0,solver='lbfgs',
                        max_iter=3000,tol=1e-6,fit_intercept=True)
                    with warnings.catch_warnings():
                        warnings.simplefilter('error',ConvergenceWarning)
                        classifier.fit(X[selected],binary)
                    q[test,j,it]=classifier.predict_proba(X[test])[:,1]
    if not np.isfinite(e).all() or np.any(e<=0) or not np.isfinite(m).all() or not np.isfinite(q).all():
        raise ValueError('Invalid out-of-fold nuisance predictions')
    if not np.allclose(e.sum(axis=1),1,atol=1e-9):raise ValueError('Out-of-fold probabilities do not sum to 1')

    ref=X.mean(axis=0)
    sd=np.sqrt(ref*(1-ref))
    use=sd>1e-12
    obs=np.empty(5);samplesd=np.empty(5);std=np.empty(5)
    base_max=np.empty(5);after_max=np.empty(5);max_gap=np.empty(5)
    ess=np.empty(5);wt_share=np.empty(5);p01=np.empty(5)
    cdf=np.empty((k,5))
    for j in range(5):
        flag=source==j+1
        w=1/e[flag,j]
        obs[j]=np.mean(score[flag]);samplesd[j]=np.std(score[flag],ddof=1)
        std[j]=np.mean(m[:,j]+flag*(score-m[:,j])/e[:,j])
        cdf[:,j]=[100*np.mean(q[:,j,it]+flag*((score<=t).astype(float)-q[:,j,it])/e[:,j]) for it,t in enumerate(ts)]
        src=X[flag].mean(axis=0)
        weighted=np.average(X[flag],axis=0,weights=w)
        base_max[j]=np.max(np.abs((src[use]-ref[use])/sd[use])) if np.any(use) else np.nan
        after_max[j]=np.max(np.abs((weighted[use]-ref[use])/sd[use])) if np.any(use) else np.nan
        max_gap[j]=100*np.max(np.abs(weighted-ref))
        ess[j]=np.sum(w)**2/np.sum(w*w)
        wt_share[j]=np.max(w)/np.sum(w)
        p01[j]=np.quantile(e[:,j],.01)
    # Distinct threshold models do not mathematically guarantee that the raw
    # AIPW curve is monotone or within [0,100]. Diagnose without silently
    # clipping or applying isotonic adjustments.
    if np.any((cdf < -1e-8) | (cdf > 100+1e-8)) or np.any(np.diff(cdf,axis=0) < -1e-8):
        warnings.warn('Estimated AIPW thresholds are not all valid CDF points; inspect diagnostics before reporting.',RuntimeWarning)
    return AggregateResult(SOURCE_NAMES,counts,obs,samplesd,std,ts,cdf,
        base_max,after_max,max_gap,ess,wt_share,p01,
        float(log_loss(source,e,labels=np.arange(1,6))),int(np.sum(~use)))
