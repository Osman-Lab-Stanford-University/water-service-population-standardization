import hashlib
import io
import json
import runpy
import sys
import warnings
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold

(
    original_script, g7_script, g7_report, workbook,
    scoring_reference, candidate, out_folder,
    original_report, support_report
) = sys.argv[1:]

out = Path(out_folder)
e04_original = json.loads(Path(original_report).read_text())
formal_support = json.loads(Path(support_report).read_text())
if not formal_support['G5_5_all_regimes_pass']:
    raise RuntimeError('Frozen G5.5 eight-covariate support did not pass')
if not formal_support['G5_5_frozen_policy_used']:
    raise RuntimeError('Frozen support implementation not verified')
if formal_support['analytical_n'] != 3066:
    raise RuntimeError('Formal support reference cohort mismatch')

saved_argv = sys.argv[:]
try:
    sys.argv = [original_script, g7_script, g7_report,
                workbook, scoring_reference, candidate, out_folder]
    with redirect_stdout(io.StringIO()):
        env = runpy.run_path(original_script, run_name='__main__')
finally:
    sys.argv = saved_argv

current_e04 = json.loads((out/'EMPIRICAL_04_REPORT.json').read_text())
if current_e04['analytical_n'] != 3066:
    raise RuntimeError('Reconstruction sample mismatch')
labels = ['PIPED','GROUNDWATER','SURFACE','TANKER','COMMERCIAL']
for label in labels:
    if not np.isclose(current_e04['by_regime'][label]['standardized_AIPW_mean'],
                      e04_original['by_regime'][label]['standardized_AIPW_mean'],
                      atol=1e-10, rtol=0):
        raise RuntimeError('E04 standardized-mean parity mismatch')

X = np.asarray(env['X'])
y = np.asarray(env['y'], dtype=int)
z = np.asarray(env['z'], dtype=str)
scores = np.asarray(env['scores'], dtype=float)
e = np.asarray(env['e'], dtype=float)
seed = int(env['seed'])
n = 3066
k = 5
if (X.shape[0] != n or y.shape != (n,) or scores.shape != (n,)
        or e.shape != (n,k) or seed != 20261022
        or list(env['labels']) != labels or not np.isfinite(X).all()
        or not np.isfinite(scores).all() or not np.isfinite(e).all()
        or np.min(scores)<0 or np.max(scores)>100
        or np.min(e)<=0 or np.max(e)>=1):
    raise RuntimeError('Primary model input or design mismatch')
if not np.allclose(e.sum(axis=1),1,atol=1e-10,rtol=0):
    raise RuntimeError('E04 propensity normalization mismatch')
mask = y[:,None] == np.arange(1,k+1)[None,:]
if not np.array_equal(mask, z[:,None] == np.asarray(labels)[None,:]):
    raise RuntimeError('Five-regime alignment mismatch')
folds = list(StratifiedKFold(n_splits=5,shuffle=True,random_state=seed).split(X,y))
if len(folds) != 5:
    raise RuntimeError('Fold count mismatch')
thresholds = np.linspace(0,100,41)
if len(thresholds)!=41 or any(t not in thresholds for t in [50,60,70]):
    raise RuntimeError('Threshold grid mismatch')

raw = np.full((len(thresholds),k),np.nan)
working_contrast_se = np.full(len(thresholds),np.nan)
oneclass = {label:0 for label in labels}

print('EMPIRICAL_15_FITTING_BEGIN',flush=True)
for ti,t in enumerate(thresholds):
    event = (scores <= t).astype(float)
    q = np.full((n,k),np.nan)
    for train,test in folds:
        for j,label in enumerate(labels):
            local_train = train[y[train] == j+1]
            target = event[local_train]
            classes = np.unique(target)
            if len(classes)==1:
                q[test,j] = float(classes[0])
                oneclass[label] += 1
            else:
                estimator = LogisticRegression(
                    penalty='l2', C=1.0, solver='lbfgs',
                    max_iter=3000, tol=1e-6, fit_intercept=True
                )
                with warnings.catch_warnings():
                    warnings.simplefilter('error', ConvergenceWarning)
                    estimator.fit(X[local_train],target)
                if not np.array_equal(estimator.classes_,np.array([0.0,1.0])):
                    raise RuntimeError('Outcome binary class ordering mismatch')
                q[test,j] = estimator.predict_proba(X[test])[:,1]
    if (not np.isfinite(q).all() or np.min(q)<0 or np.max(q)>1):
        raise RuntimeError('Invalid conditional CDF probabilities')
    psi = q + mask*(event[:,None]-q)/e
    if not np.isfinite(psi).all():
        raise RuntimeError('Nonfinite AIPW CDF scores')
    raw[ti] = np.mean(psi,axis=0)
    working_contrast_se[ti] = float(np.std(psi[:,0]-psi[:,1],ddof=1)/np.sqrt(n))
    if ti%10==0 or ti==len(thresholds)-1:
        print('GRID_PROGRESS:',int(ti+1),'OF',len(thresholds),flush=True)

if not np.allclose(raw[-1],1,atol=1e-10,rtol=0):
    raise RuntimeError('CDF upper endpoint must equal one')
# Crosscheck rounded locked E08 logistic-tail figures, tolerance reflects
# previously reported 3-decimal percentage point rounding and minor numerical settings.
anchors = {
    50: [7.469,7.764,6.149,8.690,10.553],
    60: [20.014,19.019,17.819,27.271,17.067],
    70: [38.837,35.072,33.932,45.029,43.679]
}
anchor_checks = {}
max_anchor_delta_pp = 0.0
for t, reference_pp in anchors.items():
    ti = int(round(t/2.5))
    deltas = (100*raw[ti] - np.asarray(reference_pp)).tolist()
    anchor_checks[str(t)] = {labels[j]:float(deltas[j]) for j in range(k)}
    max_anchor_delta_pp = max(max_anchor_delta_pp,max(abs(x) for x in deltas))
print('LOCKED_E08_MAX_ANCHOR_DEVIATION_PP:',round(max_anchor_delta_pp,6),flush=True)
if max_anchor_delta_pp>0.15:
    raise RuntimeError('Locked E08 logistic-tail anchors not reproduced within 0.15 percentage points')

# Raw AIPW thresholdwise estimates are not guaranteed monotone or bounded.
# For display only, a separately labelled isotonic, bounded projection is stored.
projected = np.column_stack([
    IsotonicRegression(y_min=0,y_max=1,increasing=True)
      .fit_transform(thresholds,raw[:,j])
    for j in range(k)
])
if (np.any(np.diff(projected,axis=0)<-1e-12)
    or np.any(projected<0) or np.any(projected>1)):
    raise RuntimeError('Diagnostic isotonic projection invalid')
metrics = {}
for j,label in enumerate(labels):
    deltas = np.diff(raw[:,j])
    metrics[label] = {
        'raw_steps_decreasing': int(np.sum(deltas<-1e-10)),
        'largest_raw_step_decrease_pp':float(max(0,-np.min(deltas))*100),
        'raw_out_of_range_grid_points':int(np.sum((raw[:,j]<0)|(raw[:,j]>1))),
        'maximum_isotonic_adjustment_pp':float(np.max(np.abs(projected[:,j]-raw[:,j]))*100),
        'one_class_training_fits':int(oneclass[label]),
        'area_based_mean_diagnostic':float(100-np.trapezoid(projected[:,j],thresholds))
            if hasattr(np,'trapezoid') else float(100-np.trapz(projected[:,j],thresholds))
    }

contrast = 100*(raw[:,0]-raw[:,1])
max_idx = int(np.argmax(np.abs(contrast)))
rows = []
for i,t in enumerate(thresholds):
    rows.append({
        'threshold_wsqs':float(t),
        'raw_aipw_cdf':{label:float(raw[i,j]) for j,label in enumerate(labels)},
        'display_only_isotonic_cdf':{label:float(projected[i,j]) for j,label in enumerate(labels)},
        'raw_piped_minus_groundwater_pp':float(contrast[i]),
        'working_independent_household_SE_pp':float(100*working_contrast_se[i])
    })
report = {
    'phase':'EMPIRICAL_15_PRIMARY_STANDARDIZED_CDF',
    'status':'COMPUTED_DISTRIBUTIONAL_INFERENCE_NOT_CERTIFIED',
    'analytical_n':n,
    'regimes':labels,
    'grid':{'start':0,'end':100,'increment':2.5,'n_thresholds':41},
    'crossfit':{'folds':5,'seed':seed,'propensity':'frozen E04 cross-fitted logistic',
                'threshold_outcome':'L2 logistic C=1, independently fitted at each threshold'},
    'standardization_population':'empirical distribution of six baseline covariates for all 3066 households',
    'e04_reconstruction_verified':True,
    'g5_5_formal_eight_covariate_support_reference_pass':True,
    'anchor_benchmark':'previously locked E08 rounded percentages',
    'max_anchor_difference_pp':float(max_anchor_delta_pp),
    'anchor_differences_pp':anchor_checks,
    'threshold_rows':rows,
    'cdf_diagnostics':metrics,
    'max_abs_piped_minus_groundwater_difference_pp':float(contrast[max_idx]),
    'max_abs_difference_threshold':float(thresholds[max_idx]),
    'notes':[
        'Raw AIPW CDF points retained even if outside bounds or non-monotone.',
        'Isotonic series is a separately labelled display-only projection, not an alternative fitted model.',
        'Working standard errors exclude recruitment-cluster dependence and are not publication-certified.',
        'No first-order stochastic dominance or distributional significance claim has been tested.',
        'Restricted individual outputs are never written.'
    ],
    'simulation_coverage_certified':False,
    'simultaneous_cdf_inference_certified':False,
    'source_records_exported':False,
    'predictions_exported':False,
    'public_release_authorized':False
}
(out/'EMPIRICAL_15_REPORT.json').write_text(
    json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8'
)
print('EMPIRICAL_15_RESULTS')
print('ANALYTICAL_N:',n)
print('GRID: 0_TO_100_STEP_2_5',len(thresholds),'THRESHOLDS')
print('E04_RECONSTRUCTION: VERIFIED')
print('REGISTERED_E08_ANCHOR_MAX_DELTA_PP:',round(max_anchor_delta_pp,6))
for t in [50,60,70]:
    i = int(t/2.5)
    print('RISK_LE_'+str(t), ' '.join(label+'='+str(round(100*raw[i,j],3))+'%' for j,label in enumerate(labels)))
    print('PIPED_MINUS_GROUNDWATER_PP:',round(contrast[i],4))
for label in labels:
    m = metrics[label]
    print('CDF_DIAGNOSTIC',label,'RAW_DECREASES',m['raw_steps_decreasing'],
          'RAW_OUT_OF_BOUNDS',m['raw_out_of_range_grid_points'],
          'MAX_ISOTONIC_CHANGE_PP',round(m['maximum_isotonic_adjustment_pp'],4))
print('MAX_ABS_PIPED_GROUND_CDF_DIFF_PP:',round(contrast[max_idx],3),
      'AT_WSQS:',float(thresholds[max_idx]))
print('SIMULTANEOUS_INFERENCE_CERTIFIED: NO')
print('HOUSEHOLD_LEVEL_OUTPUTS: NONE')
print('PUBLIC_RELEASE_AUTHORIZED: NO')
print('EMPIRICAL_15_RESULT: COMPLETE')
