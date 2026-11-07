#!/usr/bin/env python3
"""Six-characteristic RF and Ridge outcome-model sensitivities for the WRR study.

The empirical cohort is reconstructed in memory from the unmodified E04 program,
using its original decoder/scorer and private inputs. No microdata are written.
Requires access to the original private E04 inputs on the authorized research computer.
"""
import argparse
import ast
import csv
import hashlib
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import sklearn
from sklearn.ensemble import HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import log_loss
from sklearn.model_selection import StratifiedKFold

E04_SHA = '439fc907ff72fdc93071bb61e4ea37cd0ecff68b2948a2fd6942208bcf9d9bf1'
ENGINE_SHA = '10fafe121b64b9ff53bcb992c076d45a7b1097fd0ba76e7b1702dbb2e0cd4645'
SEED = 20261022
LABELS = ['PIPED', 'GROUNDWATER', 'SURFACE', 'TANKER', 'COMMERCIAL']
QIDS = ['Q15.17', 'Q15.13', 'Q15.15', 'Q15.14', 'Q15.10', 'Q15.5']
REFERENCE_MEANS = np.array([74.84673146, 75.44853151, 76.36525295, 71.52978731, 72.29863262])
REFERENCE_N = [1275, 1212, 234, 120, 225]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_original_e04_cohort(e04_code, source, prior, workbook, reference_root, package_root, e04_output):
    """Execute only the original E04 cohort/scoring section (before any model fits)."""
    text = e04_code.read_text(encoding='utf-8')
    tree = ast.parse(text, filename=str(e04_code))
    markers = [i for i, node in enumerate(tree.body)
               if isinstance(node, ast.Assign) and any(
                   isinstance(target, ast.Name) and target.id == 'SOURCE_WORKBOOK'
                   for target in node.targets)]
    if len(markers) != 1:
        raise RuntimeError('E04 cohort-to-model section marker changed')
    prefix = ast.Module(body=tree.body[:markers[0]], type_ignores=[])
    ast.fix_missing_locations(prefix)
    args_original = sys.argv[:]
    ns = {'__name__': '__e04_cohort_loader__'}
    try:
        # The six arguments are the original E04 CLI contract.
        sys.argv = [str(e04_code), str(source), str(prior), str(workbook),
                    str(reference_root), str(package_root), str(e04_output)]
        exec(compile(prefix, str(e04_code), 'exec'), ns)
    finally:
        sys.argv = args_original
    return ns


def fit_nuisance_models(X, y, scores, engine, seed=SEED, folds=5, trees=500, ridge_alpha=10.0, return_oof=False):
    """Keep folds/propensities fixed; compare HGBR, RF, and Ridge outcome learners."""
    n, k = len(y), len(LABELS)
    assert len(X) == len(scores) == n
    e = np.full((n, k), np.nan)
    mh = np.full((n, k), np.nan)
    mr = np.full((n, k), np.nan)
    ml = np.full((n, k), np.nan)
    fold_ids = np.full(n, -1, dtype=int)
    for fold, (train, test) in enumerate(StratifiedKFold(
            n_splits=folds, shuffle=True, random_state=seed).split(X, y)):
        if np.intersect1d(train, test).size:
            raise RuntimeError('Cross-fitting leakage')
        fold_ids[test] = fold
        prop = LogisticRegression(penalty='l2', C=1.0, solver='lbfgs',
                                  max_iter=2000, tol=1e-6, fit_intercept=True)
        with warnings.catch_warnings():
            warnings.simplefilter('error', ConvergenceWarning)
            prop.fit(X[train], y[train])
        if not np.array_equal(prop.classes_, np.arange(1, k + 1)):
            raise RuntimeError('Propensity classes changed')
        e[test, :] = prop.predict_proba(X[test])
        for j in range(k):
            z_train = train[y[train] == (j + 1)]
            if len(z_train) < 5:
                raise RuntimeError('Insufficient outcome training data')
            # EXACT primary E04 histogram gradient-boosting specification.
            hgbr = HistGradientBoostingRegressor(
                loss='squared_error', learning_rate=0.05, max_iter=200,
                max_leaf_nodes=7, max_depth=3, min_samples_leaf=10,
                l2_regularization=1.0, early_stopping=False, random_state=seed)
            # Prespecified alternative outcome learner, trained on SAME fold and data.
            rf = RandomForestRegressor(
                n_estimators=trees, criterion='squared_error',
                bootstrap=True, max_features=1.0, max_depth=None,
                min_samples_leaf=10, random_state=seed, n_jobs=1)
            # A regularized linear benchmark probes dependence on nonlinear outcome learners.
            # The fixed alpha is prespecified (no tuning on held-out outcome data).
            ridge = Ridge(alpha=ridge_alpha, fit_intercept=True, solver='svd')
            hgbr.fit(X[z_train], scores[z_train])
            rf.fit(X[z_train], scores[z_train])
            ridge.fit(X[z_train], scores[z_train])
            mh[test, j] = hgbr.predict(X[test])
            mr[test, j] = rf.predict(X[test])
            # WSQS is a bounded score. Clip only the Ridge predictions to its valid scale.
            ml[test, j] = np.clip(ridge.predict(X[test]), 0.0, 100.0)
    if not np.all(fold_ids >= 0):
        raise RuntimeError('Unassigned fold')
    if not np.isfinite(e).all() or np.min(e) <= 0 or not np.allclose(e.sum(axis=1), 1, atol=1e-9):
        raise RuntimeError('Propensity invalid')
    if not np.isfinite(mh).all() or not np.isfinite(mr).all() or not np.isfinite(ml).all():
        raise RuntimeError('Outcome predictions nonfinite')
    if np.min(mh) < -1e-8 or np.max(mh) > 100 + 1e-8 or np.min(mr) < -1e-8 or np.max(mr) > 100 + 1e-8:
        raise RuntimeError('Outcome prediction outside WSQS scale')
    z = np.asarray(LABELS, dtype=str)[np.asarray(y, dtype=int)-1]
    predictions = {'HGBR': mh, 'RF': mr, 'RIDGE': ml}
    outputs = {}
    for learner, nuisance in predictions.items():
        uncapped = engine.aipw_scores(z, scores, e, nuisance)
        cap99, _ = engine.capped_scores(z, scores, e, nuisance, 99)
        cap95, _ = engine.capped_scores(z, scores, e, nuisance, 95)
        # Score equation check is independent of engine implementation.
        check = nuisance + (z[:, None] == np.asarray(LABELS)[None,:]) * (scores[:,None]-nuisance) / e
        if not np.allclose(uncapped, check, atol=1e-8, rtol=0):
            raise RuntimeError('AIPW estimating-equation mismatch')
        outputs[learner] = {
            'uncapped': np.mean(uncapped, axis=0),
            'cap99': np.mean(cap99, axis=0),
            'cap95': np.mean(cap95, axis=0),
            'oof_rmse': [float(np.sqrt(np.mean((scores[y==j+1] - nuisance[y==j+1, j])**2)))
                         for j in range(k)]
        }
    result={'e':e,'fold_ids':fold_ids,'outcomes':outputs, 'logloss':float(log_loss(y,e,labels=[1,2,3,4,5]))}
    if return_oof:
        result['predictions_in_memory_only']=predictions
    return result



def evaluate_e02_numerical_agreement(X, y, e, logloss_current, reference02):
    """Independent anchored check against *original* E02 group-level diagnostics.

    The historical E04 code enforced 1e-9 log-loss parity in its original runtime.
    For a separately installed local environment, near-agreement is only accepted
    if both the log loss and each group-specific weighting diagnostic agree with
    registered E02 output at fixed, documented tolerances. No values are changed.
    """
    ref = reference02['model_results']['LOGISTIC_L2_C1']
    expected = float(ref['out_of_fold_logloss'])
    delta = float(logloss_current) - expected
    mode = 'exact' if abs(delta) <= 1e-9 else 'numerically_close'
    limits = {
        'n': 0,
        'ess': 0.05, # effective sample size: 0.05 household equivalent
        'ess_fraction': 0.0005,
        'largest_weight_fraction': 0.00005,
        'p01_target': 0.0002,
        'max_absolute_category_balance_gap': 0.0002,
    }
    comparisons = {}
    issues = []
    print('\nE02 NUMERICAL PROPENSITY VERIFICATION:', flush=True)
    print(f'  original log loss: {expected:.15f}', flush=True)
    print(f'  refitted log loss: {logloss_current:.15f}', flush=True)
    print(f'  signed difference: {delta:+.15e}', flush=True)
    print('  exact threshold: 1e-9 | numerical-near threshold: 1e-6', flush=True)
    if not math.isfinite(delta) or abs(delta)>1e-6:
        issues.append(f'log loss |delta|={abs(delta):.9g} exceeds 1e-6')
    for j,label in enumerate(LABELS):
        observed = np.asarray(y)==(j+1)
        weights = 1.0 / np.asarray(e)[observed,j]
        actual = {
            'n':int(observed.sum()),
            'ess':float(weights.sum()**2 / np.square(weights).sum()),
            'ess_fraction':float(weights.sum()**2 / np.square(weights).sum() / observed.sum()),
            'largest_weight_fraction':float(weights.max() / weights.sum()),
            'p01_target':float(np.quantile(e[:,j],0.01)),
            'max_absolute_category_balance_gap':float(np.max(np.abs(
                np.average(X[observed],axis=0,weights=weights) - X.mean(axis=0)
            ))),
        }
        originals = ref['regimes'][label]
        comparisons[label] = {}
        print(f'  {label}:',flush=True)
        for key,threshold in limits.items():
            old = float(originals[key])
            new = float(actual[key])
            diff = new-old
            ok = math.isfinite(diff) and abs(diff) <= threshold + 1e-12
            if not ok: issues.append(f'{label} {key} |delta|={abs(diff):.9g} exceeds {threshold}')
            comparisons[label][key] = {'E02':old,'refit':new,'difference':diff,'limit':threshold,'pass':ok}
            print(f'    {key:<34} E02={old:.9f} refit={new:.9f} delta={diff:+.9f} '+('PASS' if ok else 'FAIL'), flush=True)
    if issues:
        raise RuntimeError('E02 numerical equivalence FAILED: '+'; '.join(issues))
    if any(abs(v['difference']) > 1e-9 for group in comparisons.values() for v in group.values()):
        mode = 'numerically_close'
    print(f'PASS: E02 {mode} numerical agreement (log loss + all 30 group-level diagnostic fields)',flush=True)
    return {'level':mode,'original_logloss':expected,'refitted_logloss':float(logloss_current),
            'difference':delta,'logloss_acceptance_limit':1e-6,
            'group_comparisons':comparisons,'diagnostic_limits':limits}


def evaluate_e04_aggregate_agreement(output_dir, outcomes):
    """Compare the primary HGBR mean, two capped values and RMSE with saved E04."""
    original=Path(output_dir)/'EMPIRICAL_04_REPORT.json'
    if not original.is_file():
        raise RuntimeError(f'Missing original empirical reference report: {original}')
    report=json.loads(original.read_text())
    mapping={'uncapped':'standardized_AIPW_mean','cap99':'capped99_mean',
             'cap95':'capped95_mean','oof_rmse':'outcome_OOF_RMSE_observed_group'}
    compared={}
    issues=[]
    print('\nE04 SAVED OUTCOME RESULTS VERIFICATION:',flush=True)
    for j,label in enumerate(LABELS):
        compared[label]={}
        for current,old_key in mapping.items():
            original_value=float(report['by_regime'][label][old_key])
            new_value=float(outcomes['HGBR'][current][j])
            d=new_value-original_value
            compared[label][current]={'E04':original_value,'refit':new_value,'difference':d,'limit':0.005}
            if not math.isfinite(d) or abs(d)>0.005:
                issues.append(f'{label} {current} |delta|={abs(d):.9g} >0.005')
        deviations=[abs(v['difference']) for v in compared[label].values()]
        print(f'  {label:<12}: max |delta|={max(deviations):.9f} WSQS (mean, cap99, cap95, RMSE)',flush=True)
    if issues: raise RuntimeError('E04 HGBR result comparison FAILED: '+'; '.join(issues))
    print('PASS: all 20 E04 HGBR aggregate comparisons within 0.005 WSQS points',flush=True)
    return compared

def main():
    ap = argparse.ArgumentParser(description='Run RF and Ridge vs HGBR outcome learners on the final E04 cohort')
    ap.add_argument('--e04-code', type=Path, required=True, help='Unmodified E04/run_aipw.py')
    ap.add_argument('--decoder', type=Path, required=True, help='Original G7 survey decoder source (E04 first CLI argument)')
    ap.add_argument('--g7-report', type=Path, required=True, help='Original G7 JSON report (E04 second CLI argument)')
    ap.add_argument('--workbook', type=Path, required=True, help='Restricted original questionnaire workbook')
    ap.add_argument('--ref-root', type=Path, required=True, help='Original reference root containing scoring and crosswalk files')
    ap.add_argument('--package-root', type=Path, required=True, help='Original water_standardize package root')
    ap.add_argument('--e04-output', type=Path, required=True, help='Original E04 output with _empirical02/_empirical03_reference.json')
    ap.add_argument('--output-dir', type=Path, required=True, help='New PRIVATE output directory outside public repo')
    ap.add_argument('--trees', type=int, default=500)
    ap.add_argument('--ridge-alpha', type=float, default=10.0)
    args = ap.parse_args()
    if args.trees < 100:
        ap.error('At least 100 trees required')
    if not np.isfinite(args.ridge_alpha) or args.ridge_alpha <= 0:
        ap.error('Ridge alpha must be positive and finite')
    if sha256(args.e04_code) != E04_SHA:
        raise RuntimeError('E04 source SHA differs from registered original')
    engine_path = args.package_root / 'src/water_standardize/core.py'
    if sha256(engine_path) != ENGINE_SHA:
        raise RuntimeError('Frozen E04 AIPW engine SHA differs from registered original')
    # A separate output directory prevents changing historical stages or publishing private files.
    output_dir = args.output_dir.resolve()
    if output_dir == args.e04_output.resolve() or args.e04_output.resolve() in output_dir.parents:
        raise RuntimeError('Choose an output directory separate from the original E04 output')
    ns = load_original_e04_cohort(args.e04_code, args.decoder, args.g7_report,
                                  args.workbook, args.ref_root, args.package_root, args.e04_output)
    X, y, scores = (np.asarray(ns[v]) for v in ('X', 'y', 'scores'))
    if X.shape != (3066, 33) or len(y) != 3066 or len(scores) != 3066:
        raise RuntimeError(f'Wrong cohort shape: X={X.shape}, scores={scores.shape}')
    if [int(np.sum(y == z)) for z in range(1, 6)] != REFERENCE_N:
        raise RuntimeError('Regime counts changed')
    if list(ns['QIDS']) != QIDS:
        raise RuntimeError('Adjustment variables changed')
    sys.path.insert(0, str(args.package_root / 'src'))
    from water_standardize import core as engine
    fit = fit_nuisance_models(X, y, scores, engine, trees=args.trees, ridge_alpha=args.ridge_alpha)
    # Prevent version-drift or decoder changes from going unnoticed.
    h = fit['outcomes']['HGBR']['uncapped']
    h_diff = h - REFERENCE_MEANS
    print('\nBASELINE HGBR VERIFICATION (all five source groups):', flush=True)
    for name, got, expected, delta in zip(LABELS, h, REFERENCE_MEANS, h_diff):
        print(f'  {name:<12} reproduced={got:.9f} E04={expected:.9f} difference={delta:+.9f}', flush=True)
    if not np.allclose(h, REFERENCE_MEANS, rtol=0, atol=.005):
        raise RuntimeError('Original HGBR estimates did not reproduce (tolerance 0.005 WSQS points)')
    print('PASS: all five primary HGBR means reproduce within 0.005 WSQS points', flush=True)
    reference02_path = args.e04_output / '_empirical02_reference.json'
    reference02 = json.loads(reference02_path.read_text())
    previous_ll = float(reference02['model_results']['LOGISTIC_L2_C1']['out_of_fold_logloss'])
    logloss_delta = fit['logloss'] - previous_ll
    print('\nE02 PROPENSITY PARITY:', flush=True)
    print(f'  expected E02 OOF log loss : {previous_ll:.15f}', flush=True)
    print(f'  reproduced OOF log loss  : {fit["logloss"]:.15f}', flush=True)
    print(f'  signed difference        : {logloss_delta:+.15e}', flush=True)
    print('  E02 historical exact threshold : 1e-9; near-equivalence is separately audited', flush=True)
    print(f'  scikit-learn version     : {sklearn.__version__}', flush=True)
    print(f'  NumPy version            : {np.__version__}', flush=True)
    # The original E04 execution enforced exact 1e-9 parity; in this new runtime,
    # use a documented near-equivalence rule with six original E02 diagnostics per group.
    agreement = evaluate_e02_numerical_agreement(X,y,fit['e'],fit['logloss'],reference02)
    e04_verification = evaluate_e04_aggregate_agreement(args.e04_output,fit['outcomes'])
    diagnostic_groups = {}
    for j, label in enumerate(LABELS):
        observed = y == (j+1)
        w = 1.0 / fit['e'][observed, j]
        diagnostic_groups[label] = {
            'n': int(np.sum(observed)),
            'ess': float(w.sum()**2 / np.square(w).sum()),
            'ess_fraction': float(w.sum()**2 / np.square(w).sum() / observed.sum()),
            'largest_weight_fraction': float(w.max() / w.sum()),
            'max_absolute_category_balance_gap': float(np.max(np.abs(
                np.average(X[observed], axis=0, weights=w) - X.mean(axis=0)
            ))),
            'p01_propensity_target': float(np.quantile(fit['e'][:,j], 0.01)),
        }
    out = {'status':'completed_verified','model_contrast':'outcome learner only HGBR vs RF and Ridge',
           'n':3066, 'source_group_counts':dict(zip(LABELS, REFERENCE_N)),
           'source_diagnostics':diagnostic_groups,
           'applicable_domain_counts':{str(k):int(v) for k,v in ns['domains'].items()},
           'qids':QIDS, 'seed':SEED,'folds':5,'rf_hyperparameters':{
               'n_estimators':args.trees, 'criterion':'squared_error', 'bootstrap':True,
               'max_features':1.0,'max_depth':None,'min_samples_leaf':10,'random_state':SEED,
               'n_jobs':1},
           'ridge_hyperparameters':{'alpha':args.ridge_alpha, 'fit_intercept':True,
                                    'solver':'svd', 'prediction_clip':[0.0,100.0]},
           'baseline_hgbr_reproduced':True,'propensity_validation_passed':True,
           'propensity_parity_level':agreement['level'],
           'propensity_parity_validation':agreement,
           'e04_aggregate_comparison':e04_verification,
           'e04_source_sha256':E04_SHA,'engine_source_sha256':ENGINE_SHA,
           'sklearn_version':sklearn.__version__,'numpy_version':np.__version__,
           'propensity_oof_logloss':fit['logloss'],
           'group_means':{label: {method:{key: float(v[j]) if key!='oof_rmse' else float(v[j])
                    for key,v in fit['outcomes'][method].items()}
                    for method in ['HGBR','RF','RIDGE']}
                    for j,label in enumerate(LABELS)}}
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for j,label in enumerate(LABELS):
        rh = fit['outcomes']['HGBR']; rr = fit['outcomes']['RF']; rl = fit['outcomes']['RIDGE']
        rows.append({'source':label,'n':REFERENCE_N[j],
                     'HGBR_uncapped':f'{rh["uncapped"][j]:.6f}',
                     'RF_uncapped':f'{rr["uncapped"][j]:.6f}',
                     'RIDGE_uncapped':f'{rl["uncapped"][j]:.6f}',
                     'RF_minus_HGBR':f'{rr["uncapped"][j]-rh["uncapped"][j]:.6f}',
                     'RIDGE_minus_HGBR':f'{rl["uncapped"][j]-rh["uncapped"][j]:.6f}',
                     'HGBR_cap99':f'{rh["cap99"][j]:.6f}',
                     'RF_cap99':f'{rr["cap99"][j]:.6f}',
                     'RIDGE_cap99':f'{rl["cap99"][j]:.6f}',
                     'HGBR_cap95':f'{rh["cap95"][j]:.6f}',
                     'RF_cap95':f'{rr["cap95"][j]:.6f}',
                     'RIDGE_cap95':f'{rl["cap95"][j]:.6f}',
                     'HGBR_oof_RMSE':f'{rh["oof_rmse"][j]:.6f}',
                     'RF_oof_RMSE':f'{rr["oof_rmse"][j]:.6f}',
                     'RIDGE_oof_RMSE':f'{rl["oof_rmse"][j]:.6f}'})
    with (output_dir/'outcome_model_sensitivity_six_x_results.csv').open('w', newline='', encoding='utf-8') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    # Publishable 5-row compact comparison with no personal records or fold IDs.
    concise = [{key:r[key] for key in ('source','n','HGBR_uncapped','RF_uncapped','RIDGE_uncapped',
                'RF_minus_HGBR','RIDGE_minus_HGBR')} for r in rows]
    with (output_dir/'si_table_s8_outcome_model_sensitivity.csv').open('w', newline='', encoding='utf-8') as f:
        w=csv.DictWriter(f, fieldnames=list(concise[0]));w.writeheader();w.writerows(concise)
    # Original-provenance flag: all three use identical cross-fitted propensities and study cohort.
    out['primary_piped_minus_groundwater'] = float(
        fit['outcomes']['HGBR']['uncapped'][0] - fit['outcomes']['HGBR']['uncapped'][1])
    out['rf_piped_minus_groundwater'] = float(
        fit['outcomes']['RF']['uncapped'][0] - fit['outcomes']['RF']['uncapped'][1])
    out['ridge_piped_minus_groundwater'] = float(
        fit['outcomes']['RIDGE']['uncapped'][0] - fit['outcomes']['RIDGE']['uncapped'][1])
    (output_dir/'outcome_model_sensitivity_six_x_provenance.json').write_text(
        json.dumps(out, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('\nSOURCE | HGBR uncapped | RF uncapped | Ridge uncapped | RF−HGBR | Ridge−HGBR')
    for r in rows:
        print(' | '.join(r[key] for key in ['source','HGBR_uncapped','RF_uncapped',
               'RIDGE_uncapped','RF_minus_HGBR','RIDGE_minus_HGBR']))
    print('\nPIPED−GROUNDWATER: HGBR={:.6f}; RF={:.6f}; Ridge={:.6f}'.format(
        out['primary_piped_minus_groundwater'],out['rf_piped_minus_groundwater'],
        out['ridge_piped_minus_groundwater']))
    print('\nRESULT: PASS (original HGBR verified; E02 '+agreement['level']+' diagnostics verified)')
    print('Aggregate outputs:', output_dir)
    print('No household-level data, folds or predictions were written.')


if __name__ == '__main__':
    main()
