#!/usr/bin/env python3
"""WRR E28: learner prediction comparisons, nested HGBR tuning, working AIPW SE.

Empirical original dataset required; use original E04 cohort and E27 verified methods.
The script writes aggregate statistics only; OOF predictions, household records,
bootstrap indices, scores and fold membership remain in memory and are discarded.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import platform
from collections import Counter
from pathlib import Path

import numpy as np
import sklearn
import pandas as pd
import scipy
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import KFold, StratifiedKFold

BASE = Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('e27_extended', BASE/'E27_extended_in_memory.py')
e27=importlib.util.module_from_spec(spec)
spec.loader.exec_module(e27)

SEED=e27.SEED
LABELS=e27.LABELS
LABELS_PRETTY={'PIPED':'Municipal piped','GROUNDWATER':'Groundwater','SURFACE':'Surface water','TANKER':'Tanker delivery','COMMERCIAL':'Commercial water'}
REFERENCE_E27_MEANS={
  'HGBR':[74.846741,75.448523,76.365256,71.530038,72.298656],
  'RF':[74.854141,75.509683,76.136675,71.152282,72.272382],
  'RIDGE':[74.809342,75.483137,76.101096,71.124396,72.093996],
}
REFERENCE_E27_RMSE={
  'HGBR':[16.000470,15.737480,15.846393,14.099692,19.437763],
  'RF':[15.916857,15.554495,15.152420,13.594203,18.629721],
  'RIDGE':[15.873804,15.463228,15.309313,13.457539,18.149820],
}
# Fixed, compact candidate grid; inner-fold training only. This is a sensitivity,
# not a search for the globally optimal production learner.
GRID=[
  {'name':'original',    'learning_rate':.05,'max_iter':200,'max_leaf_nodes':7, 'max_depth':3,'min_samples_leaf':10,'l2_regularization':1.0},
  {'name':'conservative','learning_rate':.05,'max_iter':150,'max_leaf_nodes':5, 'max_depth':3,'min_samples_leaf':15,'l2_regularization':5.0},
  {'name':'flexible',    'learning_rate':.05,'max_iter':300,'max_leaf_nodes':15,'max_depth':4,'min_samples_leaf':10,'l2_regularization':.5},
  {'name':'slow_growth', 'learning_rate':.03,'max_iter':400,'max_leaf_nodes':7, 'max_depth':3,'min_samples_leaf':10,'l2_regularization':1.0},
]


def learner_from_cfg(cfg):
    return HistGradientBoostingRegressor(
        loss='squared_error',early_stopping=False,random_state=SEED,
        learning_rate=cfg['learning_rate'],max_iter=cfg['max_iter'],
        max_leaf_nodes=cfg['max_leaf_nodes'],max_depth=cfg['max_depth'],
        min_samples_leaf=cfg['min_samples_leaf'],l2_regularization=cfg['l2_regularization'])


def do_tuned_hgbr(X,y,scores, verbose=True, cfgs=GRID, inner_folds=3):
    n=len(y); k=len(LABELS)
    oof=np.full((n,k),np.nan,dtype=float)
    picks={lab:Counter() for lab in LABELS}
    # Each inner CV is confined to its source-specific subset of the outer training fold.
    grid_details=[]
    for fold,(train,test) in enumerate(StratifiedKFold(n_splits=5,shuffle=True,random_state=SEED).split(X,y),start=1):
        if verbose: print(f'Nested HGBR selection — outer fold {fold}/5',flush=True)
        for j,label in enumerate(LABELS):
            src=train[y[train]==j+1]
            if len(src)<inner_folds*20: raise RuntimeError(f'Unexpectedly small inner training cohort for {label}')
            kf=KFold(n_splits=inner_folds,shuffle=True,random_state=SEED+100*fold+j)
            split_cache=list(kf.split(src))
            model_scores=[]
            for cfg in cfgs:
                err_sum=0.0; n_hold=0
                for it,iv in split_cache:
                    model=learner_from_cfg(cfg)
                    model.fit(X[src[it]],scores[src[it]])
                    p=model.predict(X[src[iv]])
                    err_sum+=float(np.square(scores[src[iv]]-p).sum())
                    n_hold+=len(iv)
                model_scores.append(math.sqrt(err_sum/n_hold))
            selected=int(np.argmin(np.asarray(model_scores))) # predefined tie-break in order
            cfg=cfgs[selected]
            picks[label][cfg['name']]+=1
            fitted=learner_from_cfg(cfg)
            fitted.fit(X[src],scores[src])
            oof[test,j]=fitted.predict(X[test])
            grid_details.append({'source':label,'outer_fold':fold,'n_outer_train_source':len(src),
                                 'selected_name':cfg['name'], 'inner_cv_rmse':model_scores})
            if verbose: print(f'  {label:<12} {cfg["name"]:<13} inner RMSE={model_scores[selected]:.4f}',flush=True)
    if not np.isfinite(oof).all():raise RuntimeError('Tuned HGBR OOF predictions incomplete')
    if float(np.min(oof)) < -1e-6 or float(np.max(oof))>100+1e-6:
        raise RuntimeError('Tuned outcome predictions outside valid WSQS scale')
    # Source-level model selection counts only, not per-fold assignments in released output.
    selection_counts={label:{cfg['name']:int(picks[label][cfg['name']]) for cfg in cfgs} for label in LABELS}
    inner_error_averages={label:{cfg['name']:float(np.mean([r['inner_cv_rmse'][ix] for r in grid_details if r['source']==label]))
                                  for ix,cfg in enumerate(cfgs)} for label in LABELS}
    return oof,selection_counts,inner_error_averages


def oof_rmse(y,scores,preds):
    return np.asarray([np.sqrt(np.mean(np.square(scores[y==j+1]-preds[y==j+1,j])))
                       for j in range(len(LABELS))],dtype=float)


def boot_rmse_pairs(y,scores,predictions, bootstrap=2000, seed=SEED+2800):
    """Paired percentile intervals, resampling within each source with models held fixed.

    These are conditional diagnostics only: the full nuisance fitting procedure is
    not rerun for each bootstrap sample. Do not label as design-based inference.
    """
    names=['HGBR','RF','RIDGE','HGBR_TUNED']
    n=len(y); rng=np.random.default_rng(seed)
    results=[]
    bootstrap_sum={m:np.zeros(bootstrap,dtype=float) for m in names}
    error_sq={m:{} for m in names}
    for j,label in enumerate(LABELS):
        idx=np.flatnonzero(y==j+1); ng=len(idx)
        sampled=rng.integers(0,ng,size=(bootstrap,ng))
        sample_mean_sq={}
        for method in names:
            vec=np.square(scores[idx]-predictions[method][idx,j])
            error_sq[method][label]=vec
            sample_mean_sq[method]=np.mean(vec[sampled],axis=1)
            bootstrap_sum[method]+=ng*sample_mean_sq[method]
        baseline=np.sqrt(float(error_sq['HGBR'][label].mean()))
        for method in names[1:]:
            delta=float(np.sqrt(error_sq[method][label].mean())-baseline)
            draw=np.sqrt(sample_mean_sq[method])-np.sqrt(sample_mean_sq['HGBR'])
            lower,upper=np.quantile(draw,[0.025,0.975])
            results.append({'source':label,'n':ng,'comparison':method+' minus HGBR',
                            'HGBR_RMSE':baseline,'alternative_RMSE':float(np.sqrt(error_sq[method][label].mean())),
                            'RMSE_difference':delta, 'paired_bootstrap_lower':float(lower),
                            'paired_bootstrap_upper':float(upper),'bootstrap_replicates':bootstrap,
                            'conditional_on_fitted_models':True})
    baseline=float(np.sqrt(sum(float(np.sum(error_sq['HGBR'][lab])) for lab in LABELS)/n))
    for method in names[1:]:
        mse_total=sum(float(np.sum(error_sq[method][lab])) for lab in LABELS)/n
        delta=float(np.sqrt(mse_total)-baseline)
        boot_deltas=np.sqrt(bootstrap_sum[method]/n)-np.sqrt(bootstrap_sum['HGBR']/n)
        lower,upper=np.quantile(boot_deltas,[0.025,0.975])
        results.append({'source':'POOLED','n':n,'comparison':method+' minus HGBR',
                        'HGBR_RMSE':baseline,'alternative_RMSE':float(np.sqrt(mse_total)),
                        'RMSE_difference':delta,'paired_bootstrap_lower':float(lower),
                        'paired_bootstrap_upper':float(upper),'bootstrap_replicates':bootstrap,
                        'conditional_on_fitted_models':True})
    return results


def aipw_working_summary(y,scores,e,predictions,engine):
    labels=np.asarray(LABELS); zz=labels[np.asarray(y,dtype=int)-1]
    n=len(y)
    contrasts=[(i,j) for i in range(len(LABELS)) for j in range(i+1,len(LABELS))]
    means=[]; contrast_rows=[]; psi_by_method={}
    for method,pred in predictions.items():
        psi=np.asarray(engine.aipw_scores(zz,scores,e,pred),dtype=float)
        direct=pred+(zz[:,None]==labels[None,:])*(scores[:,None]-pred)/e
        if not np.allclose(psi,direct,atol=1e-8,rtol=0):raise RuntimeError(f'AIPW score equation failed for {method}')
        if not np.isfinite(psi).all():raise RuntimeError('Nonfinite AIPW scores')
        psi_by_method[method]=psi
        for j,label in enumerate(LABELS):
            values=psi[:,j]
            means.append({'learner':method,'source':label,'n':n,'mean_AIPW':float(values.mean()),
                         'working_IF_SE':float(np.std(values,ddof=1)/math.sqrt(n)),
                         'score_SD':float(np.std(values,ddof=1)),
                         'variance_method':'empirical centered AIPW scores, nuisance fits treated as fixed'})
        for j,k in contrasts:
            values=psi[:,j]-psi[:,k]
            contrast_rows.append({'learner':method,'contrast':f'{LABELS[j]}-minus-{LABELS[k]}',
                                  'n':n,'contrast_AIPW':float(values.mean()),
                                  'working_IF_SE':float(np.std(values,ddof=1)/math.sqrt(n)),
                                  'variance_method':'empirical centered contrast scores, nuisance fits treated as fixed'})
    # Comparison within same respondents. Delta in AIPW influence functions, with
    # original HGBR as reference; these are *working*, not design-based SEs.
    change_rows=[]
    for method in ['RF','RIDGE','HGBR_TUNED']:
        dp=psi_by_method[method]-psi_by_method['HGBR']
        for j,label in enumerate(LABELS):
            v=dp[:,j]
            change_rows.append({'comparison':method+' minus HGBR','estimand':label,
                                'AIPW_estimate_change':float(np.mean(v)),
                                'working_IF_SE_of_change':float(np.std(v,ddof=1)/math.sqrt(n))})
        for j,k in contrasts:
            v=dp[:,j]-dp[:,k]
            change_rows.append({'comparison':method+' minus HGBR','estimand':f'{LABELS[j]}-minus-{LABELS[k]}',
                                'AIPW_estimate_change':float(np.mean(v)),
                                'working_IF_SE_of_change':float(np.std(v,ddof=1)/math.sqrt(n))})
    return means,contrast_rows,change_rows


def save_csv(path,rows):
    if not rows:return
    keys=list(rows[0])
    with path.open('w',newline='',encoding='utf-8') as f:
        wr=csv.DictWriter(f,fieldnames=keys)
        wr.writeheader();wr.writerows(rows)


def validate_e27_anchors(y,scores,predictions,fit):
    for method in ['HGBR','RF','RIDGE']:
        means=np.asarray(fit['outcomes'][method]['uncapped'])
        rmses=oof_rmse(y,scores,predictions[method])
        ref_means=np.asarray(REFERENCE_E27_MEANS[method])
        ref_rmse=np.asarray(REFERENCE_E27_RMSE[method])
        deltas=(means-ref_means,rmses-ref_rmse)
        if max(float(np.max(abs(d))) for d in deltas)>.005:
            raise RuntimeError(f'E27 user-verified anchor mismatch {method}: {deltas}')
        print(f'  {method:<6} max mean deviation {max(abs(deltas[0])):.7f}; max RMSE deviation {max(abs(deltas[1])):.7f}',flush=True)


def validate_saved_e27_stage(reference_csv, reference_code, fit):
    """Verify this E28 replay against the actual successful E27 stage on user's Mac."""
    if not reference_csv.is_file():
        raise RuntimeError(f'Original successful E27 aggregate file missing: {reference_csv}')
    if not reference_code.is_file():
        raise RuntimeError(f'Original successful E27 code missing: {reference_code}')
    if e27.sha256(reference_code)!=e27.sha256(BASE/'E27_REFERENCE_UNMODIFIED.py'):
        raise RuntimeError('E27 stage code does not match the archived unmodified Replay4 runner')
    rows=list(csv.DictReader(reference_csv.open(newline='',encoding='utf-8')))
    if [r['source'] for r in rows]!=LABELS:
        raise RuntimeError('E27 stage source order changed')
    nchecks=0; maximum=0.0
    column_pairs={'uncapped':'uncapped','cap99':'cap99','cap95':'cap95','oof_rmse':'oof_RMSE'}
    for j,row in enumerate(rows):
        if int(row['n'])!=e27.REFERENCE_N[j]:
            raise RuntimeError('Historical E27 source group count mismatch')
        for model in ('HGBR','RF','RIDGE'):
            for key,col in column_pairs.items():
                got=float(fit['outcomes'][model][key][j]); expected=float(row[f'{model}_{col}'])
                delta=abs(got-expected)
                maximum=max(delta,maximum);nchecks+=1
                if not np.isfinite(delta) or delta>.005:
                    raise RuntimeError(f'E27 archived aggregate differs: {model}/{LABELS[j]}/{key}: {delta}')
    print(f'PASS: {nchecks} empirical E27 aggregate values reproduced; largest difference={maximum:.9f}',flush=True)
    return {'checks':nchecks,'max_delta':maximum,
            'e27_original_code_sha256':e27.sha256(reference_code),
            'e27_results_csv_sha256':e27.sha256(reference_csv)}


def create_report(out_dir,paired,means,contrasts,changes,selection,inner_avg,verification,run_cfg):
    def lookup_ms(method,source):
        return next(x for x in means if x['learner']==method and x['source']==source)
    def lookup_co(method,con):
        return next(x for x in contrasts if x['learner']==method and x['contrast']==con)
    def paired_lookup(method,source):
        return next(x for x in paired if x['comparison']==method+' minus HGBR' and x['source']==source)
    pretty={'HGBR':'HGBR original','RF':'Random Forest','RIDGE':'Ridge','HGBR_TUNED':'HGBR nested tuned'}
    lines=[]
    def line(x=''):lines.append(str(x))
    line('='*100)
    line('WRR EMPIRICAL 28 — PREDICTION, NESTED TUNING, AND AIPW WORKING VARIABILITY')
    line('='*100)
    line('STATUS: EMPIRICAL EXTENSION. E27 aggregates, E04 primary, and E02 numerical checks must all pass.')
    line('Original six-X cohort: N=3066 | 5 source strata | 33 one-hot indicators | seed=20261022')
    line('Estimated comparison: HGBR original vs RF 500 trees, Ridge alpha=10, and HGBR tuned within training folds')
    line('Bootstrap: '+str(run_cfg['bootstrap_replicates'])+' household-paired resamples within original source groups; NO learner refits')
    line('Warning: paired bootstrap intervals are descriptive conditional diagnostics, not population-level/design-based confidence intervals.')
    line('Warning: AIPW working SEs use centered empirical score variance; not publication-certified inference.')
    line('Environment: Python '+platform.python_version()+' | NumPy '+np.__version__+' | pandas '+pd.__version__+' | SciPy '+scipy.__version__+' | sklearn '+sklearn.__version__)
    line()
    line('1. PAIRED OUT-OF-FOLD PREDICTION RMSE, DIFFERENCES vs ORIGINAL HGBR (WSQS POINTS)')
    line(f"{'Source':<14}{'HGBR':>10}{'RF':>10}{'Ridge':>10}{'Tuned':>10}{'RF Δ':>10}{'Ridge Δ':>10}{'Tune Δ':>10}")
    for src in LABELS+['POOLED']:
        a=paired_lookup('RF',src); b=paired_lookup('RIDGE',src); c=paired_lookup('HGBR_TUNED',src)
        line(f'{src:<14}{a["HGBR_RMSE"]:>10.4f}{a["alternative_RMSE"]:>10.4f}{b["alternative_RMSE"]:>10.4f}{c["alternative_RMSE"]:>10.4f}'+
             f'{a["RMSE_difference"]:>+10.4f}{b["RMSE_difference"]:>+10.4f}{c["RMSE_difference"]:>+10.4f}')
    line('\n2. PAIRED BOOTSTRAP 95% PERCENTILE INTERVALS FOR RMSE DIFFERENCE (ALTERNATIVE − HGBR)')
    line(f"{'Source':<14}{'Learner':<16}{'Difference':>12}{'Lower':>12}{'Upper':>12}")
    for src in LABELS+['POOLED']:
        for method in ('RF','RIDGE','HGBR_TUNED'):
            r=paired_lookup(method,src)
            line(f'{src:<14}{method:<16}{r["RMSE_difference"]:>+12.5f}{r["paired_bootstrap_lower"]:>+12.5f}{r["paired_bootstrap_upper"]:>+12.5f}')
    line('\n3. OUTER-TRAINING-FOLD HGBR NESTED TUNING (3-FOLD INNER CV)')
    line(f"{'Source':<14}{'Original':>12}{'Conserv':>12}{'Flexible':>12}{'Slow growth':>12}")
    for src in LABELS:
        x=selection[src]
        line(f'{src:<14}{x.get("original",0):>12d}{x.get("conservative",0):>12d}{x.get("flexible",0):>12d}{x.get("slow_growth",0):>12d}')
    line('Each row sums to five outer training folds; inner CV uses only source-specific outer training observations.')
    line('\n4. AIPW STANDARDIZED MEANS AND WORKING IF STANDARD ERRORS (NOT FORMAL INFERENCE)')
    line(f"{'Source':<14}{'Method':<18}{'Mean':>12}{'Working SE':>14}{'Shift':>12}")
    for src in LABELS:
        ref=lookup_ms('HGBR',src)['mean_AIPW']
        for method in ('HGBR','RF','RIDGE','HGBR_TUNED'):
            r=lookup_ms(method,src)
            line(f'{src:<14}{pretty[method]:<18}{r["mean_AIPW"]:>12.6f}{r["working_IF_SE"]:>14.6f}{r["mean_AIPW"]-ref:>+12.6f}')
    line('\n5. PIPED−GROUNDWATER STANDARDIZED CONTRAST (WSQS POINTS)')
    line(f"{'Model':<18}{'Contrast':>14}{'Working SE':>14}{'Shift from HGBR':>18}")
    ref=lookup_co('HGBR','PIPED-minus-GROUNDWATER')['contrast_AIPW']
    for method in ('HGBR','RF','RIDGE','HGBR_TUNED'):
        x=lookup_co(method,'PIPED-minus-GROUNDWATER')
        line(f'{pretty[method]:<18}{x["contrast_AIPW"]:>14.6f}{x["working_IF_SE"]:>14.6f}{x["contrast_AIPW"]-ref:>+18.6f}')
    line('\n6. CROSS-LEARNER WORKING IF SE OF AIPW ESTIMATE DIFFERENCE')
    for method in ('RF','RIDGE','HGBR_TUNED'):
        r=next(x for x in changes if x['comparison']==method+' minus HGBR' and x['estimand']=='PIPED-minus-GROUNDWATER')
        line(f' {method:<12}: shift={r["AIPW_estimate_change"]:+.6f}; working SE shift={r["working_IF_SE_of_change"]:.6f}')
    line()
    line('SOURCE INTEGRITY: E04 code and AIPW engine SHA verified; E02, E04 and E27 parity/near-parity checks passed.')
    line('ORIGINAL DATA HANDLING: no respondent records, IDs, OOF predictions, or bootstrap indices written to disk.')
    text='\n'.join(lines)+'\n'
    (out_dir/'RESULTS_TO_PASTE_IN_CHAT.txt').write_text(text,encoding='utf-8')
    print(text,flush=True)

    # S9: one compact supplemental table with pooled prediction comparisons,
    # source-level working IF SEs and a reference contrast; no inference claim.
    sr=['**Table S9. Outcome-learner diagnostics for six-characteristic standardized WSQS estimates.**','',
       '**Panel A. Out-of-fold predictive accuracy (RMSE, WSQS points)**','',
       '| Source | HGBR | Random Forest | Ridge | Nested-tuned HGBR |','|---|---:|---:|---:|---:|']
    for src in LABELS+['POOLED']:
        a=paired_lookup('RF',src);b=paired_lookup('RIDGE',src);c=paired_lookup('HGBR_TUNED',src)
        sr.append(f'| {LABELS_PRETTY.get(src,"Pooled across sources")} | {a["HGBR_RMSE"]:.2f} | {a["alternative_RMSE"]:.2f} | {b["alternative_RMSE"]:.2f} | {c["alternative_RMSE"]:.2f} |')
    sr+=['','**Panel B. Working AIPW influence-function standard errors (WSQS points)**','',
         '| Estimate | HGBR | Random Forest | Ridge | Nested-tuned HGBR |','|---|---:|---:|---:|---:|']
    for src in LABELS:
        sr.append('| '+LABELS_PRETTY[src]+' | '+' | '.join(f'{lookup_ms(m,src)["working_IF_SE"]:.3f}' for m in ('HGBR','RF','RIDGE','HGBR_TUNED'))+' |')
    sr.append('| Piped − groundwater | '+' | '.join(f'{lookup_co(m,"PIPED-minus-GROUNDWATER")["working_IF_SE"]:.3f}' for m in ('HGBR','RF','RIDGE','HGBR_TUNED'))+' |')
    sr+=['', '**Note.** All specifications use the same original N = 3,066, five-fold outer splits, six adjustment variables (33 one-hot indicators), propensity nuisance estimates and uncapped AIPW score. The original HGBR is fixed at 200 iterations, maximum seven leaves, depth three, learning rate 0.05, minimum leaf size 10 and L2 regularization 1.0. The nested-tuned HGBR selects one of four prespecified settings using three-fold inner validation within each source-specific outer training set. Random Forest uses 500 trees with minimum leaf size 10; Ridge uses L2 penalty 10. Lower RMSE indicates better OOF prediction. Working standard errors use the empirical variance of centered AIPW scores divided by N; they are descriptive diagnostics and are not survey-design-based or publication-certified confidence intervals. Paired, source-stratified bootstrap RMSE-difference intervals and full tuning details are available in the stage CSVs. E02 numerical near-parity, not exact historical parity, was established.','']
    (out_dir/'SUPPLEMENTARY_TABLE_S9_READY_TO_PASTE.md').write_text('\n'.join(sr),encoding='utf-8')
    ar=paired_lookup('RF','POOLED'); br=paired_lookup('RIDGE','POOLED'); cr=paired_lookup('HGBR_TUNED','POOLED')
    hm=lookup_co('HGBR','PIPED-minus-GROUNDWATER')['contrast_AIPW']
    tm=lookup_co('HGBR_TUNED','PIPED-minus-GROUNDWATER')['contrast_AIPW']
    # Avoid rhetorical claim of improvement if tuned model is worse.
    tex=(
      '**Supplementary Methods (proposed).** We compared paired out-of-fold prediction errors '
      'for the original histogram gradient boosting model, Random Forest, ridge regression '
      'and a nested-tuned histogram gradient boosting model. Tuning used a prespecified '
      'four-configuration grid with three-fold validation confined to each source-specific outer training fold; '
      'the original five-fold cross-fitting, propensity model, and uncapped AIPW estimator were unchanged. '
      'We also calculated working standard errors from the empirical variance of the cross-fitted AIPW scores '
      'and examined paired, source-stratified bootstrap intervals for differences in prediction RMSE. '
      'These uncertainty diagnostics are conditional on the fitted nuisance models and do not constitute design-based inference.\n\n'
      '**Results (proposed; review before insertion).** In pooled out-of-fold prediction, the original HGBR RMSE '
      f'was {ar["HGBR_RMSE"]:.2f} WSQS points, compared with {ar["alternative_RMSE"]:.2f} for Random Forest, '
      f'{br["alternative_RMSE"]:.2f} for ridge regression and {cr["alternative_RMSE"]:.2f} for nested-tuned HGBR. '
      'The piped–groundwater standardized contrast changed from '
      f'{hm:+.2f} under the original HGBR to {tm:+.2f} under nested-tuned HGBR (Table S9). '
      'Working AIPW standard-error comparisons are reported in Table S9.\n'
    )
    (out_dir/'PROPOSED_SI_AND_RESULTS_TEXT.md').write_text(tex,encoding='utf-8')


def main():
    p=argparse.ArgumentParser(description='WRR E28 paired OOF, HGBR tuning, and working AIPW IF diagnostics')
    for flag in ('e04-code','decoder','g7-report','workbook','ref-root','package-root','e04-output','output-dir'):
        p.add_argument('--'+flag,required=True,type=Path)
    p.add_argument('--e27-reference-csv',type=Path,required=True)
    p.add_argument('--e27-reference-code',type=Path,required=True)
    p.add_argument('--bootstrap',type=int,default=2000)
    args=p.parse_args()
    if args.bootstrap<200 or args.bootstrap>10000:
        p.error('--bootstrap must be between 200 and 10000')
    if e27.sha256(args.e04_code)!=e27.E04_SHA:raise RuntimeError('Original E04 code SHA mismatch')
    enginepath=args.package_root/'src/water_standardize/core.py'
    if e27.sha256(enginepath)!=e27.ENGINE_SHA:raise RuntimeError('Frozen E04 AIPW engine SHA mismatch')
    output=args.output_dir.resolve()
    if output==args.e04_output.resolve() or args.e04_output.resolve() in output.parents:
        raise RuntimeError('Cannot write into original E04 stage')
    ns=e27.load_original_e04_cohort(args.e04_code,args.decoder,args.g7_report,args.workbook,args.ref_root,args.package_root,args.e04_output)
    X,y,scores=(np.asarray(ns[v]) for v in ('X','y','scores'))
    if X.shape!=(3066,33) or [int((y==j).sum()) for j in range(1,6)]!=e27.REFERENCE_N:
        raise RuntimeError('Cohort/group shape unexpected')
    if list(ns['QIDS'])!=e27.QIDS: raise RuntimeError('Six adjustment questions changed')
    import sys
    sys.path.insert(0,str(args.package_root/'src'))
    from water_standardize import core as engine
    print('\nREPRODUCING HISTORICAL E27 INPUTS AND FITS (NO NEW OUTPUT ACCEPTED YET)',flush=True)
    fit=e27.fit_nuisance_models(X,y,scores,engine,return_oof=True)
    h=fit['outcomes']['HGBR']['uncapped']
    if not np.allclose(h,e27.REFERENCE_MEANS,atol=.005,rtol=0):
        raise RuntimeError('Frozen E04 HGBR means failed reproduction')
    ref02=json.loads((args.e04_output/'_empirical02_reference.json').read_text())
    check02=e27.evaluate_e02_numerical_agreement(X,y,fit['e'],fit['logloss'],ref02)
    check04=e27.evaluate_e04_aggregate_agreement(args.e04_output,fit['outcomes'])
    predictions=fit['predictions_in_memory_only']
    print('\nCOMPARING WITH COMPLETED EMPIRICAL 27 AGGREGATES:',flush=True)
    validate_e27_anchors(y,scores,predictions,fit)
    saved_e27_validation=validate_saved_e27_stage(args.e27_reference_csv,args.e27_reference_code,fit)
    print('PASS: E04/E02/E27 reference checks',flush=True)
    print('\nRUNNING NESTED HGBR TRAINING-ONLY TUNING',flush=True)
    tune,selected,inner_error=do_tuned_hgbr(X,y,scores)
    predictions['HGBR_TUNED']=tune
    print('\nRUNNING PAIRED OUT-OF-FOLD RMSE BOOTSTRAPS',flush=True)
    paired=boot_rmse_pairs(y,scores,predictions,bootstrap=args.bootstrap)
    print('COMPUTING WORKING AIPW INFLUENCE-FUNCTION VARIABILITY',flush=True)
    means,contrast,changes=aipw_working_summary(y,scores,fit['e'],predictions,engine)
    # Independent anchor on AIPW mean computed by the new generic diagnostic.
    for j,label in enumerate(LABELS):
        got=next(v['mean_AIPW'] for v in means if v['learner']=='HGBR' and v['source']==label)
        if abs(got-float(h[j]))>1e-8:raise RuntimeError('E28 AIPW baseline not identical to E27 score')
    provenance={
        'status':'empirical_checked','scope':'E28 conditional diagnostics, no design-based inference',
        'historical_model_status':'original E04 HGBR unchanged; E27 RF/Ridge unchanged',
        'n':int(len(y)),'n_group':{lab:int((y==j+1).sum()) for j,lab in enumerate(LABELS)},
        'feature_count':int(X.shape[1]),'adjustment_qids':e27.QIDS,
        'seed_outer':SEED,'n_outer_folds':5,'n_inner_folds':3,'bootstrap_replicates':args.bootstrap,
        'bootstrap_resampling':'paired household resampling within source group; nuisance fits held fixed',
        'tuning_selection':'per source and outer fold, three-fold inner RMSE, four prespecified configurations',
        'hgb_tuning_grid':GRID,'hgb_tuning_selection_counts':selected,'hgb_mean_inner_cv_rmse':inner_error,
        'historical_e02_level':check02['level'],'historical_e02_logloss_delta':check02['difference'],
        'historical_e04_max_abs_difference':float(max(abs(v['difference']) for gg in check04.values() for v in gg.values())),
        'historical_e27_anchor_max_tolerance':.005,
        'archived_successful_e27_comparison':saved_e27_validation,
        'source_sha256_original_e04':e27.E04_SHA,'source_sha256_original_engine':e27.ENGINE_SHA,
        'reference_e27_code_sha256':e27.sha256(BASE/'E27_REFERENCE_UNMODIFIED.py'),
        'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,
        'scipy':scipy.__version__,'scikit_learn':sklearn.__version__,
        'restriction':'No respondent-level files, OOF predictions, folds or bootstrap samples exported',
        'inference_limit':'working empirical AIPW score SE, ignores nuisance-fitting uncertainty, nonprobability sample',
    }
    output.mkdir(parents=True,exist_ok=True)
    save_csv(output/'PAIRED_OOF_RMSE_BOOTSTRAP.csv',paired)
    save_csv(output/'AIPW_WORKING_IF_SE_MEANS.csv',means)
    save_csv(output/'AIPW_WORKING_IF_SE_CONTRASTS.csv',contrast)
    save_csv(output/'AIPW_CROSS_LEARNER_WORKING_IF_DIFFERENCES.csv',changes)
    save_csv(output/'HGBR_INNER_TUNING_SELECTION_COUNTS.csv',[
       {'source':lab,**selected[lab],**{'mean_inner_rmse_'+key:val for key,val in inner_error[lab].items()}} for lab in LABELS
    ])
    (output/'E28_VERIFIED_AGGREGATE_PROVENANCE.json').write_text(json.dumps(provenance,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    create_report(output,paired,means,contrast,changes,selected,inner_error,check02,{'bootstrap_replicates':args.bootstrap})
    print('\nE28 VERIFIED EMPIRICAL AGGREGATE RESULTS: PASS',flush=True)
    print('All shareable outputs are aggregate-only; no household-level predictions or folds written.',flush=True)


if __name__=='__main__':
    main()
