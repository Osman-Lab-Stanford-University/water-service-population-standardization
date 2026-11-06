"""E20 independent oracle-nuisance holdout validation of two FIXED candidate
full-grid band methods. No study data, fitted nuisance claims or clustered design.
New outcome DGP: correlated covariates and QUANTIZED 0-100 scores with a
logistic latent variable. The fixed inference families and methods are the
same as the E19 candidates. No method is selected from this experiment.
"""
import os,json,sys,math
from pathlib import Path
from collections import Counter
import numpy as np
from scipy.special import expit,logit

OUT=Path(sys.argv[1]);OUT.mkdir(parents=True,exist_ok=True)
N=3066;K=5;GRID=np.linspace(0.,100.,41)
REPS_NULL=int(os.environ.get('E20_NULL_REPS','240'))
REPS_ALT=int(os.environ.get('E20_ALT_REPS','120'))
MULT=int(os.environ.get('E20_MULTIPLIER_DRAWS','1000'))
assert 2<=REPS_NULL<=500 and 2<=REPS_ALT<=500 and 100<=MULT<=5000
SEED=20261030
rng=np.random.default_rng(SEED)
LABELS=['PIPED','GROUNDWATER','SURFACE','TANKER','COMMERCIAL']
PAIRS={'PIPED_GROUNDWATER':(0,1),'PIPED_TANKER':(0,3)}
FAMILIES={
 'FULL_INTERIOR_39':np.arange(1,40),
 'CENTRAL_20_TO_80':np.arange(8,33),
 'REGISTERED_50_60_70':np.array([20,24,28])}
FINE=3.125
SCALE=.75
SHIFT_ALTERNATIVE=.28
# Four cells have NONUNIFORM known mass, unlike the E17-19 DGP.
CELLS=np.array([[0,0],[0,1],[1,0],[1,1]],dtype=int)
P_B=.57
P_C_IF_B=np.array([.22,.52])
CELL_WEIGHTS=np.array([
 (P_B if b else 1-P_B)*(P_C_IF_B[b] if c else 1-P_C_IF_B[b])
 for b,c in CELLS],dtype=float)
assert np.isclose(CELL_WEIGHTS.sum(),1.)


def true_prop(b,c):
    b=np.asarray(b,dtype=float);c=np.asarray(c,dtype=float)
    eta=(np.log(np.array([.40,.39,.08,.045,.085]))[None,:]
      + b[:,None]*np.array([.48,-.38,.20,.06,-.16])[None,:]
      + c[:,None]*np.array([-.31,.28,.08,-.16,.06])[None,:]
      + (b*c)[:,None]*np.array([.11,-.07,.04,.10,-.02])[None,:])
    e=np.exp(eta-eta.max(axis=1,keepdims=True))
    return e/e.sum(axis=1,keepdims=True)


def true_mu(b,c,shift):
    base=.70+.52*b-.43*c+.35*b*c
    # Under the null groups 0,1,3 have EXACTLY equal conditional laws.
    offsets=np.array([0.,0.,.12,0.,-.08]);offsets[0]+=shift
    return base[:,None]+offsets[None,:]


def true_cdf(mu):
    # Observed Y is a latent continuous 0-100 score ROUNDED to a
    # grid with 3.125-point steps (analogous to 8-item severity scores).
    # For each threshold t, last permitted discrete value is floor(t/FINE)*FINE.
    # Its latent rounding boundary lies FINE/2 higher.
    upper=(np.floor(GRID/FINE)+.5)*FINE
    upper=np.clip(upper,1e-12,100.-1e-12)
    out=expit((logit(upper/100.)[None,:,None]-mu[:,None,:])/SCALE)
    out[:,-1,:]=1.
    if np.any(np.diff(out,axis=1)<-1e-13):
        raise RuntimeError('Synthetic conditional CDF not monotone')
    return out


def target(shift):
    c=true_cdf(true_mu(CELLS[:,0],CELLS[:,1],shift))
    return np.einsum('c,ctk->tk',CELL_WEIGHTS,c)


def rate(a):
    a=np.asarray(a,dtype=float);v=float(a.mean())
    return {'estimate':v,'monte_carlo_se':float(np.sqrt(v*(1-v)/len(a)))}


def new_record():
    return {'tanker_n':[],'tanker_zero_le_2_5':0,'tanker_zero_le_5':0,
      'families':{p:{f:{key:[] for key in ('studentized_cover',
                           'absolute_cover','studentized_reject','absolute_reject')}
                    for f in FAMILIES} for p in PAIRS},
      'pointwise70':{p:[] for p in PAIRS}}


def test_pair(psi,delta,sc,rec,pair):
    est=psi.mean(axis=0)
    centered=psi-est[None,:]
    sd=np.std(psi,axis=0,ddof=1)
    se=sd/np.sqrt(N)
    # Threshold 100 is deterministically degenerate for AIPW contrasts.
    if not np.isfinite(se).all() or np.any(se[1:40]<=1e-12):
        raise RuntimeError('Nonpositive synthetic interior standard error')
    idx0=FAMILIES['FULL_INTERIOR_39']
    standardized=centered[:,idx0]/sd[idx0][None,:]
    corr=standardized.T@standardized/(N-1)
    corr=(corr+corr.T)/2
    val,vec=np.linalg.eigh(corr)
    if np.min(val)<-1e-7: raise RuntimeError('Invalid covariance matrix')
    L=vec*np.sqrt(np.maximum(val,0.))[None,:]
    draws=rng.standard_normal((MULT,len(idx0)))@L.T
    # Draws are approximately unit-standard-normal; multiply by
    # working threshold-specific SE to form unstudentized deviations.
    for fam,idx in FAMILIES.items():
        positions=idx-1
        Z=np.abs(draws[:,positions]);err=np.abs(est[idx]-delta[idx])
        qt=float(np.quantile(np.max(Z,axis=1),.95,method='higher'))
        qu=float(np.quantile(np.max(Z*se[idx][None,:],axis=1),.95,method='higher'))
        fam_record=rec['families'][pair][fam]
        fam_record['studentized_cover'].append(int(np.max(err/se[idx])<=qt))
        fam_record['absolute_cover'].append(int(np.max(err)<=qu))
        fam_record['studentized_reject'].append(int(np.max(np.abs(est[idx]/se[idx]))>qt))
        fam_record['absolute_reject'].append(int(np.max(np.abs(est[idx]))>qu))
    rec['pointwise70'][pair].append(float(est[28]-delta[28]))


def run_one(shift,delta_truth,rec):
    b=rng.binomial(1,P_B,size=N)
    c=rng.binomial(1,P_C_IF_B[b],size=N)
    p=true_prop(b,c)
    z=(rng.random(N)[:,None]>np.cumsum(p,axis=1)).sum(axis=1)
    if not set(range(K)).issubset(set(z)):
        raise RuntimeError('Missing synthetic regime')
    tanker=(z==3);tanker_n=int(tanker.sum())
    if tanker_n<20:raise RuntimeError('Unexpected tiny tanker group')
    rec['tanker_n'].append(tanker_n)
    mu=true_mu(b,c,shift)
    latent=100*expit(mu[np.arange(N),z]+SCALE*rng.logistic(size=N))
    Y=np.clip(np.rint(latent/FINE)*FINE,0,100)
    event=(Y[:,None]<=GRID[None,:]).astype(float)
    rec['tanker_zero_le_2_5']+=int(not np.any(event[tanker,1]))
    rec['tanker_zero_le_5']+=int(not np.any(event[tanker,2]))
    q=true_cdf(mu)
    mask=(z[:,None]==np.arange(K)[None,:])
    phi=q+mask[:,None,:]*(event[:,:,None]-q)/p[:,None,:]
    if not np.isfinite(phi).all(): raise RuntimeError('Nonfinite oracle AIPW')
    for pair,(a,d) in PAIRS.items():
        test_pair(phi[:,:,a]-phi[:,:,d],delta_truth[:,a]-delta_truth[:,d],shift,rec,pair)


def finalize(rec,R):
    pairs={}
    for pair in PAIRS:
        fs={}
        for fam in FAMILIES:
            raw=rec['families'][pair][fam]
            fs[fam]={
                'studentized_joint_coverage':rate(raw['studentized_cover']),
                'absolute_joint_coverage':rate(raw['absolute_cover']),
                'studentized_equality_rejection':rate(raw['studentized_reject']),
                'absolute_equality_rejection':rate(raw['absolute_reject']),
            }
        pairs[pair]={'families':fs,'bias_at_70_pp':100*float(np.mean(rec['pointwise70'][pair]))}
    return {'replications':R,'tanker_median_n':float(np.median(rec['tanker_n'])),
       'share_tanker_no_events_le_2_5':rec['tanker_zero_le_2_5']/R,
       'share_tanker_no_events_le_5':rec['tanker_zero_le_5']/R,
       'pairs':pairs}

results={}
for scenario,shift,R in [('null',0.,REPS_NULL),('alternative',SHIFT_ALTERNATIVE,REPS_ALT)]:
    truth=target(shift)
    if scenario=='null':
        assert max(np.max(np.abs(truth[:,0]-truth[:,g])) for g in (1,3))<1e-12
    else:
        assert abs(truth[28,0]-truth[28,1])>.01
    rec=new_record()
    print('SCENARIO_START',scenario,'REPS',R,flush=True)
    for i in range(R):
        run_one(shift,truth,rec)
        if (i+1)%40==0 or i+1==R:
            print('SIMULATION_PROGRESS',scenario,i+1,'OF',R,flush=True)
    result=finalize(rec,R);results[scenario]=result
    print('SCENARIO',scenario,'REPS',R,'TANKER_MEDIAN_N',result['tanker_median_n'],
          'ZERO_LE_2_5',round(result['share_tanker_no_events_le_2_5'],4),
          'ZERO_LE_5',round(result['share_tanker_no_events_le_5'],4),flush=True)
    for pair in PAIRS:
        for family in FAMILIES:
            r=result['pairs'][pair]['families'][family]
            print('PAIR',pair,'FAMILY',family,
                  'STUDENTIZED_COVER',round(r['studentized_joint_coverage']['estimate'],4),
                  'STUDENTIZED_COVER_MCSE',round(r['studentized_joint_coverage']['monte_carlo_se'],4),
                  'ABSOLUTE_COVER',round(r['absolute_joint_coverage']['estimate'],4),
                  'ABSOLUTE_COVER_MCSE',round(r['absolute_joint_coverage']['monte_carlo_se'],4),
                  'STUDENTIZED_REJECT',round(r['studentized_equality_rejection']['estimate'],4),
                  'ABSOLUTE_REJECT',round(r['absolute_equality_rejection']['estimate'],4),flush=True)
    for pair in PAIRS:
        print('BIAS_AT_70',scenario,pair,'PP',round(result['pairs'][pair]['bias_at_70_pp'],4),flush=True)

report={
  'phase':'EMPIRICAL_20_INDEPENDENT_ORACLE_FULL_GRID_HOLDOUT',
  'design':'NEW correlated-binary-X, five-regime, quantized 3.125-step WSQS synthetic DGP; fixed regime CDF shifts',
  'comparison_candidates':['studentized_max_t','unstudentized_max_abs'],
  'not_a_prespecified_confirmatory_choice_of_method':True,
  'three_families_predeclared':{name:[float(GRID[i]) for i in ids] for name,ids in FAMILIES.items()},
  'full_41_grid_scored':True,
  'full_grid_inference_39_fixed_interior_thresholds':True,
  'n_per_replication':N,'replicate_seed':SEED,'multiplier_draws':MULT,
  'scenarios':results,
  'oracle_nuisance_only':True,'fitted_nuisance_models_validated_full_grid':False,
  'recruitment_clusters_simulated':False,'restricted_household_data_read':False,
  'case_level_data_saved':False,'publication_inference_certified':False,
  'public_release_authorized':False,'previous_stages_modified':False}
(OUT/'EMPIRICAL_20_REPORT.json').write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8')
print('EMPIRICAL_20_RESULTS')
print('TRUE_DGP_NEW_DISCRETE_WSQS: YES')
print('FITTED_NUISANCE_FULL_GRID: NO')
print('RECRUITMENT_CLUSTERING: NO')
print('HOUSEHOLD_SURVEY_DATA_READ: NO')
print('PUBLICATION_INFERENCE_CERTIFIED: NO')
print('EMPIRICAL_20_RESULT: COMPLETE')
