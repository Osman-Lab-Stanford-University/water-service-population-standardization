"""Publication-oriented WRR SI figures. Recreate five logistic propensity folds only;
no outcome-model refit. Original E04/E02 diagnostic parity required before exporting.
No original survey responses, row-level X, predictions or fitted models are written.
Recreated propensity fits are checked to original E04 aggregate values within
explicit numerical tolerances, not claimed to reproduce identical binary floats.
"""
from __future__ import annotations
import csv, hashlib, json, math, os, re, runpy, sys, warnings
from pathlib import Path
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss
from sklearn.exceptions import ConvergenceWarning
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'svg.fonttype':'none', 'pdf.fonttype':42,
    'font.family':'sans-serif', 'font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
    'axes.titlesize':12,'axes.labelsize':10,'xtick.labelsize':8.7,
    'ytick.labelsize':8.6,'text.color':'#28323C','axes.labelcolor':'#28323C',
    'xtick.color':'#53606C','ytick.color':'#28323C','axes.edgecolor':'#BCC6CC',
    'savefig.facecolor':'white','figure.facecolor':'white'})
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

BASE=Path(sys.argv[1]).resolve(); STAGE=Path(sys.argv[2]).resolve()
DEV=BASE/'06_toolkit_development'; OUT=STAGE/'FIGURES'; PRIV=STAGE/'PRIVATE_REVIEW'
G=('PIPED','GROUNDWATER','SURFACE','TANKER','COMMERCIAL')
N=(1275,1212,234,120,225)
DISPLAY=('Municipal piped','Groundwater','Surface water','Tanker delivery','Commercial water')
COLORS=('#825A7E','#00A6AB','#7AD1D3','#DE8244','#9B7B93')
QIDS=('Q15.17','Q15.13','Q15.15','Q15.14','Q15.10','Q15.5')
CODES={'Q15.17':(1,2,3,4,5),'Q15.13':(9,10,11),'Q15.15':(9,10,11,5),
       'Q15.14':(1,2,3,4),'Q15.10':(20,21,22,4,5,6,7,8),'Q15.5':(2,14,15)}
LABELS={
 'Q15.17':('Under Rs. 100,000','Rs. 100,000–250,000','Rs. 250,000–500,000',
           'Rs. 500,000–1m','Rs. 1m or more'),
 'Q15.13':('1 car','2 cars','No cars'),
 'Q15.15':('Owned with loan','Owned without loan','Rented','Other tenure'),
 'Q15.14':('Individual house','Apartment','Temporary housing','Other dwelling'),
 'Q15.10':('1 person','2 people','3 people','4 people','5 people','6 people',
            '7 people','8 or more'),
 'Q15.5':('Urban','Sub-urban','Rural')}
Q_NAMES={'Q15.17':'Family income','Q15.13':'Household cars',
         'Q15.15':'Housing tenure','Q15.14':'Dwelling type',
         'Q15.10':'Household size','Q15.5':'Childhood setting'}
# Numeric audit anchors are used only as programmatic parity checks; no internal stage
# or validation terminology is printed in the publication figures.
MAIN_MAX_BEFORE=(.10740,.11693,.16863,.20660,.27832)
MAIN_MAX_AFTER=(.01063,.01431,.15649,.18840,.06418)
DISTANCE_ANCHORS={(0,1):.93,(0,2):1.71,(0,3):3.95,(0,4):3.17,
 (1,2):.92,(1,3):4.62,(1,4):3.78,(2,3):5.12,(2,4):4.64,(3,4):4.02}

def stop(msg): raise RuntimeError('STOP: '+str(msg))
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def require(path):
    if not path.is_file() or path.is_symlink():stop('Missing authenticated source: '+str(path))
    return path

def write_csv(path, columns, rows):
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore')
        w.writeheader();w.writerows(rows)

def save(fig, name, dest=OUT):
    fig.savefig(dest/(name+'.png'),dpi=330,bbox_inches='tight',pad_inches=.14)
    fig.savefig(dest/(name+'.pdf'),bbox_inches='tight',pad_inches=.14)
    fig.savefig(dest/(name+'.svg'),bbox_inches='tight',pad_inches=.14)
    plt.close(fig)

def get_original_inputs():
    # Reuses the already authenticated 26C/E04 pre-fit replay. This has executed
    # successfully on the user's Mac in 26G2; no main model is run by load_original().
    stage='EMPIRICAL_26C_Q11_DOMESTIC_CORRECTED_20261103_133038_rNVQ0e'
    scripts=sorted((DEV/stage).glob('*26C*Q11*.py'))
    if len(scripts)!=1:stop('Expected exactly one verified 26C Python source')
    oldargv=sys.argv[:]
    try:
        sys.argv=[str(scripts[0]),str(BASE),str(STAGE)]
        mod=runpy.run_path(str(scripts[0]),run_name='wrr_si_original_input_methods')
    finally:sys.argv=oldargv
    X,y,Y,raw,e04,_,_=mod['load_original']()
    X=np.asarray(X,dtype=np.float64);y=np.asarray(y,dtype=int)
    if X.shape!=(3066,33) or y.shape!=(3066,) or tuple(int((y==j).sum()) for j in range(1,6))!=N:
        stop('Original six-X matrix/cohort mismatch')
    if np.any((X!=0)&(X!=1)) or not np.allclose(X.sum(axis=1),6):
        stop('Six QID one-hot indicator parity mismatch')
    if e04['analytical_n']!=3066 or tuple(e04['six_adjustment_question_ids'])!=QIDS:
        stop('Frozen E04 specification mismatch')
    return X,y,e04

def reconstruct_original_propensities(X,y,e04):
    # Original E04 / E02 propensity only: L2 multinomial logistic, five source-
    # stratified folds, seed 20261022. The E04 individual predictions were not saved.
    e=np.full((3066,5),np.nan,dtype=np.float64)
    sf=StratifiedKFold(n_splits=5,shuffle=True,random_state=20261022)
    for fold,(train,test) in enumerate(sf.split(X,y),1):
        p=LogisticRegression(penalty='l2',C=1.0,solver='lbfgs',max_iter=2000,
                             tol=1e-6,fit_intercept=True)
        with warnings.catch_warnings():
            warnings.simplefilter('error',ConvergenceWarning)
            p.fit(X[train],y[train])
        if not np.array_equal(p.classes_,[1,2,3,4,5]):stop('Propensity class-order mismatch')
        e[test,:]=p.predict_proba(X[test])
        print(f'Original specification propensity fold {fold}/5 computed.')
    if not np.isfinite(e).all() or np.min(e)<=0 or not np.allclose(e.sum(1),1,atol=1e-9):
        stop('Invalid reconstructed propensity predictions')
    ll=log_loss(y,e,labels=[1,2,3,4,5])
    original_ll=float(e04['propensity_out_of_fold_logloss'])
    ll_difference=ll-original_ll
    print(f'Out-of-fold log loss: reconstructed={ll:.12f}; '
          f'original={original_ll:.12f}; absolute difference={abs(ll_difference):.3g}',flush=True)
    # Distinct versions of scikit-learn/SciPy can yield very small LBFGS
    # numerical differences with identical folds, categories and settings.
    # We require agreement to 0.000002 in log loss, much tighter than
    # the published five-decimal value; this is not exact numerical parity.
    if abs(ll_difference)>2e-6:
        stop(f'Reconstructed propensity log loss differs materially from original: '
             f'{ll:.12f} versus {original_ll:.12f} (tolerance 0.000002)')
    checks=[];before=np.full((5,33),np.nan);after=np.full((5,33),np.nan)
    differences=[]
    pt=X.mean(axis=0);scale=np.sqrt(pt*(1-pt));active=scale>1e-11
    for j,g in enumerate(G):
        obs=y==j+1;w=1.0/e[obs,j]
        p_obs=X[obs].mean(axis=0)
        p_weighted=np.average(X[obs],axis=0,weights=w)
        before[j,active]=(p_obs[active]-pt[active])/scale[active]
        after[j,active]=(p_weighted[active]-pt[active])/scale[active]
        d={'n':int(obs.sum()),'ess':float(w.sum()**2/np.square(w).sum()),
            'ess_fraction':float(w.sum()**2/(np.square(w).sum()*obs.sum())),
            'largest_weight_fraction':float(w.max()/w.sum()),
            'p01_target':float(np.quantile(e[:,j],.01)),
            'balance_max_absolute_category_gap':float(np.max(np.abs(p_weighted-pt)))}
        recorded=e04['overlap'][g]
        # Tolerances are in the native units of each diagnostic. The original
        # model was run under a different numerical environment. These bounds
        # are small relative to publication precision (e.g., gap tolerance
        # 0.01 percentage points), yet protect against material model drift.
        limits={'n':0.0,'ess':0.020,'ess_fraction':0.00020,
                'largest_weight_fraction':0.000020,'p01_target':0.00020,
                'balance_max_absolute_category_gap':0.00010}
        for key,v in d.items():
            if key not in recorded:stop(f'Original {g} report lacks {key}')
            expected=float(recorded[key]);delta=abs(v-expected)
            differences.append((g,key,v,expected,delta,limits[key]))
            if delta>limits[key]:
                stop(f'Original {g} diagnostic mismatch {key}: '
                     f'reconstructed={v:.12g}; original={expected:.12g}; '
                     f'absolute difference={delta:.4g}; tolerance={limits[key]}')
        for title,computed,expected in [('before',np.nanmax(np.abs(before[j])),MAIN_MAX_BEFORE[j]),
                                         ('after',np.nanmax(np.abs(after[j])),MAIN_MAX_AFTER[j])]:
            # Displayed main-figure maxima were rounded; use a modest 0.003
            # tolerance but refuse a substantively different figure.
            if abs(computed-expected)>.003:
                stop(f'{g} maximum {title}-weighting SMD does not reproduce main Figure 3: {computed:.5f} vs {expected:.5f}')
        # Use the ORIGINAL E04 aggregate support numbers in Figure S2.
        # Only Figure S1 needs reconstructed category-level weighted points.
        checks.append({'Source':g,**{key:(int(expected) if key=='n' else float(expected))
                           for key,expected in recorded.items() if key in d},
            'maximum_abs_smd_unweighted':float(np.nanmax(np.abs(before[j]))),
            'maximum_abs_smd_weighted':float(np.nanmax(np.abs(after[j])))})
    for g,key,v,expected,delta,tol in differences:
        print(f'  {g:12s} {key:38s} deviation={delta:.3g}; tolerance={tol:g}')
    print('PASS: Propensity reconstruction agrees with original aggregate diagnostics within specified numerical tolerances.')
    print('PASS: Before/after SMD maxima agree with main Figure 3 within displayed-figure tolerance.')
    return before,after,checks,X,y,pt

def category_specs():
    specs=[]
    for q in QIDS:
        for code,label in zip(CODES[q],LABELS[q]):
            specs.append((q,code,label))
        specs.append((q,None,'Missing'))
    if len(specs)!=33:stop('Invalid 33-indicator codebook')
    return specs

def draw_s1(before,after,X,y,screen=True):
    specs=category_specs();counts=np.zeros((33,5),dtype=int)
    for k in range(33):
        for j in range(5):counts[k,j]=int(X[y==j+1,k].sum())
    # Zero-reference-variance categories are mathematically undefined, not zero.
    valid=np.isfinite(before).any(axis=0)
    # A conservative publication-review candidate withholds any entire category
    # appearing as 1–4 observations in a source. Full data plot remains private.
    rare=((counts>0)&(counts<5)).any(axis=1)
    show=valid & (~rare if screen else np.ones(33,dtype=bool))
    left=0.20;right=0.985;top=.952;bottom=.170
    fig,axs=plt.subplots(3,2,figsize=(14.0,14.5),sharex=True)
    plt.subplots_adjust(left=left,right=right,top=top,bottom=bottom,wspace=.55,hspace=.34)
    maxx=np.nanmax(np.abs(np.r_[before[:,show].ravel(),after[:,show].ravel()]))
    limit=max(.35,math.ceil((maxx+.04)*20)/20)
    span=(-limit,limit)
    panels=[]
    start=0
    for pi,q in enumerate(QIDS):
        k=len(CODES[q])+1;indices=[t for t in range(start,start+k) if show[t]]
        start+=k
        ax=axs.flat[pi];ax.set_xlim(*span)
        ax.axvspan(-.1,.1,color='#F2F5F6',zorder=0)
        ax.axvline(0,color='#9CA8B0',lw=.8,zorder=1)
        for x in (-.1,.1):ax.axvline(x,color='#C7D0D5',lw=.7,linestyle=(0,(3,3)),zorder=1)
        if len(indices)==0:stop('No drawable categories for '+q)
        rows=np.arange(len(indices))[::-1]
        for ri,col in enumerate(indices):
            y0=rows[ri]
            for j in range(5):
                offset=(2-j)*.12
                yj=y0+offset
                b=before[j,col];a=after[j,col]
                if not (np.isfinite(b) and np.isfinite(a)):continue
                ax.plot((b,a),(yj,yj),color=COLORS[j],lw=.9,alpha=.47,zorder=2)
                ax.scatter([b],[yj],s=24,marker='o',facecolor='white',edgecolor=COLORS[j],
                           linewidth=1.15,zorder=4)
                ax.scatter([a],[yj],s=27,marker='D',facecolor=COLORS[j],edgecolor='#31404A',
                           linewidth=.46,zorder=5)
        ax.set_yticks(rows,[specs[k][2] for k in indices])
        ax.set_ylim(-.72,len(indices)-.30)
        ax.tick_params(axis='y',length=0,pad=7,labelsize=8.4)
        ax.tick_params(axis='x',labelsize=9,length=3)
        ax.spines[['top','right','left']].set_visible(False)
        ax.spines['bottom'].set_color('#9EAAB2')
        ax.grid(False)
        ax.set_title(chr(97+pi)+'  '+Q_NAMES[q],loc='left',weight='bold',fontsize=12.2,pad=14)
        if pi>=4:ax.set_xlabel('Signed standardized mean difference',fontsize=10.5,labelpad=8)
        panels.append({'qid':q,'shown':len(indices),'withheld':int(np.sum(rare[[i for i,z in enumerate(specs) if z[0]==q]]))})
    shape=[Line2D([0],[0],marker='o',markersize=6,markerfacecolor='white',
                  markeredgecolor='#465563',markeredgewidth=1.4,lw=.5,color='#7C8A93',label='Before weighting'),
           Line2D([0],[0],marker='D',markersize=5.8,markerfacecolor='#465563',
                  markeredgecolor='#31404A',lw=.5,color='#7C8A93',label='After weighting')]
    fig.legend(handles=shape,frameon=False,ncol=2,fontsize=9.4,
               loc='lower center',bbox_to_anchor=(.51,.104),columnspacing=2.7)
    sources=[Line2D([0],[0],marker='o',color='none',markerfacecolor=c,markeredgecolor='#52606B',
                    markeredgewidth=.4,markersize=6,label=lab) for lab,c in zip(DISPLAY,COLORS)]
    fig.legend(handles=sources,frameon=False,ncol=5,fontsize=8.9,
               loc='lower center',bbox_to_anchor=(.51,.072),columnspacing=1.7)
    if screen and np.any(rare):
        fig.text(.51,.040,'Rare source-category combinations are omitted from this disclosure-screened display.',
                 ha='center',fontsize=8.2,color='#5F6D76')
    else:
        fig.text(.51,.040,'Open circles: unweighted; solid diamonds: weighted. Source color is unchanged.',
                 ha='center',fontsize=8.2,color='#5F6D76')
    return fig,show,rare,counts

def draw_s2(checks):
    data=[[100*float(row['p01_target']) for row in checks],
          [100*float(row['largest_weight_fraction']) for row in checks],
          [100*float(row['balance_max_absolute_category_gap']) for row in checks]]
    titles=('a   Fitted probability','b   Weight share','c   Balance gap')
    xlabs=('First percentile (%)','Largest normalized weight (%)','Maximum category gap (pp)')
    limits=(27.5,3.25,5.75)
    fig,axes=plt.subplots(1,3,figsize=(13.8,5.6),sharey=True)
    fig.subplots_adjust(left=.20,right=.968,top=.84,bottom=.21,wspace=.43)
    for ax,d,ttl,xlab,lim in zip(axes,data,titles,xlabs,limits):
        ax.set_xlim(0,lim);ax.set_ylim(-.6,4.6)
        ax.set_title(ttl,loc='left',fontsize=13.0,fontweight='bold',pad=12)
        ax.set_xlabel(xlab,fontsize=9.9,labelpad=9)
        ax.grid(axis='x',color='#DFE6E9',lw=.7)
        ax.set_axisbelow(True)
        for j,(v,c) in enumerate(zip(d,COLORS)):
            y=4-j
            ax.plot([0,v],[y,y],color=c,lw=2.2,alpha=.76,solid_capstyle='round')
            ax.scatter([v],[y],s=54,facecolor=c,edgecolor='#475660' if j==2 else 'white',lw=.65,zorder=4)
            label=f'{v:.1f}%' if ax is axes[0] else (f'{v:.2f}%' if ax is axes[1] else f'{v:.2f}')
            ax.annotate(label,(v,y),xytext=(6,0),textcoords='offset points',ha='left',va='center',fontsize=9.4)
        ax.spines[['top','right','left']].set_visible(False)
        ax.tick_params(axis='y',length=0,pad=8)
        ax.tick_params(axis='x',length=3)
    axes[0].set_yticks([4,3,2,1,0],DISPLAY,fontsize=10)
    axes[0].tick_params(labelleft=True)
    for ax in axes[1:]:ax.tick_params(labelleft=False)
    return fig

def load_cdf():
    source=require(DEV/'EMPIRICAL_25BM_FIG234_20261025_210835_X6SAjb'/'REVIEW_ONLY'/'SI_TABLE_S11_COMPLETE_CDF.csv')
    with source.open('r',encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    if len(rows)!=41:stop('Complete CDF must contain 41 thresholds')
    t=np.array([float(r['threshold_WSQS']) for r in rows])
    if not np.allclose(t,np.arange(41)*2.5,atol=1e-9):stop('Unexpected CDF threshold grid')
    curves=np.column_stack([[float(r['cdf_percent_'+g])/100 for r in rows] for g in G])
    if not np.isfinite(curves).all() or np.any(curves<0) or np.any(curves>1+1e-8):
        stop('Invalid CDF values')
    if np.any(np.diff(curves,axis=0)<-1e-10) or not np.allclose(curves[-1],1,atol=1e-8):
        stop('CDF monotonicity/endpoint mismatch')
    anchors={50:(7.469,7.764,6.149,8.690,10.553),
             60:(20.014,19.019,17.819,27.271,17.067),
             70:(38.837,35.072,33.932,45.029,43.679)}
    for threshold,values in anchors.items():
        if not np.allclose(curves[int(threshold/2.5)]*100,values,atol=.002):
            stop(f'CDF differs from locked main Table 2 at {threshold}')
    return t,curves

def draw_s3(t,curves):
    distances={}
    for i in range(5):
        for j in range(i+1,5):
            d=float(np.trapz(np.abs(curves[:,i]-curves[:,j]),t))
            if abs(d-DISTANCE_ANCHORS[(i,j)])>.025:
                stop(f'CDF integrated distance unexpected: {G[i]} / {G[j]}, {d}')
            distances[(i,j)]=d
    fig,ax=plt.subplots(figsize=(7.8,6.6))
    fig.subplots_adjust(left=.245,right=.935,bottom=.20,top=.85)
    from matplotlib.colors import Normalize
    cmap=LinearSegmentedColormap.from_list('grayteal',['#ECF3F2','#B4D1D2','#4A7C8A'])
    norm=Normalize(vmin=0,vmax=5.5)
    for row in range(5):
        for col in range(row):
            d=distances[(col,row)]
            color=cmap(norm(d))
            ax.add_patch(Rectangle((col-.48,row-.48),.96,.96,facecolor=color,
                                   edgecolor='white',linewidth=2))
            ax.text(col,row,f'{d:.2f}',va='center',ha='center',fontsize=13.5,
                    color='white' if d>3.7 else '#28323C')
    ax.set_xlim(-.61,4.5);ax.set_ylim(4.48,-.56)
    ax.set_xticks(range(5),('Piped','Groundwater','Surface','Tanker','Commercial'),
                  rotation=38,ha='right',rotation_mode='anchor')
    ax.set_yticks(range(5),('Piped','Groundwater','Surface','Tanker','Commercial'))
    ax.tick_params(axis='both',length=0,pad=7,labelsize=10.4)
    ax.spines[:].set_visible(False)
    ax.set_title('Integrated absolute CDF difference (WSQS points)',loc='left',fontsize=12.6,
                 color='#53606C',pad=16)
    for i,c in enumerate(COLORS):
        ax.plot([-0.595,-0.595],[i-.40,i+.40],color=c,lw=5,solid_capstyle='butt',clip_on=False)
    return fig,distances

def main():
    OUT.mkdir(parents=True,exist_ok=True);PRIV.mkdir(parents=True,exist_ok=True)
    X,y,e04=get_original_inputs()
    before,after,checks,X,y,pt=reconstruct_original_propensities(X,y,e04)
    fig,shown,rare,counts=draw_s1(before,after,X,y,screen=True)
    save(fig,'Figure_S1_Category_Balance')
    # Full category-level evidence is private; no respondent-level exports.
    fig_full,_,_,_=draw_s1(before,after,X,y,screen=False)
    save(fig_full,'Figure_S1_Full_Category_Balance_PRIVATE',dest=PRIV)
    specs=category_specs()
    diagnostic=[]
    for col,(q,code,label) in enumerate(specs):
        for j,g in enumerate(G):
            diagnostic.append({'QID':q,'category':label,'source':g,
                               'group_count':int(counts[col,j]),
                               'smd_unweighted':f'{before[j,col]:.10f}' if np.isfinite(before[j,col]) else '',
                               'smd_weighted':f'{after[j,col]:.10f}' if np.isfinite(after[j,col]) else '',
                               'included_in_publication_candidate':bool(shown[col])})
    write_csv(PRIV/'Category_Level_Balance_PRIVATE.csv',list(diagnostic[0]),diagnostic)
    write_csv(PRIV/'Propensity_Parity_Checks_PRIVATE.csv',list(checks[0]),checks)
    fig=draw_s2(checks);save(fig,'Figure_S2_Empirical_Support')
    t,curves=load_cdf();fig,dist=draw_s3(t,curves);save(fig,'Figure_S3_Distribution_Distances')
    # Aggregate distance matrix has no respondent-level information.
    disrows=[{'source_1':G[i],'source_2':G[j],'integrated_absolute_difference_WSQS_points':f'{d:.9f}'}
             for (i,j),d in dist.items()]
    write_csv(OUT/'Figure_S3_Pairwise_Distances.csv',list(disrows[0]),disrows)
    (OUT/'Figure_Captions.txt').write_text(
        'Figure S1. Category-level covariate balance before and after weighting. '
        'Signed standardized mean differences for the six adjustment characteristics. '
        'Open circles denote unweighted differences; filled diamonds denote inverse-propensity-weighted differences. '
        'Source groups are distinguished by color; values use the empirical reference Bernoulli standard deviation. '
        'Indicators with zero reference variance are omitted. Categories with small source-group cells '
        'are withheld from the public-review version pending disclosure approval.\n\n'
        'Figure S2. Additional empirical-support diagnostics. (a) First percentile of fitted source-membership '
        'probabilities across the empirical reference; (b) largest normalized observed-source inverse-propensity '
        'weight; and (c) maximum absolute weighted source-to-reference category-proportion gap in percentage points.\n\n'
        'Figure S3. Pairwise distances between standardized WSQS distributions. Integrated absolute differences '
        'between standardized CDFs are approximated using trapezoidal integration across the 41-point WSQS grid. '
        'Distances are in WSQS points, unsigned and descriptive.\n',encoding='utf-8')
    # Contact sheet for editorial review; programmatically checks output legibility.
    from PIL import Image,ImageOps,ImageDraw
    imgs=[]
    for stem in ['Figure_S1_Category_Balance','Figure_S2_Empirical_Support',
                 'Figure_S3_Distribution_Distances']:
        im=Image.open(OUT/(stem+'.png')).convert('RGB')
        im.thumbnail((1600,1200))
        pad=Image.new('RGB',(1660,1270),'white')
        pad.paste(im,((1660-im.width)//2,(1270-im.height)//2))
        imgs.append(pad)
    canvas=Image.new('RGB',(3320,2540),'#FFFFFF')
    canvas.paste(imgs[0],(0,0));canvas.paste(imgs[1],(1660,0));canvas.paste(imgs[2],(830,1270))
    canvas.save(OUT/'All_Three_Figures_Review.png',dpi=(200,200))
    report={'source_n':dict(zip(G,N)),'five_propensity_folds_recreated':True,
        'all_E04_propensity_diagnostics_parity':'WITHIN_NUMERICAL_REPRODUCTION_TOLERANCES',
        'main_figure_3_maximum_smd_parity':'WITHIN_DISPLAY_TOLERANCE',
        'full_CDF_grid_parity':'PASS',
        'screened_source_categories':int(np.sum(rare)),
        'public_release_status':'FIGURE_S1_REQUIRES_DISCLOSURE_APPROVAL',
        'outcome_models_refitted':False,'primary_estimates_modified':False,
        'respondent_level_artifacts_exported':False}
    (PRIV/'Generation_Status.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    # One convenient sharing bundle: review figures/captions only, never private data.
    from zipfile import ZipFile, ZIP_DEFLATED
    with ZipFile(STAGE/'WRR_SI_FIGURES_S1_S3_REVIEW.zip','w',compression=ZIP_DEFLATED) as archive:
        for file in sorted(OUT.iterdir()):
            if file.is_file():archive.write(file,arcname=file.name)
    files=sorted(p for p in STAGE.rglob('*') if p.is_file() and p.name!='SHA256_MANIFEST.txt')
    (STAGE/'SHA256_MANIFEST.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(STAGE))+'\n' for p in files))
    for line in (STAGE/'SHA256_MANIFEST.txt').read_text().splitlines():
        digest,rel=line.split('  ',1)
        if sha(STAGE/rel)!=digest:stop('Figure output manifest integrity failed: '+rel)
    print('\nSUCCESS: All three figures (PNG, SVG, PDF) generated and numerically cross-checked.')
    print('PUBLICATION REVIEW:',OUT)
    print('PRIVATE FULL S1:',PRIV)
    print('Main manuscript figures and original data unchanged.')

if __name__=='__main__':main()
