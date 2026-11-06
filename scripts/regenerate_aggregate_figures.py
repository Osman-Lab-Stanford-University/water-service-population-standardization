#!/usr/bin/env python3
"""Reconstruct comparable summary plots from released aggregate data (not locked journal figures)."""
from pathlib import Path
import csv,argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/aggregates'
COLOR=['#825A7E','#00A6AB','#7AD1D3','#DE8244','#9B7B93']
NAMES=['Municipal piped','Groundwater','Surface water','Tanker delivery','Commercial water']

def read(name):
 with (DATA/name).open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--out',default='generated_figures');args=ap.parse_args()
 out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
 m=read('main_table_1_observed_standardized_means.csv');ix=np.arange(5)
 observed=np.array([float(r['Observed mean WSQS']) for r in m]);std=np.array([float(r['Standardized mean WSQS']) for r in m])
 fig,ax=plt.subplots(figsize=(8.0,4.2))
 for j in range(5):
  ax.plot([observed[j],std[j]],[j,j],lw=2,color=COLOR[j],alpha=.65)
  ax.scatter(observed[j],j,facecolor='white',edgecolor=COLOR[j],s=60,zorder=3)
  ax.scatter(std[j],j,marker='D',color=COLOR[j],s=36,zorder=3)
 ax.set_yticks(ix,NAMES);ax.invert_yaxis();ax.set_xlabel('WSQS (0–100)');ax.set_title('Observed circles and standardized diamonds (aggregate reconstruction)');ax.grid(axis='x',alpha=.18)
 fig.tight_layout();fig.savefig(out/'reconstructed_mean_comparison.png',dpi=180);plt.close(fig)

 c=read('full_41_threshold_standardized_cdf_percent.csv')
 x=np.array([float(r['WSQS_threshold']) for r in c]);a=np.array([[float(r[k]) for k in ('CDF_PIPED_percent','CDF_GROUNDWATER_percent','CDF_SURFACE_percent','CDF_TANKER_percent','CDF_COMMERCIAL_percent')] for r in c]);
 fig,ax=plt.subplots(figsize=(8,4.7))
 for j in range(5):ax.plot(x,a[:,j],label=NAMES[j],color=COLOR[j],lw=1.8)
 ax.set(xlim=(0,100),ylim=(0,100),xlabel='WSQS threshold',ylabel='Standardized cumulative share (%)',title='41-threshold standardized CDFs (aggregate reconstruction)')
 ax.legend(ncol=3,loc='upper left',fontsize=8);ax.grid(alpha=.18)
 fig.tight_layout();fig.savefig(out/'reconstructed_full_standardized_CDFs.png',dpi=180);plt.close(fig)

 fig,ax=plt.subplots(figsize=(8,4.7))
 for j in range(1,5):ax.plot(x,a[:,0]-a[:,j],label='Piped minus '+NAMES[j],color=COLOR[j],lw=1.8)
 ax.axhline(0,color='black',lw=.65);ax.set(xlim=(0,100),xlabel='WSQS threshold',ylabel='Cumulative difference (percentage points)',title='Piped minus comparator (aggregate reconstruction)');ax.grid(alpha=.18);ax.legend(fontsize=8)
 fig.tight_layout();fig.savefig(out/'reconstructed_piped_CDF_differences.png',dpi=180);plt.close(fig)

 dist=np.zeros((5,5))
 for i in range(5):
  for j in range(5):dist[i,j]=np.sum((np.abs(a[:-1,i]-a[:-1,j])+np.abs(a[1:,i]-a[1:,j]))*np.diff(x)/200.0)
 fig,ax=plt.subplots(figsize=(6.4,5.2))
 im=ax.imshow(dist,vmin=0,vmax=6,cmap='Greys');ax.set_xticks(range(5),NAMES,rotation=35,ha='right');ax.set_yticks(range(5),NAMES)
 for i in range(5):
  for j in range(5):ax.text(j,i,f'{dist[i,j]:.2f}',ha='center',va='center',fontsize=9)
 ax.set_title('Integrated absolute CDF differences (WSQS points)');fig.colorbar(im,ax=ax,fraction=.05)
 fig.tight_layout();fig.savefig(out/'reconstructed_pairwise_CDF_distances.png',dpi=180);plt.close(fig)
 print('Aggregate-only plot reconstruction completed:',out)
