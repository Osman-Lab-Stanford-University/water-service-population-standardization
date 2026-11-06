#!/usr/bin/env python3
"""Independent, aggregate-only recreation of the ten Figure S3 CDF distances.

Outputs numerical integration of the *released finite-precision* 41-threshold
curves. These are not original person-level model fits or original journal art.
"""
from pathlib import Path
import argparse,csv,itertools

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/aggregates/full_41_threshold_standardized_cdf_percent.csv'
LABELS=['Municipal piped','Groundwater','Surface water','Tanker delivery','Commercial water']
COLS=['CDF_PIPED_percent','CDF_GROUNDWATER_percent','CDF_SURFACE_percent','CDF_TANKER_percent','CDF_COMMERCIAL_percent']

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',default='generated_figures',help='Output folder for reconstructed results')
    args=ap.parse_args()
    with SOURCE.open(newline='',encoding='utf-8') as f: rows=list(csv.DictReader(f))
    x=[float(r['WSQS_threshold']) for r in rows]
    curves=[[float(r[k]) for r in rows] for k in COLS]
    if len(x)!=41 or any(x[i+1]<=x[i] for i in range(40)):
        raise SystemExit('Invalid released threshold grid')
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    p=out/'reconstructed_figure_s3_pairwise_distances.csv'
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f)
        w.writerow(['source_1','source_2','integrated_absolute_CDF_distance_WSQS_points','calculation_scope'])
        for i,j in itertools.combinations(range(5),2):
            diff=[abs(curves[i][k]-curves[j][k])/100.0 for k in range(41)]
            value=sum((diff[k]+diff[k+1])*(x[k+1]-x[k])*0.5 for k in range(40))
            w.writerow([LABELS[i],LABELS[j],f'{value:.9f}','trapezoid on released rounded CDF grid'])
    print('Wrote 10 pairwise distances:',p)

if __name__=='__main__': main()
