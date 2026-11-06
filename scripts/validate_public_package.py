#!/usr/bin/env python3
"""Fail-closed aggregate validation; no private data needed and no dependencies."""
from pathlib import Path
import csv,math,re,itertools,sys
from decimal import Decimal as D
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/aggregates'
NAMES=['Municipal piped','Groundwater','Surface water','Tanker delivery','Commercial water']
CNT=[1275,1212,234,120,225]

def fail(message):raise AssertionError(message)

def require(condition, message):
 if not condition: fail(message)
def read(fn):
 with (DATA/fn).open(newline='',encoding='utf-8') as f:
  reader=csv.DictReader(f)
  headers=reader.fieldnames
  if not headers or len(headers)!=len(set(headers)) or any(not h for h in headers):
   fail('Missing or duplicate CSV header in '+fn)
  out=list(reader)
 if not out:fail('Empty: '+fn)
 if any(None in row or any(v is None for v in row.values()) for row in out):
  fail('CSV field-count mismatch in '+fn)
 return out

def dec(x):return D(x.strip().replace('−','-').replace('+','').replace(',',''))
def close(a,b,tol,why):
 if abs(D(str(a))-D(str(b)))>D(str(tol)):fail(f'{why}: {a} vs {b}, tol {tol}')
def interval(x):
 m=re.fullmatch(r'\[\s*([−+\-\d.]+)\s*,\s*([−+\-\d.]+)\s*\]',x.strip())
 if not m:fail('Malformed interval '+x)
 l,u=dec(m.group(1)),dec(m.group(2))
 if l>u:fail('Interval inverted '+x)
 return l,u

means=read('main_table_1_observed_standardized_means.csv')
require((len(means) == 5), 'Scientific validation failed at check 29: len(means) == 5')
require(([r['Primary drinking-water source'] for r in means] == NAMES), "Scientific validation failed at check 30: [r['Primary drinking-water source'] for r in means] == NAMES")
require(([int(r['n'].replace(',', '')) for r in means] == CNT and sum(CNT) == 3066), "Scientific validation failed at check 31: [int(r['n'].replace(',', '')) for r in means] == CNT and sum(CNT) == 3066")
for row in means:
 obs,std,delta=[dec(row[k]) for k in ('Observed mean WSQS','Standardized mean WSQS','Observed − standardized')]
 close(obs-std,delta,'0.002',row['Primary drinking-water source']+' rounded adjustment')
 require((0 <= obs <= 100 and 0 <= std <= 100), 'Scientific validation failed at check 35: 0 <= obs <= 100 and 0 <= std <= 100')
require((dec(means[0]['Standardized mean WSQS']) - dec(means[1]['Standardized mean WSQS']) == D('-0.602')), "Scientific validation failed at check 36: dec(means[0]['Standardized mean WSQS']) - dec(means[1]['Standardized mean W")

threshold=read('main_table_2_standardized_threshold_probabilities.csv')
require((len(threshold) == 5), 'Scientific validation failed at check 39: len(threshold) == 5')
require(([r['Primary drinking-water source'] for r in threshold] == NAMES), "Scientific validation failed at check 40: [r['Primary drinking-water source'] for r in threshold] == NAMES")
T=[50,60,70]
for r in threshold:
 v=[dec(r[f'WSQS ≤{t} (%)']) for t in T]
 require((0 <= v[0] <= v[1] <= v[2] <= 100), 'Scientific validation failed at check 44: 0 <= v[0] <= v[1] <= v[2] <= 100')

cdf=read('full_41_threshold_standardized_cdf_percent.csv')
require((len(cdf) == 41), 'Scientific validation failed at check 47: len(cdf) == 41')
thresholds=[dec(r['WSQS_threshold']) for r in cdf]
require((thresholds == [D('2.5') * i for i in range(41)]), "Scientific validation failed at check 49: thresholds == [D('2.5') * i for i in range(41)]")
COLS=['CDF_PIPED_percent','CDF_GROUNDWATER_percent','CDF_SURFACE_percent','CDF_TANKER_percent','CDF_COMMERCIAL_percent']
curve={name:[dec(r[k]) for r in cdf] for name,k in zip(NAMES,COLS)}
for name,vs in curve.items():
 require((all((D('0') <= v <= D('100') for v in vs))), (name, 'CDF outside [0,100]'))
 require((all((a <= b for a, b in zip(vs, vs[1:])))), (name, 'nonmonotone CDF'))
 close(vs[-1],100,'0.000001',name+' CDF endpoint')
for i,t in enumerate(T):
 idx=thresholds.index(D(t))
 for j,name in enumerate(NAMES):
  close(curve[name][idx],threshold[j][f'WSQS ≤{t} (%)'],'0.000001',f'CDF vs main Table 2 {name} t={t}')
for idx,row in enumerate(cdf):
 diff=curve[NAMES[0]][idx]-curve[NAMES[1]][idx]
 close(diff,dec(row['PIPED_minus_GROUNDWATER_pp']),'0.0011',f'full CDF diff at {row["WSQS_threshold"]}')
maxi=max(range(len(cdf)),key=lambda i:abs(curve[NAMES[0]][i]-curve[NAMES[1]][i]))
require((thresholds[maxi] == D('70')), "Scientific validation failed at check 64: thresholds[maxi] == D('70')")
close(curve[NAMES[0]][maxi]-curve[NAMES[1]][maxi],D('3.765'),'0.002','largest piped-groundwater gap')

# Reproduce published distances by numerical integration from CDF grid; no source microdata.
expected={(1,0):D('.93'),(2,0):D('1.71'),(2,1):D('.92'),(3,0):D('3.95'),(3,1):D('4.62'),(3,2):D('5.12'),(4,0):D('3.17'),(4,1):D('3.78'),(4,2):D('4.64'),(4,3):D('4.02')}
for (i,j),target in expected.items():
 a,b=curve[NAMES[i]],curve[NAMES[j]]
 v=sum((abs(a[k]-b[k])+abs(a[k+1]-b[k+1]))*(thresholds[k+1]-thresholds[k])/D('200') for k in range(40))
 close(v,target,'0.006',f'pairwise integrated CDF distance {NAMES[i]} / {NAMES[j]}')

q11=read('si_table_s3_drinking_domestic_source_comparison.csv')
require(([r['Source'] for r in q11] == NAMES), "Scientific validation failed at check 75: [r['Source'] for r in q11] == NAMES")
require(([int(r['Drinking n'].replace(',', '')) for r in q11] == CNT), "Scientific validation failed at check 76: [int(r['Drinking n'].replace(',', '')) for r in q11] == CNT")
require((sum((int(r['Domestic n'].replace(',', '')) for r in q11)) == 3044), "Scientific validation failed at check 77: sum((int(r['Domestic n'].replace(',', '')) for r in q11)) == 3044")
same=[]
for r in q11:
 n=int(r['Same source n (%)'].split()[0].replace(',',''))
 same.append(n)
 close(D(n)*100/D(r['Drinking n'].replace(',','')),D(r['Same source n (%)'].split('(')[1].rstrip(')')),'0.051','Q1.1 agreement rounded pct')
require((sum(same) == 2012 and sum(CNT) - sum(same) == 1054), 'Scientific validation failed at check 83: sum(same) == 2012 and sum(CNT) - sum(same) == 1054')

sens=read('si_table_s4_domestic_source_sensitivity.csv')
require(([r['Source'] for r in sens] == NAMES), "Scientific validation failed at check 86: [r['Source'] for r in sens] == NAMES")
for j,r in enumerate(sens):
 close(dec(r['Drinking standardized WSQS']),dec(means[j]['Standardized mean WSQS']),'0.001','domestic sensitivity primary')
 close(dec(r['Domestic standardized WSQS'])-dec(r['Drinking standardized WSQS']),dec(r['Domestic − Drinking']),'0.002','domestic delta rounded')
require(([r['Domestic ESS/n (%)'] for r in sens][3:] == ['41.2', '65.3']), "Scientific validation failed at check 90: [r['Domestic ESS/n (%)'] for r in sens][3:] == ['41.2', '65.3']")

s5=read('si_table_s5_four_contrast_working_intervals.csv')
require((len(s5) == 4), 'Scientific validation failed at check 93: len(s5) == 4')
for r in s5:
 v=dec(r['Estimate'])
 for f in ['Marginal 95% interval','Four-contrast simultaneous 95% interval']:
  l,u=interval(r[f]);require((l <= v <= u), 'Scientific validation failed at check 97: l <= v <= u')
# Reconcile four selected contrast estimates to the primary mean and CDF outputs.
close(dec(s5[0]['Estimate']), dec(means[0]['Standardized mean WSQS'])-dec(means[1]['Standardized mean WSQS']), '0.0011', 'S5 mean contrast vs Table 1')
for offset,t in enumerate(T, start=1):
    gap=dec(threshold[0][f'WSQS ≤{t} (%)'])-dec(threshold[1][f'WSQS ≤{t} (%)'])
    close(dec(s5[offset]['Estimate']),gap,'0.0011',f'S5 contrast at t={t} vs Table 2')
s6=read('si_table_s6_ten_pairwise_mean_intervals.csv')
require((len(s6) == 10), 'Scientific validation failed at check 104: len(s6) == 10')
for r in s6:
 l,u=interval(r['Simultaneous 95% interval']);v=dec(r['Difference (WSQS points)']);require((l <= 0 <= u and l <= v <= u), 'Scientific validation failed at check 106: l <= 0 <= u and l <= v <= u')
# Check all ten pairwise point estimates against Table 1, not just interval containment.
for idx, (i,j) in enumerate(itertools.combinations(range(5),2)):
    target=dec(means[i]['Standardized mean WSQS'])-dec(means[j]['Standardized mean WSQS'])
    close(dec(s6[idx]['Difference (WSQS points)']),target,'0.0021',f'S6 pairwise {i}/{j} vs Table 1')

# Published-value transcription locks (copied from author-edited SI Tables S5/S6
# and independently compared externally). This does NOT recompute interval SEs.
EXPECTED_S5 = [('Mean WSQS (points)', '[−1.894, 0.691]', '[−2.199, 0.996]'), ('WSQS ≤50 (pp)', '[−2.483, 1.892]', '[−2.999, 2.408]'), ('WSQS ≤60 (pp)', '[−2.254, 4.245]', '[−3.022, 5.012]'), ('WSQS ≤70 (pp)', '[−0.191, 7.720]', '[−1.124, 8.654]')]
EXPECTED_S6 = [('Municipal piped − Groundwater', '[−3.008, 1.805]'), ('Municipal piped − Surface water', '[−5.683, 2.646]'), ('Municipal piped − Tanker delivery', '[−2.855, 9.489]'), ('Municipal piped − Commercial water', '[−2.711, 7.807]'), ('Groundwater − Surface water', '[−5.115, 3.282]'), ('Groundwater − Tanker delivery', '[−2.288, 10.125]'), ('Groundwater − Commercial water', '[−2.144, 8.443]'), ('Surface water − Tanker delivery', '[−3.129, 12.800]'), ('Surface water − Commercial water', '[−2.985, 11.118]'), ('Tanker delivery − Commercial water', '[−9.828, 8.290]')]
require([(r['Contrast'], r['Marginal 95% interval'],r['Four-contrast simultaneous 95% interval']) for r in s5] == EXPECTED_S5,
        'Reported SI Table S5 interval endpoints differ from the author-edited manuscript')
require([(r['Contrast (first source minus second)'],r['Simultaneous 95% interval']) for r in s6] == EXPECTED_S6,
        'Reported SI Table S6 interval endpoints differ from the author-edited manuscript')

s7=read('si_table_s7_core_domain_sensitivity.csv')
require((len(s7) == 5), 'Scientific validation failed at check 112: len(s7) == 5')
for j,r in enumerate(s7):close(dec(r['CORE-8']),dec(means[j]['Standardized mean WSQS']),'0.001','CORE8 standard')

support=read('overlap_diagnostics_original_aggregate.csv');require((len(support) == 5), 'Scientific validation failed at check 115: len(support) == 5')
for j,r in enumerate(support):
 require((int(r['n']) == CNT[j]), "Scientific validation failed at check 117: int(r['n']) == CNT[j]")
 close(dec(r['ess'])/D(CNT[j]),dec(r['ess_fraction']),'0.00011','ESS fraction')
 require((D(0) < dec(r['p01_target']) <= 1), "Scientific validation failed at check 119: D(0) < dec(r['p01_target']) <= 1")
 require((D(0) < dec(r['largest_normalized_weight_share']) <= 1), "Scientific validation failed at check 120: D(0) < dec(r['largest_normalized_weight_share']) <= 1")

eight=read('eight_covariate_sensitivity.csv');require((len(eight) == 5), 'Scientific validation failed at check 122: len(eight) == 5')
for j,r in enumerate(eight):
 close(dec(r['six_X_standardized_WSQS']),dec(means[j]['Standardized mean WSQS']),'0.001','eight X primary')
 close(dec(r['eight_X_standardized_WSQS'])-dec(r['six_X_standardized_WSQS']),dec(r['eight_X_minus_six_X_WSQS']),'0.0011','eight X change rounded')
sims=read('synthetic_simulation_coverage_summary.csv');require((len(sims) == 4), 'Scientific validation failed at check 126: len(sims) == 4')
require((sorted((round(float(r['full_grid_coverage_percent']), 1) for r in sims)) == [42.5, 57.9, 95.0, 95.8]), "Scientific validation failed at check 127: sorted((round(float(r['full_grid_coverage_percent']), 1) for r in sims)) ==")
sd=read('observed_wsqs_sd_by_source.csv');require((len(sd) == 5), 'Scientific validation failed at check 128: len(sd) == 5')
require(([r['source'] for r in sd] == NAMES), "Scientific validation failed at check 129: [r['source'] for r in sd] == NAMES")
require(([int(r['n']) for r in sd] == CNT), "Scientific validation failed at check 130: [int(r['n']) for r in sd] == CNT")
for r in sd:require((D('0') < dec(r['observed_WSQS_sample_SD']) < D('50')), "Scientific validation failed at check 131: D('0') < dec(r['observed_WSQS_sample_SD']) < D('50')")
close(min(dec(r['observed_WSQS_sample_SD']) for r in sd),D('14.0452275'),'0.0000001','minimum source SD')
close(max(dec(r['observed_WSQS_sample_SD']) for r in sd),D('17.9653094'),'0.0000001','maximum source SD')
smd=read('figure3_source_level_max_absolute_smd.csv');require((len(smd) == 5), 'Scientific validation failed at check 134: len(smd) == 5')
require(([r['source'] for r in smd] == NAMES and [int(r['n']) for r in smd] == CNT), "Scientific validation failed at check 135: [r['source'] for r in smd] == NAMES and [int(r['n']) for r in smd] == CNT")
for r in smd:
 require((D('0') <= dec(r['max_abs_SMD_after']) <= dec(r['max_abs_SMD_before']) < D('0.5')), "Scientific validation failed at check 137: D('0') <= dec(r['max_abs_SMD_after']) <= dec(r['max_abs_SMD_before']) < D('")
require(([r['max_abs_SMD_after'] for r in smd] == ['0.01062516', '0.01430716', '0.15649369', '0.18846583', '0.06418441']), "Scientific validation failed at check 138: [r['max_abs_SMD_after'] for r in smd] == ['0.01062516', '0.01430716', '0.15")
require(([r['max_abs_SMD_before'] for r in smd] == ['0.10740198', '0.11693133', '0.16863009', '0.20660413', '0.27832483']), "Scientific validation failed at check 139: [r['max_abs_SMD_before'] for r in smd] == ['0.10740198', '0.11693133', '0.1")

print('PASS — four selected and ten pairwise contrasts reconcile with main table point estimates')
print('PASS — five source SDs and before/after SMD maxima verified')
print('PASS — main and SI aggregate tables: group counts, mean adjustments, and 4/10 interval families')
print('PASS — full 41-threshold by 5-group CDF grid; three main thresholds; monotonicity and [0,100] endpoint')
print('PASS — all ten Figure S3 integrated absolute CDF distances from 41-threshold grid')
print('PASS — Q1.1 correspondence and mean sensitivity; six-X/eight-X and CORE7A/B aggregate consistency')
print('PASS — support summaries and independent-synthetic scenario coverage summaries')
