#!/usr/bin/env python3
"""Independent read-only, standard-library regression tests of every released CSV."""
from pathlib import Path
from decimal import Decimal as D
import csv
import itertools
import json
import math
import re
import unittest

ROOT=Path(__file__).resolve().parents[1]
AGG=ROOT/'data'/'aggregates'
SOURCES=['Municipal piped','Groundwater','Surface water','Tanker delivery','Commercial water']
COUNTS=[1275,1212,234,120,225]
FIELDS={
 'eight_covariate_sensitivity.csv':('primary_drinking_source,six_X_standardized_WSQS,eight_X_standardized_WSQS,eight_X_minus_six_X_WSQS',5),
 'figure3_source_level_max_absolute_smd.csv':('source,n,max_abs_SMD_before,max_abs_SMD_after',5),
 'full_41_threshold_standardized_cdf_percent.csv':('WSQS_threshold,CDF_PIPED_percent,CDF_GROUNDWATER_percent,CDF_SURFACE_percent,CDF_TANKER_percent,CDF_COMMERCIAL_percent,PIPED_minus_GROUNDWATER_pp',41),
 'main_table_1_observed_standardized_means.csv':('Primary drinking-water source,n,Observed mean WSQS,Standardized mean WSQS,Observed − standardized',5),
 'main_table_2_standardized_threshold_probabilities.csv':('Primary drinking-water source,WSQS ≤50 (%),WSQS ≤60 (%),WSQS ≤70 (%)',5),
 'observed_wsqs_sd_by_source.csv':('source,n,observed_WSQS_sample_SD',5),
 'overlap_diagnostics_original_aggregate.csv':('source,n,ess,ess_fraction,largest_normalized_weight_share,max_weighted_category_prop_gap_pp,p01_target,E04_screen_pass',5),
 'si_table_s2_questionnaire_item_mapping.csv':('QID,Measure,English questionnaire wording,Response options',8),
 'si_table_s3_drinking_domestic_source_comparison.csv':('Source,Drinking n,Domestic n,Same source n (%),Drinking observed WSQS,Domestic observed WSQS',5),
 'si_table_s4_domestic_source_sensitivity.csv':('Source,Drinking standardized WSQS,Domestic standardized WSQS,Domestic − Drinking,Domestic ESS/n (%),Maximum weight share (%),Maximum balance gap (pp)',5),
 'si_table_s5_four_contrast_working_intervals.csv':('Contrast,Estimate,Marginal 95% interval,Four-contrast simultaneous 95% interval',4),
 'si_table_s6_ten_pairwise_mean_intervals.csv':('Contrast (first source minus second),Difference (WSQS points),Simultaneous 95% interval',10),
 'si_table_s7_core_domain_sensitivity.csv':('Source,CORE-8,CORE7A,CORE7B',5),
 'synthetic_simulation_coverage_summary.csv':('simulation_source,scenario,method,full_grid_coverage_percent,nominal_percent',4),
 'si_table_s8_outcome_model_sensitivity.csv':('Metric,n,HGBR_original,Random_Forest,Ridge,CV_selected_HGBR,Unit',8),
 'outcome_model_oof_rmse.csv':('Source,n,HGBR_original_RMSE,RF_RMSE,Ridge_RMSE,CV_selected_HGBR_RMSE',6),
 'outcome_model_paired_rmse_bootstrap.csv':('Source,n,Alternative_learner,RMSE_difference_vs_HGBR,Bootstrap_lower_2p5,Bootstrap_upper_97p5,Bootstrap_replicates',18),
 'outcome_model_hgbr_tuning_counts.csv':('Source,Outer_fold_count,Original,Conservative,Flexible,Slow_growth',5),
 'outcome_model_working_if_means.csv':('Source,n,Learner,AIPW_standardized_mean,Working_IF_SE',20),
 'outcome_model_working_if_piped_groundwater.csv':('Learner,AIPW_contrast,Working_IF_SE,Difference_from_original_HGBR,Working_IF_SE_of_difference',4),
}
def value(s): return D(s.strip().replace(',','').replace('−','-').replace('+',''))
def rows(name):
    with (AGG/name).open(encoding='utf-8',newline='') as f:
        d=csv.DictReader(f);result=list(d)
        if [a for a in (d.fieldnames or [])]!=FIELDS[name][0].split(','):
            raise AssertionError(f'Schema mismatch: {name}')
    if len(result)!=FIELDS[name][1] or any(None in r for r in result):
        raise AssertionError(f'Row count/shape mismatch: {name}')
    return result

def within(actual,expected,tol): return abs(value(str(actual))-value(str(expected)))<=D(str(tol))
class Contract(unittest.TestCase):
    def test_00_all_csv_rosters_and_fields(self):
        self.assertEqual(set(p.name for p in AGG.glob('*.csv')),set(FIELDS))
        for name in FIELDS: self.assertEqual(len(rows(name)),FIELDS[name][1])
    def test_01_covariate_schema(self):
        x=json.loads((ROOT/'data/metadata/six_characteristic_design_schema.json').read_text())
        self.assertEqual(x['N'],3066)
        self.assertEqual(x['source_labels_ordered'],SOURCES)
        self.assertEqual(len(x['features']),6)
        self.assertEqual(sum(y['number_of_indicator_columns'] for y in x['features']),33)
        self.assertEqual(x['constant_missing_indicators_expected'],5)
        self.assertEqual([y['question_id'] for y in x['features']],['Q15.17','Q15.13','Q15.15','Q15.14','Q15.10','Q15.5'])
    def test_02_means_and_pairs(self):
        m=rows('main_table_1_observed_standardized_means.csv')
        self.assertEqual([r['Primary drinking-water source'] for r in m],SOURCES)
        self.assertEqual([int(value(r['n'])) for r in m],COUNTS)
        self.assertEqual(sum(COUNTS),3066)
        means=[value(r['Standardized mean WSQS']) for r in m]
        for r in m:
            self.assertTrue(within(value(r['Observed mean WSQS'])-value(r['Standardized mean WSQS']),r['Observed − standardized'],.002))
        pair=rows('si_table_s6_ten_pairwise_mean_intervals.csv')
        self.assertEqual(len(pair),10)
        for r,(i,j) in zip(pair,itertools.combinations(range(5),2)):
            self.assertTrue(within(means[i]-means[j],r['Difference (WSQS points)'],.0011))
            self.assertEqual(r['Contrast (first source minus second)'],SOURCES[i]+' − '+SOURCES[j])
            low,high=self._bounds(r['Simultaneous 95% interval'])
            self.assertLessEqual(low,value(r['Difference (WSQS points)']))
            self.assertGreaterEqual(high,value(r['Difference (WSQS points)']))
    @staticmethod
    def _bounds(s):
        m=re.fullmatch(r'\[([−+0-9.]+),\s*([−+0-9.]+)\]',s)
        if not m: raise AssertionError('Malformed interval: '+s)
        a,b=value(m.group(1)),value(m.group(2))
        if a>b: raise AssertionError('Inverted interval: '+s)
        return a,b
    def test_03_cdf(self):
        c=rows('full_41_threshold_standardized_cdf_percent.csv')
        self.assertEqual([value(r['WSQS_threshold']) for r in c],[D('2.5')*i for i in range(41)])
        cols=['CDF_PIPED_percent','CDF_GROUNDWATER_percent','CDF_SURFACE_percent','CDF_TANKER_percent','CDF_COMMERCIAL_percent']
        v=[[value(r[k]) for r in c] for k in cols]
        for g in v:
            self.assertTrue(all(D(0)<=x<=D(100) for x in g))
            self.assertTrue(all(a<=b for a,b in zip(g,g[1:])))
            self.assertEqual(g[-1],D(100))
        for ix,r in enumerate(c): self.assertTrue(within(v[0][ix]-v[1][ix],r['PIPED_minus_GROUNDWATER_pp'],.0011))
        main=rows('main_table_2_standardized_threshold_probabilities.csv')
        self.assertEqual([r['Primary drinking-water source'] for r in main],SOURCES)
        for j,row in enumerate(main):
            for threshold in (50,60,70):
                self.assertEqual(v[j][threshold*2//5],value(row[f'WSQS ≤{threshold} (%)']))
        expected={(1,0):.93,(2,0):1.71,(2,1):.92,(3,0):3.95,(3,1):4.62,(3,2):5.12,(4,0):3.17,(4,1):3.78,(4,2):4.64,(4,3):4.02}
        for (i,j),target in expected.items():
            distance=sum((abs(v[i][k]-v[j][k])+abs(v[i][k+1]-v[j][k+1]))*D('2.5')/D('200') for k in range(40))
            self.assertTrue(within(distance,target,.006),f'{(i,j)}: {distance}, expected ~{target}')
    def test_04_threshold_contrasts(self):
        m=rows('main_table_1_observed_standardized_means.csv')
        t=rows('main_table_2_standardized_threshold_probabilities.csv')
        s=rows('si_table_s5_four_contrast_working_intervals.csv')
        expected=[value(m[0]['Standardized mean WSQS'])-value(m[1]['Standardized mean WSQS'])]+[
            value(t[0][f'WSQS ≤{p} (%)'])-value(t[1][f'WSQS ≤{p} (%)']) for p in (50,60,70)]
        self.assertEqual(len(expected),len(s))
        for x,row in zip(expected,s):
            self.assertTrue(within(x,row['Estimate'],.0011))
            for k in ('Marginal 95% interval','Four-contrast simultaneous 95% interval'):
                a,b=self._bounds(row[k]);self.assertLessEqual(a,value(row['Estimate']));self.assertLessEqual(value(row['Estimate']),b)
    def test_05_domestic_comparison_and_sensitivity(self):
        m=rows('main_table_1_observed_standardized_means.csv')
        a=rows('si_table_s3_drinking_domestic_source_comparison.csv')
        b=rows('si_table_s4_domestic_source_sensitivity.csv')
        self.assertEqual([r['Source'] for r in a],SOURCES)
        self.assertEqual(sum(int(value(r['Domestic n'])) for r in a),3044)
        self.assertEqual(3066-3044,22)
        self.assertEqual([int(value(r['Drinking n'])) for r in a],COUNTS)
        same=0
        for i,r in enumerate(a):
            self.assertEqual(value(r['Drinking observed WSQS']),value(m[i]['Observed mean WSQS']))
            matched=int(r['Same source n (%)'].split()[0].replace(',',''));same+=matched
            agreement=value(r['Same source n (%)'].split('(')[1].rstrip(')'))
            self.assertTrue(within(D(100)*matched/COUNTS[i],agreement,.051))
            self.assertEqual(value(b[i]['Drinking standardized WSQS']),value(m[i]['Standardized mean WSQS']))
            delta=value(b[i]['Domestic standardized WSQS'])-value(b[i]['Drinking standardized WSQS'])
            self.assertTrue(within(delta,b[i]['Domestic − Drinking'],.002))
        self.assertEqual(same,2012)
    def test_06_sensitivity_and_mapping(self):
        m=rows('main_table_1_observed_standardized_means.csv')
        eight=rows('eight_covariate_sensitivity.csv')
        dom=rows('si_table_s7_core_domain_sensitivity.csv')
        for i in range(5):
            self.assertEqual(value(m[i]['Standardized mean WSQS']),value(eight[i]['six_X_standardized_WSQS']))
            self.assertEqual(value(m[i]['Standardized mean WSQS']),value(dom[i]['CORE-8']))
            self.assertTrue(within(value(eight[i]['eight_X_standardized_WSQS'])-value(eight[i]['six_X_standardized_WSQS']),eight[i]['eight_X_minus_six_X_WSQS'],.0011))
        q=rows('si_table_s2_questionnaire_item_mapping.csv')
        self.assertEqual(len({r['QID'] for r in q}),8)
        self.assertEqual({r['QID'] for r in q},{'Q1.1','Q1.3','Q15.17','Q15.13','Q15.15','Q15.14','Q15.10','Q15.5'})
    def test_07_weights_and_sd(self):
        overlap=rows('overlap_diagnostics_original_aggregate.csv')
        smd=rows('figure3_source_level_max_absolute_smd.csv')
        sds=rows('observed_wsqs_sd_by_source.csv')
        means=rows('main_table_1_observed_standardized_means.csv')
        self.assertEqual([x['source'] for x in smd],SOURCES)
        self.assertEqual([x['source'] for x in sds],SOURCES)
        self.assertEqual([int(x['n']) for x in overlap],COUNTS)
        self.assertEqual([int(x['n']) for x in smd],COUNTS)
        self.assertEqual([int(x['n']) for x in sds],COUNTS)
        for i in range(5):
            self.assertTrue(within(value(overlap[i]['ess'])/COUNTS[i],overlap[i]['ess_fraction'],.00011))
            self.assertGreater(value(overlap[i]['ess']),0)
            self.assertLessEqual(value(overlap[i]['ess']),COUNTS[i]+1)
            self.assertLessEqual(value(smd[i]['max_abs_SMD_after']),value(smd[i]['max_abs_SMD_before']))
            self.assertGreater(value(sds[i]['observed_WSQS_sample_SD']),0)
        # ANOVA from rounded source means and precise SD: expected only approximately.
        overall=sum(D(COUNTS[i])*value(means[i]['Observed mean WSQS']) for i in range(5))/3066
        ss_between=sum(COUNTS[i]*(value(means[i]['Observed mean WSQS'])-overall)**2 for i in range(5))
        ss_within=sum((COUNTS[i]-1)*value(sds[i]['observed_WSQS_sample_SD'])**2 for i in range(5))
        percent=100*ss_between/(ss_between+ss_within)
        self.assertTrue(within(percent,D('0.489748347'),.002))
    def test_08_simulations(self):
        x=rows('synthetic_simulation_coverage_summary.csv')
        self.assertEqual({r['scenario'] for r in x},{'Synthetic null','Synthetic alternative'})
        self.assertEqual({r['method'] for r in x},{'Studentized max statistic','Unstudentized absolute max error'})
        self.assertEqual(sorted([value(r['full_grid_coverage_percent']) for r in x]),list(map(D,['42.500','57.900','95.000','95.800'])))
        self.assertTrue(all(value(r['nominal_percent'])==95 for r in x))
if __name__=='__main__':
    unittest.main(verbosity=2)
