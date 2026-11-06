#!/usr/bin/env python3
"""Exact transcription checks for original November 3 Q1.1 sensitivity aggregates."""
from pathlib import Path
import csv
from decimal import Decimal
ROOT=Path(__file__).resolve().parents[1]
AGG=ROOT/'data'/'aggregates'
def read(name):
    with (AGG/name).open(newline='',encoding='utf-8') as f:
        return list(csv.DictReader(f))
def D(x):return Decimal(x.replace('−','-').replace('+','').replace(',','').strip())
names=['Municipal piped','Groundwater','Surface water','Tanker delivery','Commercial water']
rows=read('si_table_s3_drinking_domestic_source_comparison.csv')
sens=read('si_table_s4_domestic_source_sensitivity.csv')
assert [r['Source'] for r in rows]==names==[r['Source'] for r in sens]
assert [int(D(r['Domestic n'])) for r in rows]==[1335,1449,160,50,50]
assert [int(D(r['Drinking n'])) for r in rows]==[1275,1212,234,120,225]
assert [int(r['Same source n (%)'].split('(')[0].strip()) for r in rows]==[960,952,62,14,24]
assert sum(int(r['Same source n (%)'].split('(')[0].strip()) for r in rows)==2012
assert sum(int(D(r['Domestic n'])) for r in rows)==3044 and 3066-3044==22
expected_dom=[D(x) for x in ('74.908','75.148','73.363','77.510','67.759')]
assert [D(r['Domestic standardized WSQS']) for r in sens]==expected_dom
assert [D(r['Domestic ESS/n (%)']) for r in sens]==[D(x) for x in ('96.7','97.3','89.6','41.2','65.3')]
assert [D(r['Maximum balance gap (pp)']) for r in sens]==[D(x) for x in ('0.38','0.39','2.27','17.51','8.11')]
primary=read('main_table_1_observed_standardized_means.csv')
assert [D(r['Drinking standardized WSQS']) for r in sens]==[D(r['Standardized mean WSQS']) for r in primary]
assert abs((expected_dom[0]-expected_dom[1])-D('-0.240610'))<=D('0.0011')
assert abs((D(primary[0]['Standardized mean WSQS'])-D(primary[1]['Standardized mean WSQS']))-D('-0.601800'))<=D('0.0011')
assert abs(((expected_dom[0]-expected_dom[1])-(D(primary[0]['Standardized mean WSQS'])-D(primary[1]['Standardized mean WSQS'])))-D('0.361190'))<=D('0.0021')
print('PASS — Q1.1 group counts, matched-source records, 5 standardized means and overlap diagnostics')
print('NOTICE — values independently check released transcription; not a household-data refit')
