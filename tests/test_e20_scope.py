#!/usr/bin/env python3
"""Verify the scope metadata for the archived E20 summary without fitting models."""
import csv,json,re
from pathlib import Path
root=Path(__file__).resolve().parents[1]
meta=json.loads((root/'data/metadata/e20_coverage_scope.json').read_text())
source=(root/meta['reference_program']).read_text()
with (root/meta['table']).open(newline='',encoding='utf-8') as f:rows=list(csv.DictReader(f))
assert len(rows)==4
assert meta['comparison']=='PIPED_TANKER' and "'PIPED_TANKER':(0,3)" in source
assert meta['threshold_family']=='FULL_INTERIOR_39' and "'FULL_INTERIOR_39':np.arange(1,40)" in source
assert meta['threshold_count']==39 and meta['threshold_min_wsqs']==2.5 and meta['threshold_max_wsqs']==97.5
assert meta['null_replications']==240 and "'E20_NULL_REPS','240'" in source
assert meta['alternative_replications']==120 and "'E20_ALT_REPS','120'" in source
assert meta['multiplier_draws']==1000 and "'E20_MULTIPLIER_DRAWS','1000'" in source
for r in rows:
    scenario={'Synthetic null':'null','Synthetic alternative':'alternative'}[r['scenario']]
    field='studentized_coverage_percent' if r['method']=='Studentized max statistic' else 'absolute_max_error_coverage_percent'
    assert abs(float(r['full_grid_coverage_percent'])-meta[field][scenario])<1e-9
    assert float(r['nominal_percent'])==95
print('PASS: E20 simulation summary, 39-interior-threshold scope, and original code configuration')
