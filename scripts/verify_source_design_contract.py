#!/usr/bin/env python3
"""Read E04 source declarations and compare them with the six-variable design schema."""
import ast, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
src=ROOT/'code/original/E04/run_aipw.py'
meta=json.loads((ROOT/'data/metadata/six_characteristic_design_schema.json').read_text(encoding='utf-8'))
tree=ast.parse(src.read_text(encoding='utf-8'))
required={'QIDS','categories','regime_names'}
found={}
for node in tree.body:
    if isinstance(node,ast.Assign):
        for target in node.targets:
            if isinstance(target,ast.Name) and target.id in required:
                found[target.id]=ast.literal_eval(node.value)
assert set(found)==required, f'Missing original source declarations: {required-set(found)}'
qids=tuple(item['question_id'] for item in meta['features'])
assert qids==tuple(found['QIDS']), 'Six-X question sequence differs from E04'
assert set(found['categories'])==set(qids), 'E04 category dictionary differs from schema'
assert list(found['regime_names'].values())==['PIPED','GROUNDWATER','SURFACE','TANKER','COMMERCIAL'], 'Source order differs from E04'
assert meta['source_integer_codes']==list(found['regime_names']), 'Source codes differ from E04'
for item in meta['features']:
    n=len(found['categories'][item['question_id']])+1 # explicit missing category
    assert n==item['number_of_indicator_columns']==len(item['indicator_labels']), 'Indicator count differs for '+item['question_id']
assert sum(item['number_of_indicator_columns'] for item in meta['features'])==meta['total_binary_columns']==33
assert meta['N']==3066
print('PASS: six-X question order, categorical cardinalities, 33 indicators and source code/order match E04 source definitions')
