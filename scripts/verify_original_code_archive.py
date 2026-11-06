#!/usr/bin/env python3
"""Verify source-file hashes and Python syntax without executing programs."""
from pathlib import Path
import hashlib, ast, re
ROOT=Path(__file__).resolve().parents[1]
EXPECTED={
    'code/original/E04/run_aipw.py': '439fc907ff72fdc93071bb61e4ea37cd0ecff68b2948a2fd6942208bcf9d9bf1',
'code/original/E04_engine/core.py': '10fafe121b64b9ff53bcb992c076d45a7b1097fd0ba76e7b1702dbb2e0cd4645',
'code/original/E09/logistic_inference.py': 'f12f0240cf929003dabc6eab188f52ca28aa0a78a7add186852d562cc6014222',
'code/original/E15/run_cdf.py': '3277b60e82d4d5e4714ddd7d15e4b6478797284391dd8a6b320db8404346708d',
'code/original/E20/simulate_holdout.py': '42149a134cfba1c3dc63aa594aa74f0731bcf0d5d9de7a40a082920a2737effe',
'code/original/E23/reconcile.py': 'df85dbcd7f6d6e8d0371590ddac1c9cfe8bc3cc8f136e776b306f9479056d3d2',
'code/original/S1_generator_1/build_si_figures.py': '6985a2ea05c7782cfb80bd1881f66f2e9bb529d3bc0e87adaa42848a28895e11',
'code/original/dependency_demo/demo.py': 'eec9c3e94b33085b6784ce8d65622ab9e33284d7f1c7fdc72ffb5abb8bf8d39b',
'code/original/dependency_support_policy/support_policy.py': 'd2ae5f0484e58d216dcc2deead8eb27a2452fa2b09cac562a8f34f7f3182beca',

}
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):
            h.update(chunk)
    return h.hexdigest()
actual={p.relative_to(ROOT).as_posix() for p in (ROOT/'code/original').rglob('*.py')}
assert actual==set(EXPECTED),f'Original code roster changed: {sorted(actual.symmetric_difference(EXPECTED))}'
for relative,expected in EXPECTED.items():
    p=ROOT/relative
    assert sha(p)==expected,f'Original source SHA mismatch: {relative}'
    ast.parse(p.read_text(encoding='utf-8-sig'),filename=relative)
for module,source in {
    'core.py':'E04_engine/core.py',
    'demo.py':'dependency_demo/demo.py',
    'support_policy.py':'dependency_support_policy/support_policy.py',
}.items():
    original=ROOT/'code/original'/source
    copy=ROOT/'code/runner_support/src/water_standardize'/module
    assert copy.read_bytes()==original.read_bytes(),f'Import-support module altered: {module}'
init=ROOT/'code/runner_support/src/water_standardize/__init__.py'
assert init.is_file(), 'Support package initializer missing'
t=ast.parse(init.read_text(encoding='utf-8'))
assert len(t.body)==1 and isinstance(t.body[0],ast.Expr),'Packaging shim contains executable statements'
manifest=ROOT/'code/ORIGINAL_SOURCE_SHA256SUMS.txt'
expected_manifest=''.join(f'{v}  {name}\n' for name,v in sorted(EXPECTED.items()))
assert manifest.read_text(encoding='utf-8')==expected_manifest, 'Source fingerprint register mismatch'
print(f'PASS: {len(EXPECTED)} source files verified for SHA-256 and syntax')
print('PASS: 3 support modules and package initializer checked')
print('PASS: source-file verification complete')
