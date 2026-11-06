#!/usr/bin/env python3
"""Run aggregate and code-integrity checks."""
from pathlib import Path
import subprocess,sys,os
ROOT=Path(__file__).resolve().parents[1]
env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
for name in (
    'scripts/validate_public_package.py',
    'tests/test_data_contract.py',
    'tests/test_domestic_source.py',
    'tests/test_e20_scope.py',
    'scripts/release_screen.py',
    'scripts/verify_original_code_archive.py',
    'scripts/verify_source_design_contract.py',
    'scripts/verify_file_manifest.py',
    'tests/test_rejection_cases.py',
):
    print('CHECK',name,flush=True)
    subprocess.run([sys.executable,'-B',str(ROOT/name)],cwd=ROOT,check=True,env=env)
print('PASS: aggregate results and code-integrity checks complete')
