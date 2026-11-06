#!/usr/bin/env python3
"""Check all release files (excluding SHA256SUMS.txt itself) against manifest."""
from pathlib import Path
import hashlib
import re

ROOT=Path(__file__).resolve().parents[1]
manifest=ROOT/'SHA256SUMS.txt'
if not manifest.is_file(): raise SystemExit('Missing SHA256SUMS.txt')
seen=set()
for line in manifest.read_text(encoding='utf-8').splitlines():
    parts=line.split('  ',1)
    if len(parts)!=2: raise SystemExit('Malformed SHA256SUMS entry')
    sha,rel=parts
    if not re.fullmatch(r'[a-f0-9]{64}',sha): raise SystemExit('Bad SHA256')
    p=ROOT/rel
    if rel in seen or rel.startswith('/') or '..' in Path(rel).parts or not p.is_file() or p.is_symlink():
        raise SystemExit('Invalid file manifest entry: '+rel)
    seen.add(rel)
    if hashlib.sha256(p.read_bytes()).hexdigest()!=sha:
        raise SystemExit('Checksum mismatch: '+rel)
actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and p!=manifest}
if actual!=seen: raise SystemExit('Missing/unmanifested release files: '+repr(sorted(actual^seen)))
print(f'PASS: full SHA-256 manifest ({len(seen)} files)')
