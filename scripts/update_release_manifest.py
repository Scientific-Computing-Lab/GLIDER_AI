#!/usr/bin/env python3
"""Maintainer operation: rebuild hashes after reviewing intentional release changes."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
files={}
for directory in ['experiments','checkpoints','figures','assets','glider','scripts','provenance']:
    for p in sorted((ROOT/directory).rglob('*')):
        if not p.is_file() or p.name=='README.md' or p.suffix=='.pyc' or '__pycache__' in p.parts:continue
        if p.stat().st_size>=100*1024*1024:raise RuntimeError(f'GitHub per-file limit: {p}')
        files[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
(ROOT/'RELEASE_MANIFEST.json').write_text(json.dumps(dict(format='glider-experiment-release-v2',revision='non-water contact correction 2026-10-07',files=files),indent=2)+'\n')
print('Manifest files:',len(files))
