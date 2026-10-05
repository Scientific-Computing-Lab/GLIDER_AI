#!/usr/bin/env python3
"""Reproduce either author-review control, including frozen feature extraction."""
import argparse,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
FROZEN=ROOT/'provenance/frozen_source/prospective_1'
def run(script,*args):subprocess.run([sys.executable,str(script),*map(str,args)],cwd=ROOT,check=True)
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('kind',choices=['ablation','dissociation'])
p.add_argument('--output',type=Path,required=True)
p.add_argument('--device',default='cpu')
p.add_argument('--mace-polar-m',type=Path,default=ROOT/'third_party/checkpoints/MACE-POLAR-1-M.model')
p.add_argument('--mace-polar-l',type=Path,default=ROOT/'third_party/checkpoints/MACE-POLAR-1-L.model')
p.add_argument('--distances',nargs='+',type=float,default=[3,4,5,6,8,10,15,20])
a=p.parse_args();out=a.output.resolve();diag=ROOT/'scripts/diagnostics/review_controls.py'
if a.kind=='ablation':
    run(diag,'train','--source-root',ROOT,'--output',out,'--device',a.device)
    run(diag,'build-evaluation','--source-root',ROOT,'--output',out)
    geometry=out/'evaluation_configurations.extxyz';features=out/'evaluation_features.npz';prior=out/'evaluation_prior.npz'
else:
    run(diag,'build-dissociation','--source-root',ROOT,'--output',out,'--distances',*a.distances)
    geometry=out/'configurations.extxyz';features=out/'mace_polar_m_features.npz';prior=out/'mace_polar_ml_average.npz'
run(FROZEN/'extract_mace_features.py','--configurations',geometry,'--checkpoint',a.mace_polar_m.resolve(),'--output',features,'--device',a.device)
for name,checkpoint in [('m',a.mace_polar_m),('l',a.mace_polar_l)]:
    run(FROZEN/'predict_mace_polar_cp_response.py','--configurations',geometry,'--checkpoint',checkpoint.resolve(),'--output',out/f'prior_{name}.npz','--device',a.device)
run(FROZEN/'average_qu_predictions.py','--first',out/'prior_m.npz','--second',out/'prior_l.npz','--output',prior)
run(diag,'score-ablation' if a.kind=='ablation' else 'score-dissociation','--source-root',ROOT,'--output',out)
