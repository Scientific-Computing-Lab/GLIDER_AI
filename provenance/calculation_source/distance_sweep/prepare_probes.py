#!/usr/bin/env python3
"""Freeze dense surface probes and a quantum-water subset from cached W4 densities."""
import json
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
from ase.io import read
from campaign import C, ROOT, write_json
from prepare_embedding import QM
import sys
sys.path.insert(0,str(ROOT))
from glider.inference.surface import fibonacci_sphere, VDW_ANGSTROM

SOURCE=ROOT/'results/posthoc_downstream_response_confirmatory'

def main():
    if (C/'configs/probe_manifest.json').exists(): return
    original=C/'data/original_w4'; original.mkdir(exist_ok=True)
    for name in ['base_3water.extxyz','heldout_w4.extxyz','outer_w4_w12.extxyz','configurations.csv']:
        shutil.copy2(SOURCE/'configurations'/name,original/name)
    bases=read(SOURCE/'configurations/base_3water.extxyz',':')
    variants=pd.read_csv(SOURCE/'qm/probe_variants/variant_registry.csv')
    bd=pd.read_csv(SOURCE/'qm/base_densities/density_registry.csv').set_index('case_id')
    gd=pd.read_csv(SOURCE/'predictions/glider/prediction_registry.csv').set_index('config_id')
    md=pd.read_csv(SOURCE/'predictions/mace_polar_1_l/prediction_registry.csv').set_index('config_id')
    rows=[]
    for bi,base in enumerate(bases):
        case=str(base.info['case_id']); orientations=variants[(variants.case_id==case)&(variants.water_rank==4)&(variants.variant_type=='orientation')].sort_values('variant_id')
        orientation_rows=orientations.iloc[::2].to_dict('records')
        torque_rows=variants[(variants.case_id==case)&(variants.water_rank==4)&(variants.variant_type=='torque_primary')].to_dict('records')
        cfg={'base_index':bi,'case_id':case,'molecule_id':str(base.info['molecule_id']),'base_geometry':str((original/'base_3water.extxyz').relative_to(ROOT)),'base_density':str((SOURCE/'qm/base_densities'/bd.loc[case,'density_file']).relative_to(ROOT)),'glider':str((SOURCE/'predictions/glider'/gd.loc[case,'prediction_file']).relative_to(ROOT)),'mace':str((SOURCE/'predictions/mace_polar_1_l'/md.loc[case,'prediction_file']).relative_to(ROOT)),'orientations':[{'id':r['variant_id'],'file':str((SOURCE/'qm/probe_variants'/r['density_file']).relative_to(ROOT))} for r in orientation_rows],'torque_variants':[{'id':r['variant_id'],'axis':int(r['axis']),'sign':int(r['sign']),'step_degrees':float(r['step_degrees']),'file':str((SOURCE/'qm/probe_variants'/r['density_file']).relative_to(ROOT))} for r in torque_rows],'directions':fibonacci_sphere(12).tolist(),'surface_clearances_angstrom':[0.8,1.2,1.8,2.5,3.5,5.,7.5,10.,15.],'probe_charge_model':'TIP3P only for dense atom-sampled coupling; separate exact quantum-water subset uses cached isolated water densities','minimum_probe_atom_distance_angstrom':1.2,'fd_step_angstrom':0.001,'quantum_direction_indices':[0,4,8],'quantum_distance_indices':[0,2,4,6,8],'quantum_orientation_indices':[0,3,6,9]}
        path=C/'configs'/f'A_base{bi:02d}.json'; write_json(path,cfg)
        for mode,resource,priority in [('dense','cpu',20),('quantum','gpu',20)]:
            jid=f'A_{mode}_{bi:02d}'
            write_json(C/'jobs'/f'{jid}.json',{'id':jid,'family':f'A: W4 {mode} sweep','kind':'probe','resource':resource,'priority':priority,'depends':[],'command':[QM,str(C/'scripts/run_probe_sweep.py'),'--config',str(path),'--mode',mode,'--job-id',jid],'success':f'research_campaign/results/{jid}/success.json','timeout_seconds':43200,'max_attempts':2})
        rows.append(cfg)
    write_json(C/'configs/probe_manifest.json',{'n_bases':len(rows),'dense_candidates':len(rows)*12*9*12,'quantum_candidates':len(rows)*3*5*4,'original_geometries':str(original.relative_to(ROOT)),'bases':rows,'no_new_reference_scf_required':True,'frozen_base_response':True,'distance_definition':'oxygen minimum distance to union of base atom vdW surfaces; H/atom steric filter applied separately','original_scientific_results_untouched':True})
    print(f'Queued {len(rows)} dense CPU and {len(rows)} exact quantum-water GPU sweep jobs')

def preflight():
    path=C/'configs/A_quantum_preflight.json'
    if path.exists(): return
    cfg=json.loads((C/'configs/A_base00.json').read_text())
    cfg.update(quantum_direction_indices=[0],quantum_distance_indices=[4],quantum_orientation_indices=[0],torque_variants=[])
    write_json(path,cfg)
    jid='A_quantum_preflight'
    write_json(C/'jobs'/f'{jid}.json',{'id':jid,'family':'A: quantum-water preflight','kind':'probe','resource':'gpu','priority':9,'depends':[],'command':[QM,str(C/'scripts/run_probe_sweep.py'),'--config',str(path),'--mode','quantum','--job-id',jid],'success':f'research_campaign/results/{jid}/success.json','timeout_seconds':1800,'max_attempts':2})

if __name__=='__main__':
    main()
    preflight()
