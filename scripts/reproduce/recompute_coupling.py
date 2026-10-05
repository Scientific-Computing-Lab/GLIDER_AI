#!/usr/bin/env python3
"""Recompute one held-out-water energy from geometry, QM densities and sites.

This CPU path evaluates Coulomb integrals, not SCF. It uses the same density-fit
definition as the archived GPU evaluation. Large solutes may require substantial
memory. Methane is the recommended small end-to-end example.
"""
from pathlib import Path
import argparse
import importlib.util
import sys
import json
import numpy as np
import pandas as pd
from ase.io import read
from pyscf import lib

ROOT=Path(__file__).resolve().parents[2]


def array(path):
    with np.load(path) as z:return {k:z[k] for k in z.files}


def lookup(directory,identifier,value,registry='density_registry.csv',column='density_file'):
    rows=pd.read_csv(directory/registry)
    selected=rows[rows[identifier]==value]
    if len(selected)!=1:raise ValueError((directory,identifier,value,len(selected)))
    return array(directory/selected.iloc[0][column])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solute',default='methane')
    p.add_argument('--water-rank',type=int,default=4,choices=range(4,13))
    p.add_argument('--threads',type=int,default=2)
    p.add_argument('--output',type=Path,default=ROOT/'build/coupling_check.json')
    a=p.parse_args();lib.num_threads(a.threads)
    directory=ROOT/'provenance/calculation_source/downstream_response'
    sys.path.insert(0,str(directory))
    spec=importlib.util.spec_from_file_location('archived_cpu_coupling',directory/'evaluate_couplings.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    e=ROOT/'experiments/heldout_water';case=f'confirmatory__{a.solute}__base3';probe_id=f'{case}__W{a.water_rank}'
    base=next(x for x in read(e/'geometries/base_3water.extxyz',index=':') if x.info['case_id']==case)
    probe=next(x for x in read(e/'geometries/outer_w4_w12.extxyz',index=':') if x.info['probe_id']==probe_id)
    density=lookup(e/'references/base_densities','case_id',case)
    water=lookup(e/'references/probe_densities','probe_id',probe_id)
    predictions={}
    for method in ['glider','mace_polar_1_l']:
        predictions[method]=lookup(e/'predictions'/method,'config_id',case,
            registry='prediction_registry.csv',column='prediction_file')
    computed=module.evaluate_geometry(base=base,probe_positions=probe.positions,
        probe_density=water['probe_density_matrix'],density_values=density,
        glider_values=predictions['glider'],mace_values=predictions['mace_polar_1_l'])
    stored=pd.read_csv(e/'evaluated/range_energy_per_probe.csv').set_index('probe_id').loc[probe_id]
    comparison={}
    for method in ['qm','glider','mace_polar_l']:
        key=method+'_energy_kcal_mol'
        comparison[method]=dict(recomputed=computed[key],archived=float(stored[key]),
            absolute_difference=abs(computed[key]-float(stored[key])))
    # Tolerance covers CPU/GPU density-fitting/integral differences, far below the reported errors.
    if max(x['absolute_difference'] for x in comparison.values())>1e-4:
        raise AssertionError(comparison)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(dict(probe_id=probe_id,comparison=comparison,tolerance_kcal_mol=1e-4),indent=2)+'\n')
    print(a.output.read_text())
