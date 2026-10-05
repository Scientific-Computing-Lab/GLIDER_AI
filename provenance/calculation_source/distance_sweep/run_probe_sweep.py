#!/usr/bin/env python3
"""Checkpointed atom-sampled and full quantum frozen-Coulomb probe sweeps."""
import argparse
import json
import math
import sys
import time
from pathlib import Path
import numpy as np
from ase.io import read
from pyscf import df, gto, lib
from campaign import C, ROOT, now, write_json
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'scripts/posthoc/downstream_response'))
sys.path.insert(0,str(ROOT/'scripts/posthoc/downstream_response/confirmatory'))
from common import molecule, gaussian_multipole_potential, BOHR_TO_ANGSTROM, HARTREE_TO_KCAL_MOL
from multipole_control import raw_charge_moments
from glider.inference import esp_from_sites
from glider.inference.surface import VDW_ANGSTROM

def load(path):
    with np.load(ROOT/path) as z: return {k:z[k] for k in z.files}

def surface_oxygen(base,direction,clearance):
    center=base.get_center_of_mass(); radii=np.array([VDW_ANGSTROM[s] for s in base.get_chemical_symbols()])
    lo=0.; hi=np.max(np.linalg.norm(base.positions-center,axis=1))+max(radii)+clearance+5
    for _ in range(60):
        t=(lo+hi)/2; point=center+t*direction
        if np.min(np.linalg.norm(point-base.positions,axis=1)-radii)<clearance: lo=t
        else: hi=t
    return center+hi*direction

def multipole_potential(points,origin,raw,order):
    r=(points-origin)/BOHR_TO_ANGSTROM; norm=np.linalg.norm(r,axis=1)
    potential=r@raw['dipole']/norm**3
    if order>=2:
        m=raw['second']; potential+=(3*np.einsum('pa,ab,pb->p',r,m,r)-norm**2*np.trace(m))/(2*norm**5)
    if order>=3:
        m=raw['third']; potential+=(5*np.einsum('pa,pb,pc,abc->p',r,r,r,m)-3*norm**2*np.einsum('pa,abb->p',r,m))/(2*norm**7)
    return potential

def potential_values(mol,dm,base,g,m,raw,origin,points):
    ref=np.empty(len(points))
    for start in range(0,len(points),32):
        p=points[start:start+32]/BOHR_TO_ANGSTROM
        integrals=df.incore.aux_e2(mol,gto.fakemol_for_charges(p))
        ref[start:start+len(p)]=-np.einsum('ijp,ij->p',integrals,dm)
    return {'qm':ref,'glider':esp_from_sites(base.positions,points,g['predicted_charges_e'],g['predicted_dipoles_e_bohr']),'mace_polar_l':gaussian_multipole_potential(base.positions,points,m['predicted_charges_e'],m['predicted_dipoles_e_bohr']),'exact_dipole':multipole_potential(points,origin,raw,1),'exact_dipole_quadrupole':multipole_potential(points,origin,raw,2),'exact_dipole_quadrupole_octupole':multipole_potential(points,origin,raw,3),'zero':np.zeros(len(points))}

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--config',required=True); parser.add_argument('--mode',choices=['dense','quantum']); parser.add_argument('--job-id',required=True); args=parser.parse_args()
    cfg=json.loads(Path(args.config).read_text()); jid=args.job_id; out=C/'results'/jid; out.mkdir(parents=True,exist_ok=True)
    lib.num_threads(2 if args.mode=='dense' else 4)
    base=read(ROOT/cfg['base_geometry'],cfg['base_index']); b=load(cfg['base_density']); g=load(cfg['glider']); m=load(cfg['mace'])
    mol=molecule(base.get_chemical_symbols(),base.positions); origin=base.get_center_of_mass()
    raw=raw_charge_moments(mol,b['response_density_matrix'],origin)
    np.savez_compressed(out/'base_moments.npz',origin_angstrom=origin,**raw)
    orientations=[(r,load(r['file'])) for r in cfg['orientations']]
    if args.mode=='quantum':
        from evaluate_couplings import evaluate_geometry
        torque=[(r,load(r['file'])) for r in cfg['torque_variants']]
    records=[]; started=time.perf_counter()
    for di,direction in enumerate(np.asarray(cfg['directions'])):
        if args.mode=='quantum' and di not in cfg['quantum_direction_indices']: continue
        for ci,clearance in enumerate(cfg['surface_clearances_angstrom']):
            if args.mode=='quantum' and ci not in cfg['quantum_distance_indices']: continue
            oxygen=surface_oxygen(base,direction,clearance); stem=f'd{di:02d}_r{ci:02d}'
            marker=out/f'{stem}.json'
            if marker.exists(): records.append(json.loads(marker.read_text())); continue
            accepted=[]; geometries=[]
            for oi,(meta,v) in enumerate(orientations):
                if args.mode=='quantum' and oi not in cfg['quantum_orientation_indices']: continue
                positions=v['positions_angstrom']-v['positions_angstrom'][0]+oxygen
                min_distance=float(np.min(np.linalg.norm(positions[:,None]-base.positions[None],axis=2)))
                if min_distance<cfg['minimum_probe_atom_distance_angstrom']: continue
                geometries.append(positions); accepted.append((oi,meta,v,min_distance))
            detail={'direction_index':di,'surface_clearance_angstrom':clearance,'oxygen_angstrom':oxygen.tolist(),'accepted':len(accepted),'mode':args.mode,'orientation_ids':[v[1]['id'] for v in accepted]}
            if accepted and args.mode=='dense':
                positions=np.asarray(geometries); points=positions.reshape(-1,3); h=cfg['fd_step_angstrom']
                values=potential_values(mol,b['response_density_matrix'],base,g,m,raw,origin,points)
                field={key:np.empty((len(points),3)) for key in values}
                for axis in range(3):
                    delta=np.eye(3)[axis]*h
                    plus=potential_values(mol,b['response_density_matrix'],base,g,m,raw,origin,points+delta)
                    minus=potential_values(mol,b['response_density_matrix'],base,g,m,raw,origin,points-delta)
                    for key in field: field[key][:,axis]=-(plus[key]-minus[key])/(2*h)
                result={'probe_positions_angstrom':positions,'orientation_ids':np.asarray(detail['orientation_ids']),'probe_charges_e':np.array([-.834,.417,.417]),'minimum_atom_distance_angstrom':np.array([v[3] for v in accepted])}
                for key,value in values.items():
                    esp=value.reshape(-1,3); ef=field[key].reshape(-1,3,3); force=ef*np.array([-.834,.417,.417])[None,:,None]
                    result[key+'_esp_hartree_per_e']=esp; result[key+'_field_hartree_per_e_angstrom']=ef
                    result[key+'_energy_kcal_mol']=(esp@np.array([-.834,.417,.417]))*HARTREE_TO_KCAL_MOL
                    result[key+'_force_kcal_mol_angstrom']=force.sum(1)*HARTREE_TO_KCAL_MOL
                    result[key+'_torque_kcal_mol']=np.cross(positions-oxygen,force).sum(1)*HARTREE_TO_KCAL_MOL
                np.savez_compressed(out/f'{stem}.npz',**result)
            elif accepted:
                energies=[]
                def evaluate_cached(tag,positions,density):
                    path=out/f'{stem}_{tag}.json'
                    if path.exists(): return json.loads(path.read_text())
                    val=evaluate_geometry(base=base,probe_positions=positions,probe_density=density,density_values=b,glider_values=g,mace_values=m)
                    val={k:v for k,v in val.items() if math.isfinite(v)}
                    write_json(path,val); return val
                for (oi,meta,v,_),positions in zip(accepted,geometries):
                    energies.append({'orientation_id':meta['id'],**evaluate_cached(f'o{oi:02d}',positions,v['density_matrix'])})
                # Cached +/-0.1 degree W4 densities give exact frozen-probe torque
                # at the original orientation, translated to each new oxygen position.
                derivative={}
                for meta,v in torque:
                    positions=v['positions_angstrom']-v['positions_angstrom'][0]+oxygen
                    if np.min(np.linalg.norm(positions[:,None]-base.positions[None],axis=2))<cfg['minimum_probe_atom_distance_angstrom']: continue
                    derivative[(meta['axis'],meta['sign'])]=evaluate_cached(meta['id'],positions,v['density_matrix'])
                torques={}
                for axis in range(3):
                    if (axis,-1) not in derivative or (axis,1) not in derivative: continue
                    for key in derivative[(axis,1)]:
                        if key.endswith('_energy_kcal_mol'):
                            torques.setdefault(key.replace('_energy_kcal_mol','_torque_kcal_mol'),[None]*3)[axis]=-(derivative[(axis,1)][key]-derivative[(axis,-1)][key])/(2*np.deg2rad(.1))
                detail.update(energies=energies,torques=torques)
            write_json(marker,detail); records.append(detail)
            print(f'{now()} {jid} {stem} orientations={len(accepted)}',flush=True)
    write_json(out/'success.json',{'success':True,'job_id':jid,'completed_utc':now(),'runtime_seconds':time.perf_counter()-started,'n_positions':len(records),'n_accepted_probe_geometries':sum(r['accepted'] for r in records),'source_density_reused':True,'probe_model':'TIP3P atom-sampled Coulomb' if args.mode=='dense' else 'full isolated quantum-water density frozen Coulomb','required_outputs':[str((out/'base_moments.npz').relative_to(ROOT))]})
    print('CAMPAIGN_JOB_SUCCESS '+jid,flush=True)

if __name__=='__main__': main()
