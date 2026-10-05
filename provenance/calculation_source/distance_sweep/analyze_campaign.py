#!/usr/bin/env python3
"""Post hoc analysis of frozen campaign outputs; no SCF, fitting, or job launches."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from ase import Atoms
from ase.io import read
from scipy.stats import binomtest, spearmanr

ROOT = Path(__file__).resolve().parents[2]
C = ROOT / 'research_campaign'
OUT = C / 'analysis'
FIG = C / 'figures'
B = 0.529177210903
K = 627.5094740631
D = 2.541746473
SEED = 20260915
NBOOT = 20000
METHODS = ['glider','mace_polar_l','exact_dipole','exact_dipole_quadrupole',
           'exact_dipole_quadrupole_octupole','zero']
SUMMARY = []


def jread(path):
    return json.loads(Path(path).read_text())


def jsave(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def arrays(path):
    with np.load(path, allow_pickle=False) as z:
        return {k:z[k] for k in z.files}


def region(clearance):
    return 'near' if clearance<=2.5 else 'far' if clearance>=7.5 else 'intermediate'


def rms(x):
    return float(np.sqrt(np.mean(np.asarray(x)**2)))


def vector_rms(x):
    return float(np.sqrt(np.mean(np.sum(np.asarray(x)**2, axis=-1))))


def boot_indices(n):
    if n<=4:
        return np.asarray(list(itertools.product(range(n),repeat=n)))
    return np.random.default_rng(SEED+n).integers(0,n,(NBOOT,n))


def estimate(values):
    x=np.asarray(values,float)
    assert np.isfinite(x).all() and len(x)
    lo,hi=np.quantile(x[boot_indices(len(x))].mean(1),[.025,.975])
    return float(x.mean()),float(lo),float(hi)


def record(family,subset,metric,unit,method,values,comparator='',other=None,notes=''):
    x=np.asarray(values,float); mean,lo,hi=estimate(x)
    row=dict(family=family,subset=subset,metric=metric,unit=unit,method=method,
             n_systems=len(x),estimate=mean,ci_low=lo,ci_high=hi,comparator=comparator,
             paired_delta=np.nan,delta_ci_low=np.nan,delta_ci_high=np.nan,
             comparator_mean=np.nan,relative_error_reduction=np.nan,wins=np.nan,
             ties=np.nan,sign_test_p=np.nan,paired_standardized_effect=np.nan,notes=notes)
    if other is not None:
        y=np.asarray(other,float); assert x.shape==y.shape
        delta=x-y; dm,dl,dh=estimate(delta)
        wins=int(np.sum(delta<-1e-12)); losses=int(np.sum(delta>1e-12))
        row.update(paired_delta=dm,delta_ci_low=dl,delta_ci_high=dh,comparator_mean=float(y.mean()),
                   relative_error_reduction=float(1-x.mean()/y.mean()) if y.mean()!=0 else np.nan,
                   wins=wins,ties=len(x)-wins-losses,
                   sign_test_p=binomtest(wins,wins+losses).pvalue if wins+losses else 1.,
                   paired_standardized_effect=float(delta.mean()/delta.std(ddof=1)) if len(x)>1 and delta.std(ddof=1)>0 else np.nan)
    SUMMARY.append(row)
    return row


def aggregate_cases(frame,family,metrics,group_columns=(),comparators=('zero',)):
    """Equal molecule weighting; seeds averaged within molecule, never independent blocks."""
    groups=[('all',frame)]
    for col in group_columns:
        groups += [(f'{col}={key}',part) for key,part in frame.groupby(col,sort=True)]
    for subset,part in groups:
        for metric,unit in metrics.items():
            pivot=part.groupby(['molecule_id','method'])[metric].mean().unstack('method')
            for method in pivot.columns:
                valid=pivot[method].dropna()
                record(family,subset,metric,unit,method,valid.values)
                for comparator in comparators:
                    if comparator==method or comparator not in pivot: continue
                    pair=pivot[[method,comparator]].dropna()
                    if len(pair): record(family,subset,metric,unit,method,pair[method],comparator,pair[comparator])


def multipole(points,origin,mom,order):
    r=(np.asarray(points)-origin)/B; norm=np.linalg.norm(r,axis=-1)
    v=np.einsum('...i,i->...',r,mom['dipole'])/norm**3
    if order>=2:
        m=mom['second']; v+=(3*np.einsum('...a,ab,...b->...',r,m,r)-norm**2*np.trace(m))/(2*norm**5)
    if order>=3:
        m=mom['third']; v+=(5*np.einsum('...a,...b,...c,abc->...',r,r,r,m)-3*norm**2*np.einsum('...a,abb->...',r,m))/(2*norm**7)
    return v


def site_field(positions,points,q,u):
    r=np.asarray(points)[:,None,:]-np.asarray(positions)[None,:,:]
    rr=np.linalg.norm(r,axis=-1); dot=np.einsum('pni,ni->pn',r,u)
    esp=np.sum(B*q[None,:]/rr+B*B*dot/rr**3,axis=1)
    field=np.sum(B*q[None,:,None]*r/rr[:,:,None]**3+
                 B*B*(3*dot[:,:,None]*r/rr[:,:,None]**5-u[None,:,:]/rr[:,:,None]**3),axis=1)
    return esp,field


def multipole_field(points,origin,mom,order,h=.001):
    esp=multipole(points,origin,mom,order)
    field=np.column_stack([-(multipole(points+np.eye(3)[a]*h,origin,mom,order)-
                            multipole(points-np.eye(3)[a]*h,origin,mom,order))/(2*h) for a in range(3)])
    return esp,field


def finite_archive(path):
    z=arrays(path)
    for name,v in z.items():
        if v.dtype.kind in 'fci': assert np.isfinite(v).all(),(path,name)
    return z


def audit():
    jobs=[jread(p) for p in sorted((C/'jobs').glob('*.json'))]
    rows=[]; problems=[]
    for job in jobs:
        state=jread(C/'status/jobs'/f'{job["id"]}.json')
        row={'job_id':job['id'],'family':job['family'],'state':state['state'],
             'usable':state['state']=='complete','reason':state.get('error') or ''}
        if row['usable']:
            marker=jread(ROOT/job['success'])
            assert marker.get('success') and marker['job_id']==job['id'],job['id']
            for f in marker.get('required_outputs',[]): assert (ROOT/f).is_file(),f
        rows.append(row)
    manifest=jread(C/'configs/embedding_manifest.json')
    qrows=[]
    from pyscf import gto,lib
    lib.num_threads(1)
    for geom in manifest['geometries']:
        gid=geom['geometry_id']; expected=[r for r in manifest['jobs'] if r['geometry_id']==gid]
        existing=[r for r in expected if (C/'qm'/r['job_id']/'success.json').exists()]
        if not existing: continue
        mol=gto.M(atom=list(zip(geom['symbols'],geom['positions_angstrom'])),basis=geom['basis'],verbose=0)
        overlap=mol.intor('int1e_ovlp')
        vacuum=arrays(C/'qm'/f'B_{gid}_vacuum/result.npz')
        for row in existing:
            jid=row['job_id']; z=finite_archive(C/'qm'/jid/'result.npz'); marker=jread(C/'qm'/jid/'success.json')
            log=ROOT/next(p for p in marker['required_outputs'] if p.endswith('pyscf.log'))
            assert 'converged SCF energy' in log.read_text(),jid
            assert marker['scf_converged'] and marker['theory']['basis']=='def2-tzvpd'
            assert marker['theory']['xc']=='wb97x-d3bj' and marker['theory']['grid_level']==4
            assert abs(float(z['energy_hartree']-z['vacuum_energy_hartree']-z['energy_difference_hartree']))<1e-10
            assert np.allclose(z['response_density_matrix'],z['density_matrix']-vacuum['density_matrix'],rtol=0,atol=1e-12)
            assert np.allclose(z['response_esp_hartree_per_e'],z['esp_hartree_per_e']-vacuum['esp_hartree_per_e'],rtol=0,atol=1e-10)
            electrons=float(np.einsum('ij,ji->',z['density_matrix'],overlap))
            dq=float(np.einsum('ij,ji->',z['response_density_matrix'],overlap))
            assert abs(electrons-mol.nelectron)<1e-6 and abs(dq)<1e-6
            moments=finite_archive(C/'results'/f'C_{gid}'/f'{jid}.npz')
            assert np.allclose(moments['response_raw_rank1']*D,z['induced_dipole_debye'],atol=1e-7)
            assert abs(float(moments['response_raw_rank0']))<1e-6
            cfg=jread(ROOT/row['configuration'])
            assert np.allclose(z['positions_angstrom'],cfg['positions_angstrom'],atol=1e-10)
            assert np.allclose(z['points_angstrom'],vacuum['points_angstrom'],atol=1e-10)
            qrows.append({**row,'usable':True,'n_atoms':len(z['symbols']),'n_esp_points':len(z['points_angstrom']),
                          'electron_error':electrons-mol.nelectron,'response_electrons':dq,
                          'esp_signal_mEh_per_e':rms(z['response_esp_hartree_per_e'])*1000,
                          'induced_dipole_D':float(np.linalg.norm(z['induced_dipole_debye'])),
                          'energy_difference_kcal_mol':float(z['energy_difference_hartree'])*K})
        print('Validated QM',gid,flush=True)
    pd.DataFrame(rows).to_csv(OUT/'job_inventory.csv',index=False)
    pd.DataFrame(qrows).to_csv(OUT/'qm_inventory.csv',index=False)
    splits=jread(C/'configs/splits.json')['molecule_split']
    formulas={}
    for geom in manifest['geometries']:
        formulas.setdefault(geom['molecule_id'],Atoms(geom['symbols']).get_chemical_formula())
    assert len(set(formulas.values()))==len(formulas),'Check chemical aliases: duplicate molecular formula'
    usable={r['molecule_id'] for r in qrows}
    counts={s:len({m for m in usable if splits[m]==s}) for s in ['train','validation','test']}
    detail={'jobs':pd.Series([r['state'] for r in rows]).value_counts().to_dict(),
            'qm_success':len(qrows),'qm_embedded':sum(r['mode']!='vacuum' for r in qrows),
            'usable_molecules':len(usable),'usable_geometries':len({r['geometry_id'] for r in qrows}),
            'usable_molecule_splits':counts,'molecular_formulas':formulas,
            'excluded_molecule':'prospective_diverse_02','excluded_geometry':'g005',
            'excluded_qm_ids':[r['job_id'] for r in manifest['jobs'] if r['geometry_id']=='g005'],
            'test_molecules':[m for m in splits if splits[m]=='test'],'problems':problems}
    jsave(OUT/'validation_summary.json',detail)
    return detail


def w4():
    from glider.inference.surface import VDW_ANGSTROM
    dense=[]; quantum=[]; coverage=[]
    bases=read(C/'data/original_w4/base_3water.extxyz',':')
    for bi in range(10):
        cfg=jread(C/f'configs/A_base{bi:02d}.json'); mid=cfg['molecule_id']; base=bases[bi]
        raw=arrays(C/f'results/A_dense_{bi:02d}/base_moments.npz'); origin=raw['origin_angstrom']
        extent=max(np.linalg.norm(base.positions-origin,axis=1)+np.array([VDW_ANGSTROM[s] for s in base.get_chemical_symbols()]))
        for family,target in [('dense',dense),('quantum',quantum)]:
            directory=C/f'results/A_{family}_{bi:02d}'
            paths=sorted(directory.glob('d??_r??.json'))
            assert len(paths)==(108 if family=='dense' else 15)
            accepted=0
            for path in paths:
                meta=jread(path); accepted+=meta['accepted']; dist=meta['surface_clearance_angstrom']
                if not meta['accepted']: continue
                common=dict(molecule_id=mid,base_index=bi,clearance_A=dist,region=region(dist),
                            direction=meta['direction_index'],n_orientations=meta['accepted'])
                if family=='dense':
                    z=finite_archive(path.with_suffix('.npz'))
                    assert z['probe_positions_angstrom'].shape==(meta['accepted'],3,3)
                    assert (z['minimum_atom_distance_angstrom']>=1.2-1e-10).all()
                    assert list(z['orientation_ids'])==meta['orientation_ids']
                    inside=float(np.mean(np.linalg.norm(z['probe_positions_angstrom']-origin,axis=-1)<extent))
                    ref=z['qm_energy_kcal_mol']; std=float(np.std(ref)); denom=rms(ref)
                    for method in METHODS:
                        energy=z[method+'_energy_kcal_mol']; esp=z[method+'_esp_hartree_per_e']
                        assert np.allclose(energy,esp@z['probe_charges_e']*K,rtol=1e-9,atol=1e-9)
                        row={**common,'method':method,'inside_source_fraction':inside,
                             'energy_mae':float(np.mean(np.abs(energy-ref))),
                             'energy_rmse':rms(energy-ref),'energy_signal':denom,
                             'energy_nrmse':rms(energy-ref)/max(denom,1e-15),
                             'centered_orientation_rmse':rms((energy-energy.mean())-(ref-ref.mean())),
                             'orientation_signal':std,'anisotropy_fraction':std/max(denom,1e-15),
                             'esp_rmse_mEh':rms(esp-z['qm_esp_hartree_per_e'])*1000,
                             'field_vrmse':vector_rms(z[method+'_field_hartree_per_e_angstrom']-z['qm_field_hartree_per_e_angstrom'])*K,
                             'force_vrmse':vector_rms(z[method+'_force_kcal_mol_angstrom']-z['qm_force_kcal_mol_angstrom']),
                             'torque_vrmse':vector_rms(z[method+'_torque_kcal_mol']-z['qm_torque_kcal_mol']),
                             'torque_signal':vector_rms(z['qm_torque_kcal_mol'])}
                        pred_min=np.flatnonzero(np.isclose(energy,energy.min(),rtol=0,atol=1e-12))
                        row['preferred_orientation_accuracy']=float(int(np.argmin(ref)) in pred_min)/len(pred_min)
                        row['centered_orientation_nrmse']=row['centered_orientation_rmse']/max(std,1e-15)
                        target.append(row)
                else:
                    ref=np.array([r['qm_energy_kcal_mol'] for r in meta['energies']]); assert len(ref)==meta['accepted']
                    methods=['glider','mace_polar_l','exact_dipole_base_com','exact_dipole_solute_com','exact_dipole_nuclear_center','zero']
                    for method in methods:
                        energy=np.array([r[method+'_energy_kcal_mol'] for r in meta['energies']])
                        assert np.isfinite(energy).all()
                        rt=np.array(meta['torques'].get('qm_torque_kcal_mol',[np.nan]*3),float)
                        pt=np.array(meta['torques'].get(method+'_torque_kcal_mol',[np.nan]*3),float)
                        row={**common,'method':method,'energy_mae':float(np.mean(np.abs(energy-ref))),
                             'energy_rmse':rms(energy-ref),'energy_signal':rms(ref),
                             'energy_nrmse':rms(energy-ref)/max(rms(ref),1e-15),
                             'centered_orientation_rmse':rms((energy-energy.mean())-(ref-ref.mean())),
                             'orientation_signal':float(np.std(ref)),
                             'torque_vrmse':float(np.linalg.norm(pt-rt)) if np.isfinite(rt).all() and np.isfinite(pt).all() else np.nan,
                             'torque_signal':float(np.linalg.norm(rt)) if np.isfinite(rt).all() else np.nan}
                        row['centered_orientation_nrmse']=row['centered_orientation_rmse']/max(row['orientation_signal'],1e-15)
                        target.append(row)
            marker=jread(directory/'success.json'); assert marker['n_accepted_probe_geometries']==accepted
            coverage.append(dict(family=family,molecule_id=mid,n_positions=len(paths),accepted=accepted))
        print('Analyzed W4',bi,flush=True)
    pd.DataFrame(coverage).to_csv(OUT/'w4_coverage.csv',index=False)
    for family,rows in [('dense',dense),('quantum',quantum)]:
        frame=pd.DataFrame(rows); frame.to_csv(OUT/f'w4_{family}_positions.csv',index=False)
        metrics={'energy_mae':'kcal/mol','energy_rmse':'kcal/mol','energy_nrmse':'1',
                 'centered_orientation_rmse':'kcal/mol','centered_orientation_nrmse':'1','torque_vrmse':'kcal/mol'}
        if family=='dense': metrics.update(esp_rmse_mEh='mEh/e',field_vrmse='kcal/mol/e/A',force_vrmse='kcal/mol/A')
        aggregate_cases(frame,'W4_'+family,metrics,('region','clearance_A'),
                        ('zero','mace_polar_l','exact_dipole' if family=='dense' else 'exact_dipole_base_com'))
        # Compare normalized advantage between near and far, not merely decaying absolute errors.
        for comparator in ['mace_polar_l','exact_dipole' if family=='dense' else 'exact_dipole_base_com']:
            for metric in ['energy_nrmse','centered_orientation_nrmse']:
                pivot=frame.groupby(['molecule_id','region','method'])[metric].mean().unstack('method')
                diff=(pivot['glider']-pivot[comparator]).unstack('region')
                record('mechanism_'+family,'near-minus-far',metric,'1','glider',diff['near'],
                       comparator,diff['far'],'Negative contrast means stronger normalized near-field advantage')
    # Angular modulation fraction is reference-only; each distance analyzed separately.
    df=pd.DataFrame(dense)
    for dist,part in df.groupby('clearance_A'):
        pivot=part.pivot(index=['molecule_id','direction'],columns='method',values='centered_orientation_nrmse')
        anis=part[part.method=='glider'].set_index(['molecule_id','direction'])['anisotropy_fraction'].reindex(pivot.index)
        mids=np.array(pivot.index.get_level_values(0)); unique=np.unique(mids)
        for comparator in ['mace_polar_l','exact_dipole']:
            gain=pivot[comparator]-pivot.glider
            rho=float(spearmanr(anis,gain).statistic)
            rng=np.random.default_rng(SEED); values=[]
            for _ in range(2000):
                idx=np.concatenate([np.flatnonzero(mids==m) for m in rng.choice(unique,len(unique),replace=True)])
                values.append(float(spearmanr(anis.iloc[idx],gain.iloc[idx]).statistic))
            lo,hi=np.quantile(values,[.025,.975])
            SUMMARY.append(dict(family='anisotropy',subset=f'clearance_A={dist}',metric='Spearman_QM_modulation_vs_normalized_gain',
                                unit='1',method='glider',comparator=comparator,n_systems=len(unique),estimate=rho,
                                ci_low=float(lo),ci_high=float(hi),notes='Exploratory; molecule-blocked, fixed-distance angular correlation'))
    return pd.DataFrame(dense),pd.DataFrame(quantum)


def qm_mm_reference():
    """Contract saved densities with AO Coulomb integrals. No SCF or new geometries."""
    from pyscf import df,gto,lib
    lib.num_threads(1)
    manifest=jread(C/'configs/embedding_manifest.json'); directory=OUT/'qm_mm_fields'; directory.mkdir(exist_ok=True)
    bookkeeping=[]
    for geom in manifest['geometries']:
        gid=geom['geometry_id']; vp=C/'qm'/f'B_{gid}_vacuum/result.npz'
        if not vp.exists(): continue
        vacuum=arrays(vp); mol=gto.M(atom=list(zip(geom['symbols'],geom['positions_angstrom'])),basis=geom['basis'],verbose=0)
        atoms=Atoms(geom['symbols'],positions=geom['positions_angstrom'])
        for mode in ['source','flipped','outward_1p5','outward_4p0']:
            selected=[r for r in manifest['jobs'] if r['geometry_id']==gid and r['mode']==mode]
            first=arrays(C/'qm'/selected[0]['job_id']/'result.npz'); points=first['mm_positions_angstrom']; h=.001
            cache=directory/f'{gid}_{mode}_integrals.npz'
            # Integral arrays are short-lived; the small derived fields are the durable cache.
            integral=df.incore.aux_e2(mol,gto.fakemol_for_charges(points/B))
            gradient=[]
            for axis in range(3):
                delta=np.eye(3)[axis]*h
                plus=df.incore.aux_e2(mol,gto.fakemol_for_charges((points+delta)/B))
                minus=df.incore.aux_e2(mol,gto.fakemol_for_charges((points-delta)/B))
                gradient.append((plus-minus)/(2*h))
            gradient=np.asarray(gradient)
            displacement=points[:,None,:]-atoms.positions[None,:,:]; radius=np.linalg.norm(displacement,axis=-1)
            nuclear_esp=(B/radius)@atoms.numbers
            nuclear_field=np.sum(B*atoms.numbers[None,:,None]*displacement/radius[:,:,None]**3,axis=1)
            vesp=nuclear_esp-np.einsum('ijp,ji->p',integral,vacuum['density_matrix'])
            vfield=nuclear_field+np.einsum('aijp,ji->pa',gradient,vacuum['density_matrix'])
            for row in selected:
                jid=row['job_id']; z=arrays(C/'qm'/jid/'result.npz')
                assert np.allclose(points,z['mm_positions_angstrom'])
                resp=-np.einsum('ijp,ji->p',integral,z['response_density_matrix'])
                field=np.einsum('aijp,ji->pa',gradient,z['response_density_matrix'])
                q=z['mm_charges_e']; frozen_vacuum=float(q@vesp*K)
                coupling=float(q@resp*K); relaxation=float(z['energy_difference_hartree']*K-frozen_vacuum)
                np.savez_compressed(directory/f'{jid}.npz',response_esp=resp,response_field=field,
                                    vacuum_esp=vesp,vacuum_field=vfield,mm_positions=points,mm_charges=q,
                                    response_coupling_kcal_mol=coupling,relaxation_energy_kcal_mol=relaxation,
                                    vacuum_coulomb_kcal_mol=frozen_vacuum)
                bookkeeping.append({**row,'response_coupling_kcal_mol':coupling,
                                    'relaxation_energy_kcal_mol':relaxation,'vacuum_coulomb_kcal_mol':frozen_vacuum,
                                    'embedded_minus_vacuum_kcal_mol':float(z['energy_difference_hartree'])*K,
                                    'half_response_coupling_kcal_mol':coupling/2,
                                    'half_coupling_error_kcal_mol':coupling/2-relaxation,
                                    'qm_mm_field_vrms':vector_rms(field)*K,
                                    'minimum_mm_atom_distance_A':float(radius.min())})
        print('Derived frozen-density MM fields',gid,flush=True)
    frame=pd.DataFrame(bookkeeping); frame.to_csv(OUT/'embedding_energy_bookkeeping.csv',index=False)
    return frame


def embedding():
    bookkeeping=qm_mm_reference()
    dataset=jread(C/'configs/training_dataset.json'); records=dataset['records']
    baseline=[]; trained=[]; validation=[]; linearity=[]
    runs=sorted((C/'results').glob('F_*'))
    test_ids={r['job_id'] for r in records if r['split']=='test'}
    for directory in runs:
        meta=jread(directory/'success.json'); archived=jread(directory/'test_metrics.json')
        assert {r['job_id'] for r in archived}==test_ids
        assert meta['n_test']==len(test_ids) and meta['epochs']==64
        validation.append(dict(run=directory.name,n_train=meta['n_train'],n_train_molecules=len(meta['training_molecules']),
                               nominal_fraction=meta['fraction'],actual_label_fraction=meta['n_train']/240,
                               actual_molecule_fraction=len(meta['training_molecules'])/12,
                               training_molecules=';'.join(meta['training_molecules']),
                               n_test=meta['n_test'],seed=meta['seed'],mode=meta['mode'],runtime_seconds=meta['runtime_seconds']))
    pd.DataFrame(validation).to_csv(OUT/'training_runs.csv',index=False)
    for row in records:
        jid=row['job_id']; z=arrays(C/'qm'/jid/'result.npz'); data=arrays(ROOT/row['file'])
        mm=arrays(OUT/'qm_mm_fields'/f'{jid}.npz'); mom=arrays(C/'results'/f'C_{row["geometry_id"]}'/f'{jid}.npz')
        origin=mom['origin_angstrom']; raw={'dipole':mom['response_raw_rank1'],'second':mom['response_raw_rank2'],'third':mom['response_raw_rank3']}
        ref=z['response_esp_hartree_per_e']; signal=rms(ref); refdip=z['induced_dipole_debye']
        mmref=mm['response_esp']; fieldref=mm['response_field']; coupling=float(mm['response_coupling_kcal_mol'])
        metadata={k:row[k] for k in ['job_id','geometry_id','molecule_id','split','mode','alpha']}
        basecommon={**metadata,'signal_mEh':signal*1000,'coupling_signal_kcal_mol':abs(coupling),
                    'field_signal':vector_rms(fieldref)*K}
        def metrics(method,pred,phi,field,dipole):
            return {**basecommon,'method':method,'esp_rmse_mEh':rms(pred-ref)*1000,
                    'esp_nrmse':rms(pred-ref)/max(signal,1e-12),
                    'dipole_error_D':float(np.linalg.norm(dipole-refdip)),
                    'mm_esp_rmse_mEh':rms(phi-mmref)*1000,
                    'mm_field_vrmse':vector_rms(field-fieldref)*K,
                    'mm_field_nrmse':vector_rms(field-fieldref)/max(vector_rms(fieldref),1e-12),
                    'coupling_error_kcal_mol':abs(float(mm['mm_charges']@phi*K)-coupling),
                    'coupling_prediction_kcal_mol':float(mm['mm_charges']@phi*K),
                    'coupling_reference_kcal_mol':coupling,
                    'half_coupling_relaxation_error_kcal_mol':abs(float(mm['mm_charges']@phi*K/2)-float(mm['relaxation_energy_kcal_mol']))}
        baseline.append(metrics('zero',np.zeros_like(ref),np.zeros_like(mmref),np.zeros_like(fieldref),np.zeros(3)))
        for order,label in [(1,'exact_dipole'),(2,'exact_dipole_quadrupole'),(3,'exact_dipole_quadrupole_octupole')]:
            pred=multipole(z['points_angstrom'],origin,raw,order); phi,field=multipole_field(mm['mm_positions'],origin,raw,order)
            baseline.append(metrics(label,pred,phi,field,refdip))
        q,u=data['response_q'],data['response_u']; pred=data['design']@np.r_[q,u.ravel()]
        phi,field=site_field(data['positions'],mm['mm_positions'],q,u)
        baseline.append(metrics('QM_fitted_local_sites_oracle',pred,phi,field,refdip))
        if row['split']=='test':
            for directory in runs:
                meta=jread(directory/'success.json'); p=finite_archive(directory/f'{jid}_prediction.npz')
                np.testing.assert_allclose(p['predicted_response_esp'],p['predicted_esp'] if meta['target']=='response' else p['predicted_esp']-z['vacuum_esp_hartree_per_e'],rtol=1e-5,atol=1e-7)
                assert abs(float(p['predicted_q'].sum()))<1e-5
                dipole=(np.sum(p['predicted_q'][:,None]*data['positions']/B,axis=0)+np.sum(p['predicted_u'],axis=0))*D
                np.testing.assert_allclose(dipole,p['predicted_dipole_debye'],rtol=2e-4,atol=2e-5)
                phi,field=site_field(data['positions'],mm['mm_positions'],p['predicted_q'],p['predicted_u'])
                if meta['target']=='total':
                    phi-=mm['vacuum_esp']; field-=mm['vacuum_field']; dipole-=z['dipole_debye']-refdip
                method=meta['mode'] if meta['fraction']==1 else directory.name
                result=metrics(method,p['predicted_response_esp'],phi,field,dipole)
                result.update(run=directory.name,seed=meta['seed'],fraction=meta['fraction'])
                archived=next(r for r in jread(directory/'test_metrics.json') if r['job_id']==jid)
                assert abs(result['esp_nrmse']-archived['response_nrmse'])<1e-7
                trained.append(result)
    base=pd.DataFrame(baseline); train=pd.DataFrame(trained)
    base.to_csv(OUT/'embedding_baseline_cases.csv',index=False); train.to_csv(OUT/'embedding_model_cases.csv',index=False)
    metrics={'esp_rmse_mEh':'mEh/e','esp_nrmse':'1','dipole_error_D':'D',
             'mm_esp_rmse_mEh':'mEh/e','mm_field_vrmse':'kcal/mol/e/A','mm_field_nrmse':'1',
             'coupling_error_kcal_mol':'kcal/mol','half_coupling_relaxation_error_kcal_mol':'kcal/mol'}
    aggregate_cases(base,'embedding_QM_compression',metrics,('split','mode','alpha'),('zero','exact_dipole'))
    combined=pd.concat([base[base.split=='test'],train],ignore_index=True)
    aggregate_cases(combined,'embedding_test',metrics,('mode','alpha'),('zero','exact_dipole','total_B','total_C'))
    # Seed-specific outcomes retained, while primary comparisons average seeds within molecule.
    seedframe=train.copy(); seedframe['method']=seedframe['run']
    aggregate_cases(seedframe,'embedding_test_each_run',metrics,(),())
    # Charge homogeneity: all four scales remain in the same block and split.
    for (gid,mode),part in bookkeeping.groupby(['geometry_id','mode']):
        row1=part[part.alpha==1].iloc[0]; basez=arrays(C/'qm'/row1.job_id/'result.npz')
        fields=[]; alphas=[]
        for row in part.sort_values('alpha').itertuples():
            z=arrays(C/'qm'/row.job_id/'result.npz'); fields.append(z['response_esp_hartree_per_e']); alphas.append(row.alpha)
        fields=np.asarray(fields); alphas=np.asarray(alphas)
        residual=fields-alphas[:,None]*basez['response_esp_hartree_per_e'][None,:]
        blind=fields-fields.mean(0,keepdims=True)
        linearity.append(dict(geometry_id=gid,molecule_id=row1.molecule_id,split=row1.split,mode=mode,
                              relative_nonlinearity=rms(residual)/max(rms(fields),1e-15),
                              alpha_blind_pooled_nrmse_lower_bound=rms(blind)/max(rms(fields),1e-15)))
    pd.DataFrame(linearity).to_csv(OUT/'charge_linearity.csv',index=False)
    frozen_inference_analysis()
    return combined


def frozen_inference_analysis():
    """Sanity/disagreement only: no target-mismatched accuracy leaderboard."""
    sys.path.insert(0,str(ROOT/'scripts/posthoc/downstream_response'))
    from common import gaussian_multipole_potential
    from glider.inference import esp_from_sites
    directory=C/'results/D_frozen_inference'; meta=jread(directory/'success.json')
    gd=ROOT/meta['glider_directory']; registry=pd.read_csv(gd/'prediction_registry.csv')
    frames=read(C/'data/compatible_explicit_waters.extxyz',':'); fmap={str(a.info['config_id']):a for a in frames}
    priors={name:finite_archive(directory/f'water_prior_{name}.npz') for name in ['m','l']}
    rows=[]
    for r in registry.itertuples():
        a=fmap[r.config_id]; pred=finite_archive(gd/r.prediction_file)
        assert hashlib.sha256((gd/r.prediction_file).read_bytes()).hexdigest()==r.prediction_sha256
        assert abs(float(pred['predicted_charges_e'].sum()))<1e-5
        phi=esp_from_sites(a.positions,pred['points_angstrom'],pred['predicted_charges_e'],pred['predicted_dipoles_e_bohr'])
        np.testing.assert_allclose(phi,pred['predicted_esp_hartree_per_e'],rtol=1e-5,atol=1e-7)
        for label,z in priors.items():
            ix=int(np.flatnonzero(z['config_ids']==r.config_id)[0]); first,last=z['offsets'][ix:ix+2]
            q,u=z['charges_e'][first:last],z['dipoles_e_bohr'][first:last]
            assert len(q)==len(a) and abs(float(q.sum()))<1e-5
            other=gaussian_multipole_potential(a.positions,pred['points_angstrom'],q,u)
            dip=(np.sum(q[:,None]*a.positions/B,axis=0)+u.sum(0))*D
            rows.append(dict(config_id=r.config_id,molecule_id=r.molecule_id,method='mace_polar_'+label,
                             shared_geometry_has_pointcharge_QM=r.molecule_id!='prospective_diverse_02',
                             glider_mace_field_disagreement_mEh=rms(phi-other)*1000,
                             glider_mace_dipole_disagreement_D=float(np.linalg.norm(pred['predicted_dipole_debye']-dip)),
                             mace_response_net_charge=float(q.sum()),glider_response_net_charge=float(pred['predicted_charges_e'].sum())))
    assert len(registry)==96
    pd.DataFrame(rows).to_csv(OUT/'frozen_model_disagreement_NOT_accuracy.csv',index=False)


def finish():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    old=pd.concat([pd.read_csv(OUT/'summary_initial.csv'),pd.read_csv(OUT/'summary_embedding.csv')],ignore_index=True)
    base=pd.read_csv(OUT/'embedding_baseline_cases.csv'); models=pd.read_csv(OUT/'embedding_model_cases.csv')
    combined=pd.concat([base[base.split=='test'],models],ignore_index=True)
    metrics={'esp_rmse_mEh':'mEh/e','esp_nrmse':'1','mm_field_vrmse':'kcal/mol/e/A',
             'mm_field_nrmse':'1','coupling_error_kcal_mol':'kcal/mol'}
    aggregate_cases(combined,'embedding_extra_comparisons',metrics,(),('exact_dipole_quadrupole','exact_dipole_quadrupole_octupole','QM_fitted_local_sites_oracle'))
    dense=pd.read_csv(OUT/'w4_dense_positions.csv')
    aggregate_cases(dense,'W4_discrete_orientation',{'preferred_orientation_accuracy':'fraction'},('region','clearance_A'),('zero','mace_polar_l','exact_dipole'))
    # The independent unit remains the molecule, including charge-scale consistency checks.
    curves=[]
    dataset=jread(C/'configs/training_dataset.json')['records']
    test=[r for r in dataset if r['split']=='test']
    for directory in sorted((C/'results').glob('F_*')):
        meta=jread(directory/'success.json'); method=meta['mode'] if meta['fraction']==1 else directory.name
        for gid,mode in sorted({(r['geometry_id'],r['mode']) for r in test}):
            records=sorted([r for r in test if r['geometry_id']==gid and r['mode']==mode],key=lambda r:r['alpha'])
            pred=np.asarray([arrays(directory/f'{r["job_id"]}_prediction.npz')['predicted_response_esp'] for r in records])
            ref=np.asarray([arrays(C/'qm'/r['job_id']/'result.npz')['response_esp_hartree_per_e'] for r in records])
            alpha=np.array([r['alpha'] for r in records]); one=int(np.flatnonzero(alpha==1)[0])
            curves.append(dict(method=method,run=directory.name,molecule_id=records[0]['molecule_id'],geometry_id=gid,mode=mode,
                               response_homogeneity_error_relative_to_QM=rms(pred-alpha[:,None]*pred[one])/rms(ref),
                               QM_homogeneity_error=rms(ref-alpha[:,None]*ref[one])/rms(ref)))
    pd.DataFrame(curves).to_csv(OUT/'model_charge_homogeneity.csv',index=False)
    aggregate_cases(pd.DataFrame(curves),'charge_homogeneity',{'response_homogeneity_error_relative_to_QM':'1'},(),())
    null=pd.read_csv(OUT/'null_field_inference.csv')
    aggregate_cases(null,'null_field',{'null_response_rms_mEh':'mEh/e'},(),())
    replay=pd.read_csv(OUT/'checkpoint_replay_metrics.csv')
    aggregate_cases(replay,'checkpoint_replay',{'esp_nrmse':'1','own_target_nrmse':'1'},('partition',),())
    summary=pd.concat([old,pd.DataFrame(SUMMARY)],ignore_index=True)
    summary.to_csv(C/'preliminary_results_table.csv',index=False)
    system=combined.groupby(['method','molecule_id'])[list(metrics)].mean().reset_index()
    system.to_csv(OUT/'embedding_per_molecule.csv',index=False)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'savefig.dpi':180,'figure.dpi':120,'axes.grid':True,'grid.alpha':.18})
    labels={'glider':'GLIDER','mace_polar_l':'MACE-POLAR-L','zero':'Zero response','exact_dipole':'Exact dipole',
            'exact_dipole_quadrupole':'Exact D+Q','exact_dipole_quadrupole_octupole':'Exact D+Q+O',
            'response_C':'GLIDER-EE','total_B':'Direct total B','total_C':'Direct total C',
            'QM_fitted_local_sites_oracle':'QM-fitted sites (oracle)'}
    colors={'glider':'#007a8a','mace_polar_l':'#da7c30','zero':'#777777','exact_dipole':'#8a508f',
            'exact_dipole_quadrupole':'#cc6677','exact_dipole_quadrupole_octupole':'#4477aa',
            'response_C':'#007a8a','total_B':'#da7c30','total_C':'#aa4499','QM_fitted_local_sites_oracle':'#228833'}
    def save(fig,name):
        fig.savefig(FIG/(name+'.png'),bbox_inches='tight'); fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight'); plt.close(fig)
    def get(family,subset,metric,method,comparator=None):
        filt=(summary.family==family)&(summary['subset']==subset)&(summary.metric==metric)&(summary.method==method)
        filt &= summary.comparator.isna() if comparator is None else summary.comparator==comparator
        got=summary[filt]; assert len(got)==1,(family,subset,metric,method,comparator,len(got))
        return got.iloc[0]
    def distance_curves(ax,family,metric,methods,ylabel,log=False):
        for method in methods:
            part=summary[(summary.family==family)&(summary.metric==metric)&(summary.method==method)&summary.comparator.isna()&summary['subset'].str.startswith('clearance_A=')].copy()
            part['x']=part['subset'].str.split('=').str[1].astype(float); part=part.sort_values('x')
            ax.plot(part.x,part.estimate,'o-',markersize=4,label=labels.get(method,method),color=colors.get(method))
            if method in ['glider','mace_polar_l']: ax.fill_between(part.x,part.ci_low,part.ci_high,alpha=.12,color=colors[method])
        if log: ax.set_yscale('log')
        ax.set_xlabel('Probe oxygen clearance from base vdW surface (Å)'); ax.set_ylabel(ylabel)
    fig,axs=plt.subplots(2,2,figsize=(12,8.5),constrained_layout=True)
    distance_curves(axs[0,0],'W4_dense','energy_mae',['glider','mace_polar_l','zero'],'Mean energy MAE (kcal/mol)',True)
    axs[0,0].set_title('a  TIP3P probe: modest advantage over MACE'); axs[0,0].legend(fontsize=8)
    distance_curves(axs[0,1],'W4_dense','energy_nrmse',['glider','exact_dipole','exact_dipole_quadrupole','exact_dipole_quadrupole_octupole'],'Mean profile energy NRMSE',True)
    axs[0,1].axhline(1,color='#777777',ls=':',lw=1); axs[0,1].set_title('b  Higher moments recover in the exterior regime'); axs[0,1].legend(fontsize=8)
    distance_curves(axs[1,0],'W4_quantum','energy_mae',['glider','mace_polar_l','zero'],'Mean energy MAE (kcal/mol)',True)
    axs[1,0].set_title('c  Full quantum-water probe (same ten base systems)')
    for family,color,label in [('W4_dense','#007a8a','TIP3P probe'),('W4_quantum','#da7c30','Quantum-water probe')]:
        part=summary[(summary.family==family)&(summary.method=='glider')&(summary.comparator=='mace_polar_l')&(summary.metric=='energy_nrmse')&summary['subset'].str.startswith('clearance_A=')].copy()
        part['x']=part['subset'].str.split('=').str[1].astype(float); part=part.sort_values('x')
        axs[1,1].errorbar(part.x,part.paired_delta,yerr=[part.paired_delta-part.delta_ci_low,part.delta_ci_high-part.paired_delta],fmt='o-',color=color,capsize=3,label=label)
    axs[1,1].axhline(0,color='black',lw=1); axs[1,1].set_xlabel('Probe oxygen clearance (Å)'); axs[1,1].set_ylabel('GLIDER − MACE energy NRMSE'); axs[1,1].set_title('d  No consistent normalized near-field gain over MACE'); axs[1,1].legend(fontsize=8)
    save(fig,'01_w4_distance_and_multipoles')
    fig,axs=plt.subplots(2,2,figsize=(12,8.5),constrained_layout=True)
    mids=['MNSOL-CONN-XSQUKJJJFZCRTK','MNSOL-CONN-BYEAHWXPCBROCE','pros_azaaryl_sulfoxide']
    names=['Urea','Hexafluoroisopropanol','Fluorinated azaaryl\nsulfoxide']; x=np.arange(3)
    for ix,method in enumerate(['response_C','total_B','total_C','exact_dipole']):
        values=system[system.method==method].set_index('molecule_id').loc[mids,'esp_nrmse']
        axs[0,0].bar(x+(ix-1.5)*.18,values,width=.17,color=colors[method],label=labels[method])
    axs[0,0].axhline(1,color='black',ls=':',lw=1); axs[0,0].set_yscale('log'); axs[0,0].set_xticks(x,names); axs[0,0].set_ylabel('Response-ESP NRMSE (log scale)'); axs[0,0].set_title('a  Three held-out molecules; no large validation panel'); axs[0,0].legend(fontsize=8)
    ordered=['zero','exact_dipole','exact_dipole_quadrupole_octupole','QM_fitted_local_sites_oracle','response_C','total_B','total_C']
    for ax,metric,ylabel,title in [(axs[0,1],'coupling_error_kcal_mol','Frozen response-coupling MAE (kcal/mol)','b  Surface-ESP gains do not guarantee coupling superiority'),
                                    (axs[1,0],'mm_field_nrmse','MM-site field NRMSE','c  Exact D+Q+O is a stronger field control')]:
        for ix,method in enumerate(ordered):
            row=get('embedding_test','all',metric,method)
            ax.errorbar(ix,row.estimate,yerr=[[row.estimate-row.ci_low],[row.ci_high-row.estimate]],fmt='o',color=colors.get(method),capsize=4)
        ax.set_xticks(range(len(ordered)),[labels[k].replace(' ','\n',1) for k in ordered],fontsize=8)
        if metric=='mm_field_nrmse': ax.set_yscale('log'); ax.axhline(1,color='black',ls=':',lw=1)
        ax.set_ylabel(ylabel); ax.set_title(title)
    modes=['source','flipped','outward_1p5','outward_4p0']
    for method in ['response_C','exact_dipole','zero']:
        values=[get('embedding_test','mode='+m,'esp_nrmse',method).estimate for m in modes]
        axs[1,1].plot(range(4),values,'o-',color=colors[method],label=labels[method])
    axs[1,1].set_xticks(range(4),['Source','Flipped','Outward 1.5 Å','Outward 4 Å']); axs[1,1].set_ylabel('Response-ESP NRMSE'); axs[1,1].set_title('d  Weak-field cases expose a failure mode'); axs[1,1].legend(fontsize=8)
    save(fig,'02_embedding_controls_and_physics')
    fig,axs=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    order=['F_efficiency_05','F_efficiency_10','F_efficiency_25','F_efficiency_50','response_C']; subset_x=np.arange(5)
    for ax,metric,label,title in [(axs[0],'esp_nrmse','Response-ESP NRMSE','a  More labels do not monotonically improve NRMSE'),
                                  (axs[1],'coupling_error_kcal_mol','Frozen coupling MAE (kcal/mol)','b  Coupling error does improve with more labels')]:
        rows=[get('embedding_test','all',metric,method) for method in order]
        mean=np.array([r.estimate for r in rows]); low=np.array([r.ci_low for r in rows]); high=np.array([r.ci_high for r in rows])
        ax.errorbar(subset_x,mean,yerr=[mean-low,high-mean],fmt='o-',color='#007a8a',capsize=4)
        ax.axhline(get('embedding_test','all',metric,'zero').estimate,color='black',ls=':',label='Zero response')
        ax.set_xticks(subset_x,['1\n32 labels','2\n64 labels','3\n96 labels','6\n144 labels','12\n240 labels']); ax.set_xlabel('Training subset: molecules / response labels'); ax.set_ylabel(label); ax.set_title(title)
        for seed in [2026091401,2026091402]:
            r=get('embedding_test_each_run','all',metric,f'F_response_C_{seed}')
            ax.plot(4,r.estimate,'x',color='#aa4499')
        ax.legend(fontsize=8)
    save(fig,'03_data_efficiency')
    print('Wrote result table and three figures',flush=True)


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--phase',choices=['initial','embedding','finish'],default='initial'); args=parser.parse_args()
    OUT.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)
    if args.phase=='initial':
        jsave(OUT/'analysis_plan.json',{'created_utc':pd.Timestamp.now(tz='UTC').isoformat(),'post_hoc':True,
            'bins_angstrom':{'near':[.8,1.2,1.8,2.5],'intermediate':[3.5,5.],'far':[7.5,10.,15.]},
            'bootstrap_blocks':'molecule; all configurations, geometries and seeds stay together',
            'bootstrap_replicates':NBOOT,'seed':SEED,'small_n':'exact empirical bootstrap enumeration at n<=4; report exact two-sided sign test; cannot establish population reliability from n=3',
            'primary_embedding_partition':'fixed held-out test only; no model/subset selected by test performance',
            'all_existing_models_included':True,'energy_boundary':'frozen Coulomb response coupling is distinct from variational polarization energy and total interaction energy'})
        audit(); w4()
        pd.DataFrame(SUMMARY).to_csv(OUT/'summary_initial.csv',index=False)
    elif args.phase=='embedding':
        embedding()
        pd.DataFrame(SUMMARY).to_csv(OUT/'summary_embedding.csv',index=False)
    elif args.phase=='finish':
        finish()


if __name__=='__main__':
    sys.path.insert(0,str(ROOT))
    main()
