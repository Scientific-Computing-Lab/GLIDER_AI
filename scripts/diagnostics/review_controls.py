#!/usr/bin/env python3
"""Post hoc controls requested during author review on 4 October 2026.

The published GLIDER checkpoint is never changed. Training uses only the original
48 examples and the original five seeds, 64 epochs, objective and optimizer.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import random
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import torch
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[2]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def frozen(root):
    options = [root / 'provenance/frozen_source/prospective_1',
               root / 'evidence/frozen_source/prospective_1']
    directory = next(p for p in options if p.exists())
    spec = importlib.util.spec_from_file_location('review_frozen_head', directory / 'response_learning.py')
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module, directory


def paths(root):
    new = root / 'experiments/training'
    if new.exists():
        return new / 'geometries/configurations.extxyz', new / 'references', new / 'features'
    old = root / 'data/development'
    return old / 'combined_development.extxyz', old / 'observables', old / 'frozen_training'


def ablated_class(module):
    class WithoutGlobal(module.ResponseHead):
        def forward(self, sample):
            # The B forward branch is precisely the shared local mapping, before
            # reconciliation and prior addition. Restore the C loss weight and
            # add the same prior explicitly. All parameters were initialized as C.
            self.head = 'B'
            try:
                q, u = super().forward(sample)
            finally:
                self.head = 'C'
            q = q + sample.base_q
            return q - q.mean(), u + sample.base_u
    return WithoutGlobal


def train(args):
    m, directory = frozen(args.source_root)
    geometry, refs, feature_dir = paths(args.source_root)
    feature = feature_dir / 'mace_polar_m_features.npz'
    prior = feature_dir / 'mace_polar_ml_average.npz'
    samples = m.load_samples(SimpleNamespace(configurations=[geometry], observables=[refs],
        features=feature, base_predictions=prior, encoder='mace-polar-m'))
    args.output.mkdir(parents=True, exist_ok=False)
    protocol = dict(date='2026-10-04', status='post hoc author-review diagnostic',
        n_training_configurations=48, n_training_solutes=14, epochs=64,
        seeds=list(range(2026081300, 2026081305)), learning_rate=0.002,
        weight_decay=0.00002, dipole_loss_weight=2.0, latent_loss_weight=0.05,
        selection='original fixed training schedule; no test-based model selection',
        intervention='remove global dipole readout and moment reconciliation; retain local mapping and averaged prior',
        source_sha256=digest(directory/'response_learning.py'),
        feature_sha256=digest(feature), prior_sha256=digest(prior),
        published_checkpoint_sha256=digest(args.source_root/'checkpoints/glider_site_response_ensemble.pt'),
        torch_version=torch.__version__, device=args.device)
    (args.output/'protocol.json').write_text(json.dumps(protocol, indent=2)+'\n')
    device = torch.device(args.device)
    data = [m.move(s, device) for s in samples]
    rows = []
    for variant, cls in [('matched_full', m.ResponseHead), ('without_global', ablated_class(m))]:
        states = []
        for seed in protocol['seeds']:
            random.seed(seed)
            np.random.seed(seed)
            torch.manual_seed(seed)
            model = cls(data, 'C').to(device)
            optimizer = torch.optim.AdamW(model.parameters(), lr=0.002, weight_decay=0.00002)
            scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, 64, eta_min=0.00002)
            for epoch in range(64):
                model.train()
                optimizer.zero_grad(set_to_none=True)
                loss = m.loss_for(model, data, True)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
                optimizer.step()
                scheduler.step()
            model.eval()
            with torch.no_grad():
                final_loss = float(m.loss_for(model, data, False))
            states.append({k:v.detach().cpu() for k,v in model.state_dict().items()})
            rows.append(dict(variant=variant, seed=seed, loss=final_loss))
            print(variant, seed, final_loss, flush=True)
        torch.save(dict(format='posthoc-global-branch-control-v1', head='C',
            variant=variant, epochs=64, seeds=protocol['seeds'], states=states), args.output/f'{variant}.pt')
    pd.DataFrame(rows).to_csv(args.output/'training_losses.csv', index=False)


def build_dissociation(args):
    m, _ = frozen(args.source_root)
    geometry, _, features = paths(args.source_root)
    ids = set(np.load(features/'mace_polar_m_features.npz')['config_ids'].astype(str))
    frames = [a for a in read(geometry, index=':') if str(a.info['config_id']) in ids]
    parents = {}
    for a in sorted(frames, key=lambda x: str(x.info['config_id'])):
        if a.info.get('regime',a.info.get('motif')) in ['equilibrium','compact']:
            parents[str(a.info['molecule_id'])] = a
    if len(parents) != 14:
        raise RuntimeError(f'Expected one compact/equilibrium parent per training solute, got {len(parents)}')
    args.output.mkdir(parents=True, exist_ok=False)
    geometries, rows, probes = [], [], {}
    # Radial translation away from the solute, with internal geometry fixed.
    # Distance is the minimum solute/environment atom distance, not probe clearance.
    distances = args.distances
    for key, parent in parents.items():
        ns = int(parent.info['n_solute_atoms'])
        solute = parent[:ns]
        centre = solute.positions.mean(0)
        from glider.inference import molecular_surface
        probes[key] = molecular_surface(solute.get_chemical_symbols(), solute.positions)
        for system in ['single_water','whole_environment']:
            end = ns+3 if system == 'single_water' else len(parent)
            base = parent[:end]
            environment = base.positions[ns:].copy()
            direction = environment[0] - centre
            direction /= np.linalg.norm(direction)
            def minimum(t):
                return np.linalg.norm((environment+t*direction)[:,None]-solute.positions[None],axis=2).min()
            for distance in distances:
                low, high = 0., max(100.,2*distance)
                while minimum(high)<distance:
                    high*=2
                if minimum(0) > distance:
                    shift=0.
                else:
                    for _ in range(70):
                        mid=(low+high)/2
                        if minimum(mid) < distance: low=mid
                        else: high=mid
                    shift=high
                a=base.copy()
                a.positions[ns:] = environment + shift*direction
                cid=f'dissociation__{key}__{system}__{distance:g}A'
                a.info=dict(config_id=cid,molecule_id=key,regime=system,n_solute_atoms=ns)
                geometries.append(a)
                rows.append(dict(config_id=cid,molecule_id=key,regime=system,
                    n_solute_atoms=ns,n_waters=(len(a)-ns)//3,
                    minimum_fragment_distance_A=float(minimum(shift)),requested_distance_A=distance))
    write(args.output/'configurations.extxyz',geometries)
    pd.DataFrame(rows).to_csv(args.output/'configuration_registry.csv',index=False)
    np.savez_compressed(args.output/'solute_probe_points.npz',**probes)
    (args.output/'protocol.json').write_text(json.dumps(dict(date='2026-10-04',
        status='post hoc dissociation diagnostic, no new QM labels',
        parents='all 14 original training solutes, compact/equilibrium configuration',
        systems=['solute plus first water','solute plus intact original environment'],
        separation='minimum interfragment atom distance in angstrom',
        observation='fixed surface anchored to stationary solute',
        distances_A=distances, n_configurations=len(rows),
        expected_limit='complex-minus-fragments response approaches zero as the fragments separate',
        checkpoints='unchanged original GLIDER and frozen M/L prior'),indent=2)+'\n')
    print('Dissociation configurations:',len(rows),flush=True)


def score_dissociation(args):
    m,_=frozen(args.source_root)
    root=args.output
    frames={str(a.info['config_id']):a for a in read(root/'configurations.extxyz',index=':')}
    registry=pd.read_csv(root/'configuration_registry.csv').set_index('config_id')
    points=np.load(root/'solute_probe_points.npz')
    feature=np.load(root/'mace_polar_m_features.npz')
    prior=np.load(root/'mace_polar_ml_average.npz')
    if not np.array_equal(feature['config_ids'],prior['config_ids']) or not np.array_equal(feature['offsets'],prior['offsets']):
        raise RuntimeError('Feature/prior IDs or atom offsets differ')
    checkpoint=torch.load(args.source_root/'checkpoints/glider_site_response_ensemble.pt',map_location='cpu',weights_only=False)
    models=[m.ResponseHead.from_state_dict(s,'C').eval() for s in checkpoint['states']]
    records=[]; arrays={}
    for i,cid in enumerate(feature['config_ids'].astype(str)):
        a=frames[cid];lo,hi=feature['offsets'][i:i+2];lo,hi=int(lo),int(hi)
        prod=feature['product_1'][lo:hi]; channels=prod.shape[1]//4
        sample=SimpleNamespace(scalars=torch.tensor(prod[:,:channels],dtype=torch.float32),
            vectors=torch.tensor(prod[:,channels:].reshape(len(a),channels,3)[:,:,[2,0,1]],dtype=torch.float32),
            geometry_vectors=torch.tensor(m.geometry_basis(a.positions,a.numbers)),
            positions=torch.tensor(a.positions,dtype=torch.float32),
            base_q=torch.tensor(prior['charges_e'][lo:hi],dtype=torch.float32),
            base_u=torch.tensor(prior['dipoles_e_bohr'][lo:hi],dtype=torch.float32))
        with torch.no_grad(): pred=[model(sample) for model in models]
        q=torch.stack([x[0] for x in pred]).mean(0).numpy();q-=q.mean()
        u=torch.stack([x[1] for x in pred]).mean(0).numpy()
        key=str(a.info['molecule_id']); grid=points[key]; design=m.esp_design(a.positions,grid)
        for method,charge,dip in [('glider',q,u),('averaged_prior',sample.base_q.numpy(),sample.base_u.numpy())]:
            esp=design@np.r_[charge,dip.ravel()]
            mu=(np.sum(charge[:,None]*a.positions/m.BOHR_TO_ANGSTROM,axis=0)+dip.sum(0))*m.DEBYE_PER_E_BOHR
            records.append(dict(config_id=cid,molecule_id=key,system=a.info['regime'],method=method,
                distance_A=float(registry.loc[cid,'minimum_fragment_distance_A']),
                esp_rms_mEh_per_e=float(np.sqrt(np.mean(esp**2))*1000),
                dipole_norm_D=float(np.linalg.norm(mu)),net_charge_e=float(charge.sum()),
                solute_response_charge_e=float(charge[:int(a.info['n_solute_atoms'])].sum())))
            arrays[f'{cid}__{method}']=esp
            arrays[f'{cid}__{method}__charges_e']=charge
            arrays[f'{cid}__{method}__site_dipoles_e_bohr']=dip
    pd.DataFrame(records).to_csv(root/'results.csv',index=False)
    np.savez_compressed(root/'predicted_probe_potentials.npz',**arrays)
    summary=pd.DataFrame(records).groupby(['system','method', 'distance_A'],as_index=False).mean(numeric_only=True)
    # Requested grid avoids tiny bisection rounding differences in grouping.
    df=pd.DataFrame(records);df['distance_A']=df.distance_A.round(6)
    df.groupby(['system','method','distance_A'],as_index=False)[['esp_rms_mEh_per_e','dipole_norm_D']].mean().to_csv(root/'summary.csv',index=False)
    print(df[df.distance_A==df.distance_A.max()].groupby(['system','method'])[['esp_rms_mEh_per_e','dipole_norm_D']].mean(),flush=True)


def evaluation_panel(root):
    frames=[]; refs={}; panels={}
    for n in [1,2,3]:
        new=root/f'experiments/panel_{n}'
        old=root/(f'data/prospective_{n}' if n<3 else 'canonical_panel_3')
        base=new if new.exists() else old
        geometry=(base/'geometries/configurations.extxyz' if new.exists()
                  else base/('configurations.extxyz' if n<3 else 'data/configurations.extxyz'))
        for a in read(geometry,index=':'):
            frames.append(a); panels[str(a.info['config_id'])]=n
        registry=pd.read_csv(base/'references/observable_registry.csv')
        for row in registry.itertuples(): refs[str(row.config_id)]=base/'references'/row.observable_file
    return frames,refs,panels


def build_evaluation(args):
    frames,_,_=evaluation_panel(args.source_root)
    args.output.mkdir(parents=True,exist_ok=True)
    write(args.output/'evaluation_configurations.extxyz',frames)
    print('Evaluation configurations:',len(frames),flush=True)


def score_ablation(args):
    m,_=frozen(args.source_root)
    frames,refs,panels=evaluation_panel(args.source_root)
    frame_map={str(a.info['config_id']):a for a in frames}
    feature=np.load(args.output/'evaluation_features.npz')
    prior=np.load(args.output/'evaluation_prior.npz')
    if not np.array_equal(feature['config_ids'],prior['config_ids']): raise RuntimeError('IDs differ')
    ensembles={}
    for variant,checkpoint,cls in [
        ('published_glider',args.source_root/'checkpoints/glider_site_response_ensemble.pt',m.ResponseHead),
        ('matched_full',args.output/'matched_full.pt',m.ResponseHead),
        ('without_global',args.output/'without_global.pt',ablated_class(m))]:
        payload=torch.load(checkpoint,map_location='cpu',weights_only=False)
        ensembles[variant]=[cls.from_state_dict(state,'C').eval() for state in payload['states']]
    pred_dir=args.output/'predictions';pred_dir.mkdir(exist_ok=True)
    records=[]; max_difference=0.
    for i,cid in enumerate(feature['config_ids'].astype(str)):
        a=frame_map[cid];lo,hi=map(int,feature['offsets'][i:i+2]);prod=feature['product_1'][lo:hi]
        channels=prod.shape[1]//4
        sample=SimpleNamespace(scalars=torch.tensor(prod[:,:channels],dtype=torch.float32),
            vectors=torch.tensor(prod[:,channels:].reshape(len(a),channels,3)[:,:,[2,0,1]],dtype=torch.float32),
            geometry_vectors=torch.tensor(m.geometry_basis(a.positions,a.numbers)),
            positions=torch.tensor(a.positions,dtype=torch.float32),
            base_q=torch.tensor(prior['charges_e'][lo:hi],dtype=torch.float32),
            base_u=torch.tensor(prior['dipoles_e_bohr'][lo:hi],dtype=torch.float32))
        ref=np.load(refs[cid]); points=ref['points_angstrom'];target=ref['delta_esp_hartree_per_e']
        design=m.esp_design(a.positions,points);payload={'points_angstrom':points}
        values={}
        for variant,models in ensembles.items():
            with torch.no_grad(): prediction=[model(sample) for model in models]
            q=torch.stack([x[0] for x in prediction]).mean(0).numpy();q-=q.mean()
            u=torch.stack([x[1] for x in prediction]).mean(0).numpy()
            esp=design@np.r_[q,u.ravel()]
            mu=(np.sum(q[:,None]*a.positions/m.BOHR_TO_ANGSTROM,axis=0)+u.sum(0))*m.DEBYE_PER_E_BOHR
            values[variant]=esp
            payload[variant+'_esp_hartree_per_e']=esp
            payload[variant+'_dipole_debye']=mu
            records.append(dict(panel=panels[cid],config_id=cid,molecule_id=str(a.info['molecule_id']),
                method=variant,esp_nrmse=float(np.sqrt(np.mean((esp-target)**2)/np.mean(target**2))),
                esp_rmse_mEh_per_e=float(np.sqrt(np.mean((esp-target)**2))*1000),
                dipole_mse_D2=float(np.mean((mu-ref['delta_dipole_debye'])**2))))
        max_difference=max(max_difference,float(np.max(np.abs(values['matched_full']-values['published_glider']))))
        np.savez_compressed(pred_dir/(hashlib.sha256(cid.encode()).hexdigest()[:20]+'.npz'),**payload)
    df=pd.DataFrame(records);df.to_csv(args.output/'results.csv',index=False)
    solute=df.groupby(['panel','molecule_id','method'],as_index=False)[['esp_nrmse','esp_rmse_mEh_per_e','dipole_mse_D2']].mean()
    solute.to_csv(args.output/'solute_results.csv',index=False)
    summary=solute.groupby(['panel','method'],as_index=False)[['esp_nrmse','esp_rmse_mEh_per_e','dipole_mse_D2']].mean()
    summary['dipole_rmse_D']=np.sqrt(summary.pop('dipole_mse_D2'))
    summary.to_csv(args.output/'summary.csv',index=False)
    pair=solute.pivot(index=['panel','molecule_id'],columns='method',values='esp_nrmse')
    delta=(pair.without_global-pair.published_glider).to_numpy()
    rng=np.random.default_rng(20261004)
    bootstrap=delta[rng.integers(0,len(delta),size=(100000,len(delta)))].mean(1)
    report=dict(n_test_solutes=len(delta),n_test_configurations=len(frames),
        max_retrained_full_vs_published_esp_difference_hartree_per_e=max_difference,
        esp_delta_without_global_minus_published=float(delta.mean()),
        paired_solute_bootstrap_95ci=np.quantile(bootstrap,[.025,.975]).tolist(),
        n_solutes_without_global_better=int((delta<0).sum()),
        scope='post hoc fixed-protocol intervention, no new prospective claim')
    (args.output/'evaluation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(summary.to_string(index=False),flush=True);print(json.dumps(report,indent=2),flush=True)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('operation',choices=['train','build-dissociation','score-dissociation','build-evaluation','score-ablation'])
    p.add_argument('--source-root',type=Path,default=ROOT)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--device',default='cuda')
    p.add_argument('--distances',type=float,nargs='+',default=[3.,4.,5.,6.,8.,10.,15.,20.])
    a=p.parse_args()
    {'train':train,'build-dissociation':build_dissociation,'score-dissociation':score_dissociation,
     'build-evaluation':build_evaluation,'score-ablation':score_ablation}[a.operation](a)


if __name__=='__main__':
    main()
