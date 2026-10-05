#!/usr/bin/env python3
"""Recompute response errors directly from released reference/prediction arrays."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
EXPERIMENTS=('panel_1','panel_2','panel_3','nonwater','liquid','shell_size')


def score(name, output):
    base=ROOT/'experiments'/name
    registry=next((base/'references').glob('*REGISTRY.csv'),None)
    if registry is None: registry=base/'references/observable_registry.csv'
    references=pd.read_csv(registry)
    geometry=pd.read_csv(base/'geometries/configuration_registry.csv').set_index('config_id')
    records=[]
    for directory in sorted((base/'predictions').iterdir()):
        if not directory.is_dir() or directory.name=='intermediate':continue
        pred_registry=directory/'prediction_registry.csv'
        prediction_map={}
        if pred_registry.exists():
            for row in pd.read_csv(pred_registry).itertuples():
                prediction_map[str(row.config_id)]=directory/str(row.prediction_file)
        for row in references.itertuples():
            cid=str(row.config_id)
            path=prediction_map.get(cid,directory/row.observable_file)
            if not path.exists():
                if list(directory.glob('*.npz')): raise FileNotFoundError(path)
                continue
            reference=np.load(base/'references'/row.observable_file)
            prediction=np.load(path)
            dipole_error=prediction['predicted_dipole_debye']-reference['delta_dipole_debye']
            record=dict(config_id=cid,molecule_id=str(row.molecule_id),method=directory.name,
                dipole_mse_D2=float(np.mean(dipole_error**2)))
            if 'predicted_esp_hartree_per_e' in prediction:
                ref=reference['delta_esp_hartree_per_e'];pred=prediction['predicted_esp_hartree_per_e']
                if 'points_angstrom' in prediction:
                    np.testing.assert_allclose(prediction['points_angstrom'],reference['points_angstrom'],atol=1e-9,rtol=0)
                if pred.shape!=ref.shape:raise RuntimeError(f'Probe shape mismatch: {path}')
                record.update(esp_nrmse=float(np.sqrt(np.mean((pred-ref)**2)/np.mean(ref**2))),
                    esp_rmse_mEh_per_e=float(np.sqrt(np.mean((pred-ref)**2))*1000),
                    response_rms_mEh_per_e=float(np.sqrt(np.mean(ref**2))*1000))
            for key in ['regime','n_waters','perturbant','parent_snapshot_index']:
                if key in geometry.columns:record[key]=geometry.loc[cid,key]
            records.append(record)
    df=pd.DataFrame(records)
    numeric=['esp_nrmse','esp_rmse_mEh_per_e','response_rms_mEh_per_e','dipole_mse_D2']
    expected=base/'results.csv'
    if expected.exists():
        archived=pd.read_csv(expected).sort_values(['config_id','method']).reset_index(drop=True)
        current=df.sort_values(['config_id','method']).reset_index(drop=True)
        pd.testing.assert_frame_equal(current,archived,check_dtype=False,atol=1e-12,rtol=1e-12)
    solute=df.groupby(['method','molecule_id'],as_index=False)[numeric].mean()
    summary=solute.groupby('method',as_index=False)[numeric].mean()
    summary['dipole_rmse_D']=np.sqrt(summary.pop('dipole_mse_D2'))
    output.mkdir(parents=True,exist_ok=True)
    df.to_csv(output/'results.csv',index=False)
    solute.to_csv(output/'solute_results.csv',index=False)
    summary.to_csv(output/'summary.csv',index=False)
    print(name, 'configurations',len(references),'solutes',references.molecule_id.nunique(),flush=True)
    print(summary[summary.method.isin(['glider','mace_polar_l','mace_polar_m'])].to_string(index=False),flush=True)
    return dict(experiment=name,configurations=len(references),solutes=int(references.molecule_id.nunique()),methods=summary.method.tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--experiment',choices=[*EXPERIMENTS,'all'],default='all')
    p.add_argument('--output',type=Path,default=ROOT/'build/recomputed')
    a=p.parse_args()
    names=EXPERIMENTS if a.experiment=='all' else [a.experiment]
    reports=[score(n,a.output/n) for n in names]
    (a.output/'scoring_report.json').write_text(json.dumps(reports,indent=2)+'\n')
