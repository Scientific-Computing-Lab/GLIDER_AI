#!/usr/bin/env python3
"""Verify training membership, raw prediction constraints and release hashes."""
from pathlib import Path
import csv,hashlib,json,shlex
import numpy as np
import pandas as pd
from glider.inference import dipole_from_sites
ROOT=Path(__file__).resolve().parents[2]

def frames(path):
    with path.open() as f:
        while line:=f.readline():
            n=int(line);info=dict(x.split('=',1) for x in shlex.split(f.readline()) if '=' in x)
            rows=[f.readline().split() for _ in range(n)]
            yield info,np.array([[float(v) for v in row[1:4]] for row in rows])

def main():
    manifest=json.loads((ROOT/'RELEASE_MANIFEST.json').read_text())
    for name,expected in manifest['files'].items():
        path=ROOT/name
        if not path.is_file():raise FileNotFoundError(path)
        actual=hashlib.sha256(path.read_bytes()).hexdigest()
        if actual!=expected:raise AssertionError(f'Hash mismatch: {name}')
    training=ROOT/'experiments/training'
    ids=set(np.load(training/'features/mace_polar_m_features.npz')['config_ids'].astype(str))
    geometry=list(frames(training/'geometries/configurations.extxyz'))
    assert len(geometry)==48 and {info['config_id'] for info,_ in geometry}==ids
    waters=[(len(x)-int(info['n_solute_atoms']))//3 for info,x in geometry]
    assert waters.count(3)==24 and waters.count(4)==24
    assert len({info['molecule_id'] for info,_ in geometry})==14
    checkpoint=hashlib.sha256((ROOT/'checkpoints/glider_site_response_ensemble.pt').read_bytes()).hexdigest()
    assert checkpoint=='a288afa285128e13cb7a05ba9459f1dce07e79b8638c6e272c2509ae68d7a288'
    for n,count in [(1,48),(2,96),(3,80)]:
        base=ROOT/f'experiments/panel_{n}'
        positions={info['config_id']:x for info,x in frames(base/'geometries/configurations.extxyz')}
        refs=pd.read_csv(base/'references/observable_registry.csv')
        assert len(positions)==len(refs)==count
        charge=moment=0.
        for row in refs.itertuples():
            q=np.load(base/'predictions/glider'/row.observable_file)
            charge=max(charge,abs(float(q['predicted_charges_e'].sum())))
            rebuilt=dipole_from_sites(positions[row.config_id],q['predicted_charges_e'],q['predicted_dipoles_e_bohr'])
            moment=max(moment,float(np.abs(rebuilt-q['predicted_dipole_debye']).max()))
        assert charge<1e-5 and moment<1e-6,(n,charge,moment)
        print(f'Panel {n}: {count} configurations, charge {charge:.2e} e, moment {moment:.2e} D')
    for cfg in (ROOT/'experiments/distance_sweep/configurations').glob('*.json'):
        def inspect(value):
            if isinstance(value,dict):
                for k,v in value.items():
                    if k in ['base_geometry','base_density','glider','mace','file'] and isinstance(v,str):
                        assert (ROOT/v).is_file(),(cfg,k,v)
                    inspect(v)
            elif isinstance(value,list):
                for item in value:inspect(item)
        inspect(json.loads(cfg.read_text()))
    contact=ROOT/'experiments/nonwater_contact'
    freeze=json.loads((contact/'prediction_freeze.json').read_text())
    reference_manifest=json.loads((contact/'references/REFERENCE_MANIFEST.json').read_text())
    digest=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    design=json.loads((contact/'design_manifest.json').read_text())
    assert freeze['new_qm_reference_labels_accessed'] is False
    assert digest(contact/'configurations.extxyz')==freeze['geometry_sha256']
    assert digest(contact/'design_manifest.json')==freeze['design_manifest_sha256']
    assert digest(ROOT/'scripts/diagnostics/generate_nonwater_contact_panel.py')==design['generator_sha256']
    assert reference_manifest['prediction_freeze_sha256']==digest(contact/'prediction_freeze.json')
    assert reference_manifest['n_configurations']==36
    with (contact/'configuration_registry.csv').open(newline='') as stream:
        contact_cases=list(csv.DictReader(stream))
    assert len(contact_cases)==36
    assert len({row['molecule_id'] for row in contact_cases})==12
    assert {row['perturbant'] for row in contact_cases}=={'NH3','CH3OH','CH3CN'}
    with (contact/'references/REFERENCE_REGISTRY.csv').open(newline='') as stream:
        contact_references=list(csv.DictReader(stream))
    assert len(contact_references)==36
    assert {row['config_id'] for row in contact_cases}=={row['config_id'] for row in contact_references}
    assert digest(contact/'references/REFERENCE_REGISTRY.csv')==reference_manifest['reference_registry_sha256']
    reference_records=[json.loads(path.read_text()) for path in (contact/'references').glob('*.json')
                       if path.name!='REFERENCE_MANIFEST.json']
    assert len(reference_records)==36
    assert {row['config_id'] for row in reference_records}=={row['config_id'] for row in contact_cases}
    assert all(row['all_components_converged'] for row in reference_records)
    for method,key in [('glider','glider_predictions'),('mace_polar_l','mace_polar_l_predictions')]:
        prediction_dir=contact/'predictions'/method
        registry=prediction_dir/'prediction_registry.csv'
        assert digest(registry)==freeze[key]['registry_sha256']
        with registry.open(newline='') as stream:
            predictions=list(csv.DictReader(stream))
        assert len(predictions)==36
        for row in predictions:
            assert digest(prediction_dir/row['prediction_file'])==freeze[key]['case_file_sha256'][row['config_id']]
    for row in contact_references:
        assert digest(contact/'references'/row['observable_file'])==row['observable_sha256']
    print('Contact follow-up: 36 complete frozen geometries, predictions and QM references')
    print(f"PASS: {len(manifest['files'])} content hashes; original 24 + 24 training split; canonical input paths")
if __name__=='__main__':main()
