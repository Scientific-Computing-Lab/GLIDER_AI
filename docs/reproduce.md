# Reproduce a result

Run commands from the repository root with Python 3.11 or newer. Outputs go into `build/` unless explicitly specified.

## 1. Verify the released inputs

```bash
python -m pip install -e .
python scripts/reproduce/verify_release.py
python scripts/verify_companion.py
```

Every required file in the current release manifest must exist and match its hash. The verifier also checks the exact 24/24 water-count split, frozen checkpoint and charge/moment constraints on all 224 primary predictions.

## 2. Re-score the raw prediction arrays

```bash
python scripts/reproduce/score_experiment.py --experiment all
python scripts/reproduce/recompute_statistics.py --output build/statistics
```

The first command reconstructs per-configuration and equal-solute response errors for the three panels, non-water, liquid and shell-size tests. The second computes the primary paired solute-bootstrap intervals from the released solute tables. Scoring results can be compared with `experiments/*/results.csv`.

The non-water panel's subsequent fixed-geometry sanity check uses the stored
counterpoise component energies and geometry registry:

```bash
uv run --no-project --with 'ase==3.29.0' python scripts/reproduce/audit_nonwater_contacts.py
python scripts/figures/plot_nonwater_contact_audit.py
```

It reproduces [the 72-case contact audit](../experiments/nonwater/) without
discarding or changing any archived prediction or QM reference. The audit is
post hoc and changes the scope of the transfer claim.

The separate, geometry-corrected follow-up can be rescored directly from
its released 36 QM references and both frozen prediction sets:

```bash
python scripts/diagnostics/score_nonwater_contact_panel.py --base experiments/nonwater_contact
```

This command checks the pre-QM prediction freeze and every case-file hash,
then recomputes the case, equal-solute and paired-bootstrap summaries. The
follow-up was designed after the original contact audit. It does not replace
the 72-case score or retroactively make that test a realistic-contact panel.
To rerun the geometry construction or expensive reference SCFs, install
`.[contact]` and follow the chronological recipe in
[the diagnostic scripts](../scripts/diagnostics/README.md).

## 3. Recompute a physical coupling from QM densities

```bash
python -m pip install -e '.[qm]'
python scripts/reproduce/recompute_coupling.py --solute methane --water-rank 4
```

This performs fresh CPU Coulomb integrals using the released density matrices and predicted response sites. It checks agreement with the archived physical values. It is a density-to-energy test, not a new SCF calculation.

## 4. Predict a new response

Follow the [geometry inference guide](../scripts/benchmark/README.md). It specifies the MACE commit, official checkpoints and complete invocation. The example output includes surface coordinates, potential, dipole and atom-centred response sites.

## 5. Recompute the distance sweep

Use the [distance-sweep command](../experiments/distance_sweep/README.md). Configurations point to released base and probe arrays. The CPU path is slower for large bases and for full quantum-water integration. All original accepted probe placements and evaluated arrays are already available for independent rescoring.

## Author-review controls

After the model and official M/L dependencies are installed:

```bash
python scripts/diagnostics/run_controls.py ablation --output build/ablation --device cuda
python scripts/diagnostics/run_controls.py dissociation --output build/separation --device cuda
```

The wrapper runs the entire sequence: fixed-schedule training for the ablation, geometry construction, M-feature extraction, M/L prior prediction, averaging and scoring. For the separation test, use `--distances 50 100` with a different output folder to reproduce the extension. Output folders must be new. The [underlying operations](../scripts/diagnostics/review_controls.py) can also be run separately. CPU is supported.

## Direct QM check of one separation trajectory

The [20-case `dev_cyclic_carbamate` package](../experiments/dissociation_qm/) contains the geometries, 512 fixed probes, frozen predictions, complete complex and ghost-fragment potential arrays, component energies and one result row per geometry. The standard `scripts/verify_companion.py` command rederives all response and error RMS values from these arrays and checks the source hashes. It needs no SCF run.

To repeat the expensive reference calculations on CPU, install the QM and figure dependencies and write into `build/`:

```bash
python -m pip install -e '.[contact,figures]'
python scripts/diagnostics/run_cyclic_carbamate_separation_qm.py \
  --output build/separation_qm_dev_cyclic_carbamate
python scripts/diagnostics/summarize_cyclic_carbamate_separation_qm.py \
  build/separation_qm_dev_cyclic_carbamate
```

The script uses the unchanged archived geometry and prediction arrays from `dissociation/` and `dissociation_extended/`. One full calculation has 60 component SCFs and may take substantial CPU time. The published arrays were generated with the same reference Hamiltonian on GPU; independent CPU reruns of 13 cases are compared in the [experiment guide](../experiments/dissociation_qm/).

## Rebuild the gallery

```bash
python -m pip install -e '.[figures]'
python scripts/render_figures.py
```

The gallery uses the current figure-numbered CSVs. Source manuscript files are provided separately to the authors; this repository distributes scientific code, models, data and documentation.

## Numerical reproducibility of fresh inference

Rescoring the stored arrays reproduces the published numbers. Re-extracting neural features on a fresh runtime need not be bitwise identical to the original run. In the released 28-atom example, fresh CUDA inference gave ESP NRMSE 0.319834 versus archived 0.319858, with maximum potential difference 3.20 × 10⁻⁶ Eh/e and maximum dipole-component difference 1.54 × 10⁻⁴ D. The matched ablation compares variants on one common newly extracted cache; its full-head check is exact on that cache. [Measured differences](../provenance/validation/fresh_inference_comparison.json).
