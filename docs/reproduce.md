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

## Rebuild the gallery

```bash
python -m pip install -e '.[figures]'
python scripts/render_figures.py
```

The gallery uses the current figure-numbered CSVs. Source manuscript files are provided separately to the authors; this repository distributes scientific code, models, data and documentation.

## Numerical reproducibility of fresh inference

Rescoring the stored arrays reproduces the published numbers. Re-extracting neural features on a fresh runtime need not be bitwise identical to the original run. In the released 28-atom example, fresh CUDA inference gave ESP NRMSE 0.319834 versus archived 0.319858, with maximum potential difference 3.20 × 10⁻⁶ Eh/e and maximum dipole-component difference 1.54 × 10⁻⁴ D. The matched ablation compares variants on one common newly extracted cache; its full-head check is exact on that cache. [Measured differences](../provenance/validation/fresh_inference_comparison.json).
