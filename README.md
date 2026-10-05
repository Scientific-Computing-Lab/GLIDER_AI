<div align="center">

# GLIDER
### Learning the electronic response between molecules

**Sparse supervision · spatial fields · local physical reuse**

Gal Oren · Boris Fain · Michael Levitt<br>
Accepted at the **2nd SIMBIOCHEM Workshop at NeurIPS 2026**

[Understand the model](#from-molecular-geometry-to-a-response-field) · [Explore an experiment](#find-the-result-you-want) · [Reproduce a result](#start-with-a-reproducible-result) · [Model limits](#what-the-current-model-can-and-cannot-do)

<img src="assets/spatial-response.png" alt="Two QM response-potential maps: similar dipole magnitudes, different spatial patterns" width="100%">

*Two water arrangements. Almost the same response-dipole magnitude. Different places for a neighbour to interact.*

</div>

A molecule changes the electronic environment around it when another molecule approaches. A molecular dipole summarizes part of that change. A **response field** tells us where it occurs.

GLIDER learns the electrostatic potential associated with this rearrangement. It starts from pretrained polar representations and learns from **48 labelled configurations of 14 solutes**: 24 have three water neighbours and 24 have four. The response head is then frozen before the primary transfer tests.

## From molecular geometry to a response field

<img src="assets/architecture.svg" alt="Frozen MACE features and an averaged M/L response prior feed learned local and global corrections" width="100%">

The target is a fixed-geometry difference: evaluate the complex and its separated fragments in the same quantum-chemical basis, then subtract.

$$\Delta V_{\mathrm{resp}}(\mathbf r)=V_{SE}(\mathbf r)-V_S^{E\text{-ghost}}(\mathbf r)-V_E^{S\text{-ghost}}(\mathbf r).$$

This potential describes the **mutual electronic response** of solute and environment. GLIDER represents it with atom-centred charges and dipoles. The site charges sum to zero, and their combined moment agrees with a separately predicted response dipole.

MACE-POLAR-1-M supplies the frozen geometry features. An equal-weight average of MACE-POLAR-1-M and -1-L complex-minus-fragments predictions supplies the starting response. In comparisons, independently evaluated MACE-POLAR-1-L is the **polar baseline**. It is one contributor to the prior, rather than an independent physics model. [How training and inference work →](docs/model.md)

## What the experiments show

| Question | Evidence | Read the data |
|:--|:--|:--|
| Does it transfer to new solutes? | **32–40% lower response-ESP error** across 56 solutes and 224 configurations. | [Panels I–III](experiments/README.md#prospective-solute-transfer) |
| Does it extend beyond water? | **27.0% lower response-ESP error** on 12 solutes with three neutral neighbour species. | [Non-water environments](experiments/nonwater/) |
| Can another molecule use the field? | Held-out-water response-coupling MAE: **0.114 → 0.049 kcal mol⁻¹**. | [Held-out water](experiments/heldout_water/) |
| What does the global branch add? | A matched ablation finds a **modest average benefit**, with differences across solutes. | [Global-branch control](experiments/global_branch/) |
| Does the response vanish at separation? | **No.** A residual remains at 20 Å, including in the frozen prior. | [Dissociation diagnostic](experiments/dissociation/) |

<img src="assets/prospective-transfer.svg" alt="Three prospective panels show lower spatial-response error with GLIDER" width="100%">

The strongest evidence is for spatial response. Global-dipole gains are less consistent. Panel I missed its broader preregistered joint ESP-and-dipole gate, and all its cases remain available. [Selection, label reuse and retained failures →](docs/data_lineage.md)

## Put the field to use

Compute the response of a solute and its first three waters, then bring in a fourth water. That fourth water enters neither the response-model input nor the base-reference calculation. Its interaction with the frozen field tests whether the spatial prediction is useful beyond its fitting observable.

<img src="assets/frozen-coupling.svg" alt="Frozen coupling to held-out water and reuse from W4 through W12" width="100%">

An exact QM response dipole has an exact global moment, but its single-origin field can be inaccurate nearby. Higher multipoles recover with distance and eventually outperform GLIDER. This experiment measures a **frozen response contribution to Coulomb coupling**. It does not rank complete force fields by total interaction energy. [Explore distance and orientation →](experiments/distance_sweep/)

## Find the result you want

The release is organized around experiments. Each guide points directly to its geometries, references, predictions and results. The [paper map](docs/paper_map.md) connects the current figure and table numbers to these folders.

```text
experiments/
  training/             Original 48 configurations, labels and frozen features
  panel_1/              12 solutes · 48 configurations
  panel_2/              24 solutes · 96 configurations
  panel_3/              20 solutes · 80 configurations
  nonwater/             12 solutes · 72 environments
  liquid/               6 solutes · 24 liquid-derived clusters
  shell_size/           3 solutes · nested 1/3/6/12-water clusters
  heldout_water/        10 solutes · densities, probe states and coupling
  heldout_water_pilot/  Separate original six-solute pilot
  distance_sweep/       Probe positions, fields, densities and raw evaluations
  global_branch/       Matched training ablation and predictions
  dissociation/        Fragment-separation test at 3–20 Å
```

[Experiment catalogue](experiments/README.md) · [Array formats and units](docs/data_format.md) · [Current figure data](figures/) · [Frozen checkpoint](checkpoints/) · [Historical provenance](provenance/)

## Start with a reproducible result

```bash
git clone https://github.com/Scientific-Computing-Lab/GLIDER_AI.git
cd GLIDER_AI
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
python scripts/reproduce/verify_release.py
python scripts/reproduce/score_experiment.py --experiment all
python scripts/reproduce/recompute_statistics.py --output build/statistics
python scripts/verify_companion.py
```

This path needs no GPU or new quantum-chemistry calculation. It verifies file hashes and constraints and computes errors from the released reference and prediction arrays. Results go into `build/`, leaving the archive unchanged.

For a physical calculation from density matrices, install `.[qm]` and run:

```bash
python scripts/reproduce/recompute_coupling.py --solute methane --water-rank 4
```

The CPU recomputation agrees with the archived QM, GLIDER and baseline coupling values to better than **5 × 10⁻¹¹ kcal mol⁻¹**. [Reproduction levels and commands →](docs/reproduce.md)

For prediction on a new geometry, follow the [inference guide](scripts/benchmark/README.md). It installs the paper-matched MACE implementation and downloads the official M/L checkpoints from their publisher.

## What the current model can and cannot do

The released checkpoint supports analysis and frozen local reuse for the compact, neutral systems tested here. It conserves net response charge and matches its predicted global moment. **It does not enforce fragment separability:** moving the inducing water far away leaves a residual response. The [new diagnostic](experiments/dissociation/) records that failure, and an [extension to 50–100 Å](experiments/dissociation_extended/) checks that increasing separation does not resolve it.

A general molecular-simulation component needs a controlled dissociation limit and self-consistent coupling. Total-energy accuracy, ions and arbitrary condensed phases have not been established by these experiments.

## Cite and contact

**Gal Oren, Boris Fain and Michael Levitt.** *Sparse Supervision Turns Polar Pretraining into Transferable Molecular Response Fields.* 2nd SIMBIOCHEM Workshop at NeurIPS 2026. [Machine-readable citation](CITATION.cff).

Correspondence: [galoren@stanford.edu](mailto:galoren@stanford.edu) · [levittm@stanford.edu](mailto:levittm@stanford.edu)

Original release code is [MIT licensed](LICENSE). Third-party software, weights and source structures retain their own licences. [Sources and attribution](docs/sources.md). Generative AI assisted writing, code development, data processing and visual presentation. Scientific conclusions are tied to the calculations and archived records linked above.
