# Final prospective campaign report — frozen failure

## Decision

**The preregistered prospective SOTA gate failed. This result must not be called
SOTA.**

On the prediction-before-label panel of 12 chemically unseen solutes and 48
structured three-water environments, the frozen site-resolved response model had
the lowest nominal error among all matched executable methods for both co-primary
observables:

| Observable | Frozen candidate | Strongest comparator | Relative reduction | Paired molecule-bootstrap 95% CI for candidate − comparator |
|---|---:|---:|---:|---:|
| Response-ESP NRMSE | **32.180%** | MACE-POLAR-1-L: 53.569% | **39.93%** | **[-26.365, -16.433] percentage points** |
| Induced-dipole vector RMSE | **0.18936 D** | MACE-POLAR-1-M: 0.19651 D | **3.64%** | **[-0.06078, +0.02965] D** |

The model improved both metrics for 11 of 12 molecules and was 156.5× faster
than the complete counterpoise QM reference workload. Nevertheless, three
mandatory conditions failed: absolute ESP NRMSE exceeded 30%, dipole improvement
was below 15%, and the dipole confidence interval crossed zero. No refitting,
metric change, exclusion, or second prospective attempt was performed.

## Exact failed gates

| Gate | Required | Observed | Result |
|---|---:|---:|---|
| Absolute response-ESP NRMSE | ≤30% | 32.180% | **FAIL** |
| Dipole reduction vs strongest comparator | ≥15% | 3.64% | **FAIL** |
| Paired dipole 95% CI | upper bound <0 | [-0.06078, +0.02965] D | **FAIL** |

All other frozen conditions passed: lowest nominal ESP and dipole errors; ≥20%
ESP reduction; ESP bootstrap support; 11/12 molecules improved on both metrics;
no orientation or compressed-stratum regression; charge conservation; rotation
stability; and ≥100× speedup.

## Comparator table

| Method | Molecule-mean ESP NRMSE | Dipole vector RMSE (D) |
|---|---:|---:|
| **Frozen candidate** | **32.180%** | **0.18936** |
| MACE-POLAR-1-L | 53.569% | 0.21556 |
| MACE-POLAR-1-M | 62.082% | 0.19651 |
| MACE-POLAR M/L mean | 55.517% | 0.19937 |
| Zero response | 100.000% | 0.47441 |
| SALTER-style MACE kernel | 111.271% | 0.39689 |
| Ridge, MACE features | 115.707% | 0.59526 |
| Ridge, ViSNet features | 117.347% | 0.57198 |
| Linear, ViSNet features | 282.137% | 1.62793 |
| Linear, MACE features | 874.161% | 5.77365 |
| Static q/u + Thole control | 20122.976% | 413.92839 |

The frozen stationary-response model was not a fair geometry-only comparator: it
requires per-solute QM observables at inference and does not cover the full panel.
It was therefore preregistered as ineligible, not removed after scoring.

## Regime-resolved result

The comparator is MACE-POLAR-1-L for ESP and MACE-POLAR-1-M for dipole.

| Regime | Candidate ESP | Comparator ESP | Candidate dipole (D) | Comparator dipole (D) |
|---|---:|---:|---:|---:|
| compressed | 25.274% | 48.683% | 0.05619 | 0.10466 |
| cooperative | 27.975% | 41.262% | 0.12124 | 0.14607 |
| equilibrium | 31.444% | 48.671% | 0.23496 | 0.21080 |
| orientation | 44.027% | 75.660% | 0.26527 | 0.27882 |

The candidate passed the explicit compressed and orientation no-regression gates.
Its dipole error regressed in the equilibrium stratum, explaining why the aggregate
dipole margin was small despite 11 molecule-level joint wins. The sole joint
molecule-level loss was the randomized max-min-selected iodoethanol chemistry
(`prospective_diverse_02`), whose candidate ESP NRMSE was 74.419% and whose dipole
RMSE contribution exceeded MACE-POLAR-1-M.

## Frozen method

The candidate is a five-seed ensemble comprising a frozen MACE-POLAR-1-M
product-basis encoder, an unfitted equal-weight MACE-POLAR-1-M/L distributed q/u
response prior, a local anisotropic site-resolved induced-charge/induced-dipole
correction, exact total-charge projection, and a directly supervised extensive
equivariant molecular-dipole readout. Training used only counterpoise-consistent
response ESP and induced dipoles for 14 development chemistries and 48 environments.
No energy, force, hydration-free-energy, prospective-label, or encoder-fine-tuning
supervision was used.

Checkpoint SHA-256:
`a288afa285128e13cb7a05ba9459f1dce07e79b8638c6e272c2509ae68d7a288`.

Maximum prospective net-charge deviation was `9.40e-9 e`. The target-free rotation
audit passed with maximum errors of `6.20e-5 e` in charge, `2.40e-4 e bohr` in local
dipole, and `0.00580 D` in reconstructed molecular dipole.

## Prospective protocol and integrity

- Method/code freeze: SHA-256
  `a762e97273c98e8cba80f21f765e38b8e6e73b68e1f5f648cf29b2cf131ff366`.
- All 550 candidate/baseline prediction artifacts were frozen before reference QM:
  SHA-256
  `8407beb6760d060a2c96c46f135ef3933febfc27b9a873d875bf052c3892730f`.
- Prospective reference freeze: SHA-256
  `22e49dcd0366044109837bd1af7ee5a7241c90f26de92a0c68853c0627f9e63a`.
- Reference: counterpoise-consistent ωB97X-D3BJ/def2-TZVPD, grid level 4,
  SCF convergence `1e-10`.
- Twelve chemistry blocks, four fixed regimes per block; no exclusions.
- The five new chemistries were selected by the frozen max-min Morgan/Mordred
  procedure. A post-freeze Pandas Boolean-to-float casting error was corrected
  before selection and before any label existed; the algorithm, ranking, model,
  and thresholds were unchanged and the erratum is preserved.
- One iodine chemistry required an exact SCF solver continuation. Final densities
  were converged under the unchanged Hamiltonian, basis, grid, zero level shift,
  and tolerance; its full recovery cost was retained in the timing.
- `EXPERIMENTAL_HYDRATION_TARGETS_ACCESSED = FALSE`.
- `PROSPECTIVE_LABELS_ACCESSED_BEFORE_PREDICTION_FREEZE = FALSE`.
- `MODEL_OR_METRIC_CHANGED_AFTER_PROSPECTIVE_LABEL_ACCESS = FALSE`.
- `MOLECULES_OR_CONFIGURATIONS_EXCLUDED_AFTER_SCORING = FALSE`.

The acquisition utility's legacy `manifest.json` retained a stale
“development-only” scope string. It was not rewritten after label generation; the
authoritative prospective scope and hashes are recorded in
`audit/prospective_reference_freeze_manifest.json`.

## What remains scientifically valid

Direct observable supervision produced a robust prospective improvement in the
spatial electronic-response field: 39.9% lower ESP NRMSE than the strongest matched
foundation-model comparator, with a molecule-blocked interval wholly favoring the
candidate and improvement across all four physical regimes. It also produced the
lowest nominal dipole error and 11/12 joint molecule wins. This supports the narrower
finding that site-resolved direct response supervision transfers useful electronic
field information under chemical shift.

It does **not** establish the preregistered prospective capability SOTA because the
absolute field accuracy and dipole effect size/statistical evidence were insufficient.
Per the campaign hard stop, there was no post-test tuning and no publication package
claiming SOTA.

## Reproducibility artifacts

The content-addressed negative release is in `frozen_failure_release/`. It contains
the checkpoint, all candidate and comparator predictions, all 48 reference
observables, configuration and molecule metrics, the panel, and the full audit
chronology. The authoritative machine-readable decision is
`frozen_failure_release/results/gate_result.json`.
