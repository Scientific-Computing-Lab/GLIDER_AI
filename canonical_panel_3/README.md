# Prospective Panel III · independent-source solute identities

This panel evaluates **20 solutes and 80 three-water configurations** selected from an identity-only extraction of pinned FreeSolv v0.52. FreeSolv supplied molecular identities, **not hydration free energies**. Geometry, GLIDER and comparator predictions were frozen before QM response references were produced. The model was not retrained. The [preregistration](PREREGISTRATION.md), [selection manifest](selection/SELECTION_MANIFEST.json), [prediction freeze](predictions/ALL_PREDICTIONS_FREEZE_MANIFEST.json) and [reference manifest](references/REFERENCE_MANIFEST.json) retain the chronology and exclusions.

| Navigate | Contents |
|:--|:--|
| [`source/`](source/README.md), [`selection/`](selection/README.md) | Identity-only source, eligibility, deterministic selection and exclusions. |
| [`data/`](data/README.md) | Frozen geometries and configuration registry. |
| [`predictions/`](predictions/README.md) | GLIDER, public comparators and pre-QM freeze manifest. |
| [`references/`](references/README.md) | Later QM response observables and recovery audit. |
| [`audits/`](audits/README.md), [`scripts/`](scripts/README.md) | Rotation/comparator checks and original panel procedures. |

In the camera-ready paper (Figure 3), GLIDER response-ESP NRMSE is **0.427** against **0.630** for the MACE-L polar baseline, with **19/20** solutes improved. The global response-dipole ranking is less conclusive: its closest paired interval includes zero. Read [`FINAL_REPORT.md`](FINAL_REPORT.md) with the [current paper map](../results/README.md); report language reflects the frozen panel's original interpretation and does not expand the workshop paper's claims. Neither this panel nor its source data constitute a hydration-free-energy benchmark.
