# Data and archived observables

Start with the [paper-to-data map](../results/README.md). This directory contains the **published scientific inputs and results**, not a manuscript package.

| Subdirectory | What it contains | Use |
|:--|:--|:--|
| [`development/`](development/README.md) | Original 14-chemistry development geometries, QM observables, frozen feature/prior caches. | Inspect the sparse supervision and repeat the head fit with the frozen inputs. |
| [`prospective_1/`](prospective_1/README.md) | Panel I: 12 solutes, 48 configurations, predicted and QM observables. | Audit the first prospective panel. |
| [`prospective_2/`](prospective_2/README.md) | Panel II: 24 solutes, 96 configurations, predicted and QM observables. | Audit the second prospective panel. |
| [`figure_data/`](figure_data/README.md) | Compact CSVs read by the camera-ready plots, plus downstream plotting tables. | Replot or verify a figure. |
| [`source_data/`](source_data/README.md) | Per-figure underlying numeric data and spreadsheet. | Inspect points and paired comparisons. |

Panel III uses a separate label-blind identity source and lives at [`canonical_panel_3/`](../canonical_panel_3/README.md). The 56 prospective solutes are absent from GLIDER *response supervision*; foundation-model exposure is not established. Raw per-configuration predictions are compressed NumPy archives. The relevant registry maps each configuration ID to its file. Do not aggregate configurations as independent samples when reporting molecule-level uncertainty.
