# Panel III · twenty FreeSolv identities

80 configurations: 20 solutes × four three-water arrangements. Selection used identities and connectivity from a pinned FreeSolv source without accessing hydration targets. All eligible predictions were frozen before QM. The difficult retained outlier is included.

**Paper:** Fig. 3a–c. [Full paper map](../../docs/paper_map.md).

| Open | Contents |
|:--|:--|
| [Geometries](geometries/) | ExtXYZ coordinates, solute identities, water counts and configuration registry. |
| [QM references](references/) | Full per-configuration reference potentials, probe coordinates and dipoles. |
| [Predictions](predictions/) | Frozen outputs, grouped by method, on the same probe coordinates. |
| [Results](results.csv) | Recomputed per-configuration response errors. |
| [Solute results](solute_results.csv) | Configurations averaged within solute. |
| [Summary](summary.csv) | Equal-solute aggregate metrics. |

```bash
python scripts/reproduce/score_experiment.py --experiment panel_3
```

Run from the repository root. The script reads the raw reference and prediction arrays and writes a fresh result under `build/recomputed/panel_3/`. [Units and NPZ fields](../../docs/data_format.md).

[Selection and label chronology](../../docs/data_lineage.md) · [Back to experiments](../README.md)
