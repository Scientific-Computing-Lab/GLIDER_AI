# Panel II · twenty-four new solutes

96 configurations: 24 solutes × four three-water arrangements. The published GLIDER checkpoint and comparator predictions were frozen before response QM. The `later_global_prior` prediction folder is a separate later model trained with Panel-I labels. It is retained for comparison, not substituted for GLIDER.

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
python scripts/reproduce/score_experiment.py --experiment panel_2
```

Run from the repository root. The script reads the raw reference and prediction arrays and writes a fresh result under `build/recomputed/panel_2/`. [Units and NPZ fields](../../docs/data_format.md).

[Selection and label chronology](../../docs/data_lineage.md) · [Back to experiments](../README.md)
