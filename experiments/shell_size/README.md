# Nested water environments

36 configurations: three solutes × three parent snapshots × 1, 3, 6 and 12 waters. Each size uses the same solute and nested neighbour list for its parent. The original response training contains both three- and four-water environments. Six and twelve waters exceed that maximum.

**Paper:** Fig. S17a,b. [Full paper map](../../docs/paper_map.md).

| Open | Contents |
|:--|:--|
| [Geometries](geometries/) | ExtXYZ coordinates, solute identities, water counts and configuration registry. |
| [QM references](references/) | Full per-configuration reference potentials, probe coordinates and dipoles. |
| [Predictions](predictions/) | Frozen outputs, grouped by method, on the same probe coordinates. |
| [Results](results.csv) | Recomputed per-configuration response errors. |
| [Solute results](solute_results.csv) | Configurations averaged within solute. |
| [Summary](summary.csv) | Equal-solute aggregate metrics. |

```bash
python scripts/reproduce/score_experiment.py --experiment shell_size
```

Run from the repository root. The script reads the raw reference and prediction arrays and writes a fresh result under `build/recomputed/shell_size/`. [Units and NPZ fields](../../docs/data_format.md).

[Selection and label chronology](../../docs/data_lineage.md) · [Back to experiments](../README.md)
