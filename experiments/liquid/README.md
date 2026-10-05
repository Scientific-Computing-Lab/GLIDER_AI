# Local clusters from explicit-water sampling

24 configurations: six solutes × four temporally spaced snapshots, retaining three waters. The six identities are outside the original 14-solute response-training set. The source simulations were part of an earlier energy-model development campaign. These records establish absence from response-head fitting, not complete historical independence from every research decision.

**Paper:** Fig. S16a–d. [Full paper map](../../docs/paper_map.md).

| Open | Contents |
|:--|:--|
| [Geometries](geometries/) | ExtXYZ coordinates, solute identities, water counts and configuration registry. |
| [QM references](references/) | Full per-configuration reference potentials, probe coordinates and dipoles. |
| [Predictions](predictions/) | Frozen outputs, grouped by method, on the same probe coordinates. |
| [Results](results.csv) | Recomputed per-configuration response errors. |
| [Solute results](solute_results.csv) | Configurations averaged within solute. |
| [Summary](summary.csv) | Equal-solute aggregate metrics. |

```bash
python scripts/reproduce/score_experiment.py --experiment liquid
```

Run from the repository root. The script reads the raw reference and prediction arrays and writes a fresh result under `build/recomputed/liquid/`. [Units and NPZ fields](../../docs/data_format.md).

[Selection and label chronology](../../docs/data_lineage.md) · [Back to experiments](../README.md)
