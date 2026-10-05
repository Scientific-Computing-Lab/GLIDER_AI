# Panel I · twelve new solutes

48 configurations: 12 solutes × four three-water arrangements. Seven identities came from an earlier sealed design list, and five were selected by chemistry-space diversity after the GLIDER freeze. All GLIDER predictions preceded response QM. Some public comparators were added later without response fitting. The broader preregistered joint gate failed.

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
python scripts/reproduce/score_experiment.py --experiment panel_1
```

Run from the repository root. The script reads the raw reference and prediction arrays and writes a fresh result under `build/recomputed/panel_1/`. [Units and NPZ fields](../../docs/data_format.md).

[Selection and label chronology](../../docs/data_lineage.md) · [Back to experiments](../README.md)
