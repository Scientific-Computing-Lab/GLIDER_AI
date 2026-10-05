# Changing the molecular neighbours

72 configurations: 12 Panel-I solutes × three neutral species (NH₃, CH₃OH, CH₃CN) × two arrangements. The water-trained checkpoint is unchanged. Neighbour number and identity both change. This is a combined environment-transfer test.

**Paper:** Fig. S15b. [Full paper map](../../docs/paper_map.md).

| Open | Contents |
|:--|:--|
| [Geometries](geometries/) | ExtXYZ coordinates, solute identities, water counts and configuration registry. |
| [QM references](references/) | Full per-configuration reference potentials, probe coordinates and dipoles. |
| [Predictions](predictions/) | Frozen outputs, grouped by method, on the same probe coordinates. |
| [Results](results.csv) | Recomputed per-configuration response errors. |
| [Solute results](solute_results.csv) | Configurations averaged within solute. |
| [Summary](summary.csv) | Equal-solute aggregate metrics. |

```bash
python scripts/reproduce/score_experiment.py --experiment nonwater
```

Run from the repository root. The script reads the raw reference and prediction arrays and writes a fresh result under `build/recomputed/nonwater/`. [Units and NPZ fields](../../docs/data_format.md).

[Selection and label chronology](../../docs/data_lineage.md) · [Back to experiments](../README.md)
