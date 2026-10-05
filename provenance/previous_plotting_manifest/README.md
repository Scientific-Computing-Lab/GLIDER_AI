# Camera-ready figure tables

These small CSVs are the exact numeric inputs to plotted values. They are copied byte-for-byte from the verified workshop camera-ready source. The [paper-to-data map](../../results/README.md) resolves the historical filenames (`Fig4_*` for camera-ready Figure 3, for example).

| File family | Main use |
|:--|:--|
| `Fig2_controls.csv` | Sparse-supervision development controls in Figure 2b,c. |
| `Fig4_*.csv` | Three prospective panel aggregates, eligible comparators and paired Figure 3 comparisons. |
| `Fig5_*.csv` | Chemistry, regime and new-neighbour extensions. |
| `Fig6_*`, `Fig7_*` | Further controls and ancillary plotted comparisons retained in the appendix. |
| `FigS*.csv` | Supplementary figures, including representation error and distance tests. |
| [`downstream/`](downstream/README.md) | W4 and W4–W12 held-out-water plotting tables. |

The [SVG renderer](../../scripts/render_figures.py) reads selected tables here and makes the repository's visual overview. It does not recalculate the benchmark or change any source CSV.
