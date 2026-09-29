# Visual guide

The first visual in the [main README](../README.md) comes from the paper's **actual QM Figure 1c**. Its four [transparent source layers](figure1_source/README.md) are overlaid without changing the scientific map; `python scripts/render_readme_field.py` adds labels and framing to make [`spatial-response.png`](spatial-response.png).

The remaining gallery is vector artwork. Run `python scripts/render_figures.py` with the `figures` extra to recreate it from the [numeric source tables](../data/figure_data/README.md).

| Asset | Kind | What it shows |
|:--|:--|:--|
| [`architecture.svg`](architecture.svg) | Conceptual | Frozen polar inputs, learned response corrections and constrained spatial output. |
| [`prospective-transfer.svg`](prospective-transfer.svg) | Quantitative | Response-ESP NRMSE in three prospective panels, from `Fig4_glider.csv`, `Fig4_mace_polar_l.csv` and `Fig4_paired.csv`. |
| [`frozen-coupling.svg`](frozen-coupling.svg) | Quantitative | W4 energy MAE and the W4–W12 rank series from the released downstream CSV. |
| [`distance-and-moments.svg`](distance-and-moments.svg) | Quantitative | Fixed-charge far-field control from `FigS15_dense_energy.csv`. |

The gallery is a reading aid, not a replacement for the paper's complete panels, uncertainty analysis or captions. Use the [paper-to-data map](../results/README.md) for an exact panel-to-source route.
