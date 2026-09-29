# Visual guide

All assets here are SVG, so the gallery stays crisp on GitHub and in presentations. Run `python scripts/render_figures.py` with the `figures` extra to regenerate them from the [source tables](../data/figure_data/README.md).

| Asset | Kind | Data or role |
|:--|:--|:--|
| [`glider-banner.svg`](glider-banner.svg) | Conceptual | Decorative molecule and field; not a QM or model output. |
| [`architecture.svg`](architecture.svg) | Conceptual | An explanatory flow from frozen MACE inputs to constrained response. |
| [`prospective-transfer.svg`](prospective-transfer.svg) | Quantitative | Exact `Fig4_glider.csv`, `Fig4_mace_polar_l.csv` and `Fig4_paired.csv` values. |
| [`frozen-coupling.svg`](frozen-coupling.svg) | Quantitative | Paper W4 aggregate values and `downstream/energy_rank.csv`. |
| [`distance-and-moments.svg`](distance-and-moments.svg) | Quantitative | `FigS15_dense_energy.csv` fixed-charge sweep. |

The SVGs are a repository-specific reading aid. They do not replace the paper's full panels, uncertainty analysis or captions; follow the [paper map](../results/README.md) for those.
