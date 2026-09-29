# Scripts

| Path | Task |
|:--|:--|
| [`benchmark/`](benchmark/README.md) | Run the frozen model on a geometry or cached features. |
| [`reproduce/`](reproduce/README.md) | Fetch hash-checked external weights, verify archived predictions and recompute panel statistics. |
| [`render_figures.py`](render_figures.py) | Rebuild the repository's SVG reading guide from published CSV values. |
| [`verify_companion.py`](verify_companion.py) | Check companion completeness, original hashes, headline values and links. |

The model-related scripts preserve the original scientific implementation and use the same checkpoint. The gallery script draws an explanatory overview; it is not part of the benchmark scoring pipeline.
