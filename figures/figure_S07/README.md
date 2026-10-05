# Figure S7 · Predicted response as the inducing fragments separate

The [plot values](separation.csv) combine the archived 3–20 Å and 50–100 Å summaries. Each value is the mean over 14 solutes of the predicted response-potential RMS on a fixed solute-centred probe surface. The curves compare GLIDER and its frozen averaged prior; they are **not** finite-distance errors against new QM labels. The one-water protocol contains only the solute and one water. The intact-environment protocol contains the solute and its original three- or four-water cluster.

Full geometries, probe points, predictions and per-case results are in the [3–20 Å experiment](../../experiments/dissociation/) and [50–100 Å extension](../../experiments/dissociation_extended/). The extension also archives the [per-case reconciliation diagnostic](../../experiments/dissociation_extended/posthoc_components.csv) cited in Appendix K.4. Figure S7 in the revised paper displays all ten distances. [Paper map](../../docs/paper_map.md).

To rebuild this table from the two archived summaries, run `python scripts/reproduce/build_figure_s07.py` at the repository root. To regenerate the SVG reading aid from the table, run `python -c 'from scripts.render_figures import separation; separation()'`.
