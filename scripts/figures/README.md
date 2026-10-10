# Rebuild the experiment graphics

Run these commands from the repository root. Each script reads released
experiment data and writes an editable SVG beside the corresponding
experiment; none reruns quantum chemistry or changes the predictions.

| Result | Command | Output guide |
|:--|:--|:--|
| Journal Fig. 9 reading aid, including the ARROW-covered comparison and shell-size bias | `python scripts/figures/plot_arrow_comparison.py` | [ARROW comparison](../../experiments/arrow_comparison/) |
| Original 72-case non-water geometry audit | `python scripts/figures/plot_nonwater_contact_audit.py` | [Original non-water test](../../experiments/nonwater/) |
| Separate 36-contact non-water follow-up | `python scripts/figures/plot_nonwater_contact_followup.py` | [Contact follow-up](../../experiments/nonwater_contact/) |

The printed paper figure is built from its native LaTeX source; these SVGs
make the same numerical stories easy to inspect in the data companion.
