# Paper → data

This map uses the SIMBIOCHEM camera-ready numbering after the author-review
revision. The three prospective panels are [Panel I](../experiments/panel_1/),
[Panel II](../experiments/panel_2/) and [Panel III](../experiments/panel_3/).
The updated journal manuscript numbers the external-probe sweep as
**Supplementary Fig. S15, Table S19 and Section 14**. Its Fig. S15 draws
from [workshop Fig. 5](../figures/figure_05/) and
[workshop Fig. S13](../figures/figure_S13/), with complete inputs under
[distance_sweep](../experiments/distance_sweep/). The S15 row below refers
to the workshop's *different* transfer figure.

| Figure | What it answers | Plot values | Raw experiment |
|:--|:--|:--|:--|
| 1 | Physical response, target and representation | [CSV data](../figures/figure_01/) | [training](../experiments/training/) |
| 2 | Architecture and development controls | [CSV data](../figures/figure_02/) | [training](../experiments/training/) |
| 3 | Prospective response-ESP and dipole comparisons | [CSV data](../figures/figure_03/) | [panel_1](../experiments/panel_1/) |
| 4 | Held-out-water energy, torque, orientation and rank | [CSV data](../figures/figure_04/) | [heldout_water](../experiments/heldout_water/) |
| 5 | Fixed-charge probe distance and exact moments | [CSV data](../figures/figure_05/) | [distance_sweep](../experiments/distance_sweep/) |
| S1 | Reference construction and probe surface | [CSV data](../figures/figure_S01/) | [training](../experiments/training/) |
| S2 | Representation reconstruction error floor | [CSV data](../figures/figure_S02/) | [training](../experiments/training/) |
| S3 | All 56 prospective solutes | [CSV data](../figures/figure_S03/) | [panel_1](../experiments/panel_1/) |
| S4 | Bootstrap and leave-one-solute sensitivity | [CSV data](../figures/figure_S04/) | [panel_1](../experiments/panel_1/) |
| S5 | Predicted and QM spatial fields | See experiment / provenance | [panel_3](../experiments/panel_3/) |
| S6 | Comparator observable coverage | [CSV data](../figures/figure_S06/) | [panel_1](../experiments/panel_1/) |
| S7 | Response as the inducing fragments separate, 3–100 Å | [CSV data](../figures/figure_S07/) | [3–20 Å](../experiments/dissociation/) and [50–100 Å](../experiments/dissociation_extended/) |
| S8 | Complementary W4 error and distance diagnostics | [CSV data](../figures/figure_S08/) | [heldout_water](../experiments/heldout_water/) |
| S9 | Exact multipole hierarchy | [CSV data](../figures/figure_S09/) | [heldout_water](../experiments/heldout_water/) |
| S10 | Torque across water rank | [CSV data](../figures/figure_S10/) | [heldout_water](../experiments/heldout_water/) |
| S11 | Cumulative frozen coupling | [CSV data](../figures/figure_S11/) | [heldout_water](../experiments/heldout_water/) |
| S12 | Response fraction of frozen electrostatics | [CSV data](../figures/figure_S12/) | [heldout_water](../experiments/heldout_water/) |
| S13 | Quantum-water sweep and paired uncertainty | [CSV data](../figures/figure_S13/) | [distance_sweep](../experiments/distance_sweep/) |
| S14 | Prediction chronology | See experiment / provenance | [panel_1](../experiments/panel_1/) |
| S15 | Chemistry, regime and non-water transfer | [CSV data](../figures/figure_S15/) | [nonwater](../experiments/nonwater/) |
| S16 | Liquid-derived local environments | [CSV data](../figures/figure_S16/) | [liquid](../experiments/liquid/) |
| S17 | Size and representation controls | [CSV data](../figures/figure_S17/) | [shell_size](../experiments/shell_size/) |

## Tables

| Table | Contents | Data and definitions |
|:--|:--|:--|
| 1 | Extension summaries | [nonwater](../experiments/nonwater/) |
| S1 | Absolute response scale | [panel_1](../experiments/panel_1/) |
| S2 | Original 48 training configurations | [training](../experiments/training/) |
| S3 | Prospective identities | [panel_3](../experiments/panel_3/) |
| S4 | Comparator definitions and provenance | [panel_1](../experiments/panel_1/) |
| S5 | Complete prospective leaderboards | [panel_1](../experiments/panel_1/) |
| S6 | Regime-specific prospective and liquid errors | [liquid](../experiments/liquid/) |
| S7–S8 | Non-water aggregate and species comparisons | [nonwater](../experiments/nonwater/) |
| S9 | Liquid identities | [liquid](../experiments/liquid/) |
| S10 | Same-supervision AIMNet2 control | [training](../experiments/training/) |
| S11–S12 | Nested size results | [shell_size](../experiments/shell_size/) |
| S13 | Controlled global-branch ablation | [global_branch](../experiments/global_branch/) |
| S14 | Held-out-water identities and structures | [heldout_water](../experiments/heldout_water/) |
| S15 | W4 coupling | [heldout_water](../experiments/heldout_water/) |
| S16 | Exact multipole control | [heldout_water](../experiments/heldout_water/) |
| S17 | Frozen ensemble-member robustness | [panel_1](../experiments/panel_1/) |
| S18 | Reported artifact hashes | [training](../experiments/training/) |
| S19 | Distance–orientation sweep | [distance_sweep](../experiments/distance_sweep/) |

Unique historical derivation tables, including detailed development-control and ensemble-member analyses, are retained in [provenance](../provenance/previous_release_tables/). They are supporting analysis records rather than a second set of authoritative experiment data. The [filename map](../provenance/figure_filename_map.json) explains old plotting identifiers.
