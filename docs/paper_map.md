# Paper → data

The figure-numbered gallery below uses the SIMBIOCHEM camera-ready numbering
after the author-review revision. The three prospective panels are [Panel I](../experiments/panel_1/),
[Panel II](../experiments/panel_2/) and [Panel III](../experiments/panel_3/).
The updated journal manuscript numbers the external-probe sweep as
**Supplementary Fig. S15, Table S19 and Section 14**. Its Fig. S15 draws
from [workshop Fig. 5](../figures/figure_05/) and
[workshop Fig. S13](../figures/figure_S13/), with complete inputs under
[distance_sweep](../experiments/distance_sweep/). The S15 row below refers
to the workshop's *different* transfer figure.
The journal manuscript's **Supplementary Figs. S16–S17** show the subsequent
[non-water geometry audit](../experiments/nonwater/contact_audit.svg) and
[separate contact follow-up](../experiments/nonwater_contact/contact_followup.svg).
They preserve the original 72-case response score and add 36 new contacts.
Its **Supplementary Fig. S18** checks one fragment-separation trajectory
against [new QM references](../experiments/dissociation_qm/).

## Expanded journal manuscript

| Journal display | Scientific question | Released experiment |
|:--|:--|:--|
| Fig. 1 | Why a response field contains more spatial information than a dipole | [Panel II example](../experiments/panel_2/) and [response definition](../docs/model.md) |
| Fig. 2 | What the constrained model learns from 48 configurations | [Training and development](../experiments/training/) |
| Fig. 3 | How the three prospective panels were selected and frozen | [Panel I](../experiments/panel_1/) · [II](../experiments/panel_2/) · [III](../experiments/panel_3/) |
| Fig. 4 | Prospective response-ESP and dipole accuracy | [Panels I–III](../experiments/) |
| Fig. 5a,c,d | Chemistry, regime and per-solute breadth of the water-panel result | [Panels I–III](../experiments/) |
| Fig. 5b | Non-water contact follow-up, one geometry per solute and neighbour species | [Frozen 36-case follow-up](../experiments/nonwater_contact/) |
| Fig. 6 | Liquid-derived local environments | [Liquid bridge](../experiments/liquid/) |
| Fig. 7 | Frozen coupling to a held-out water | [Held-out water](../experiments/heldout_water/) |
| Fig. 8 | Nested water size and inducing-fragment separation | [Shell size](../experiments/shell_size/) · [dissociation](../experiments/dissociation/) · [50–100 Å extension](../experiments/dissociation_extended/) |
| Supplementary Fig. S16 | Why the original 72 replacements were a geometry stress test | [Complete 72-case audit](../experiments/nonwater/) |
| Supplementary Fig. S17 | Contact energies and all case scores in the corrected follow-up | [Frozen 36-case follow-up](../experiments/nonwater_contact/) |
| Supplementary Fig. S18 | Does the one-solute QM response vanish as the inducing water separates? | [20 QM references and frozen predictions](../experiments/dissociation_qm/) |

The original 72-case non-water scores remain in journal Supplementary Tables
S8–S9 and the [original experiment folder](../experiments/nonwater/). The
follow-up was designed after that panel's geometry audit and does not replace
its predeclared result.

## SIMBIOCHEM camera-ready gallery

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
| S15 | Chemistry, water regimes and the 36-contact follow-up in panel b | [CSV data](../figures/figure_S15/) | [nonwater_contact](../experiments/nonwater_contact/) · [water panels](../experiments/) |
| S16 | Liquid-derived local environments | [CSV data](../figures/figure_S16/) | [liquid](../experiments/liquid/) |
| S17 | Size and representation controls | [CSV data](../figures/figure_S17/) | [shell_size](../experiments/shell_size/) |
| S18 | Contact-energy audit of the original 72 replacements | [Case-level audit](../experiments/nonwater/contact_audit.csv) | [nonwater](../experiments/nonwater/) |
| S19 | All energies and per-case response scores in the 36-contact follow-up | [Case scores](../experiments/nonwater_contact/results.csv) | [nonwater_contact](../experiments/nonwater_contact/) |

## SIMBIOCHEM camera-ready tables

| Table | Contents | Data and definitions |
|:--|:--|:--|
| 1 | Extension summaries, including the original non-water challenge and separate contact follow-up | [nonwater](../experiments/nonwater/) · [nonwater_contact](../experiments/nonwater_contact/) |
| S1 | Absolute response scale | [panel_1](../experiments/panel_1/) |
| S2 | Original 48 training configurations | [training](../experiments/training/) |
| S3 | Prospective identities | [panel_3](../experiments/panel_3/) |
| S4 | Comparator definitions and provenance | [panel_1](../experiments/panel_1/) |
| S5 | Complete prospective leaderboards | [panel_1](../experiments/panel_1/) |
| S6 | Regime-specific prospective and liquid errors | [liquid](../experiments/liquid/) |
| S7 | Post hoc water-contact sensitivity of the prospective panels | [water_contact_audit](../experiments/water_contact_audit/) |
| S8–S9 | Original 72 unrelaxed non-water replacements: aggregate and species comparisons | [nonwater](../experiments/nonwater/) |
| S10 | Liquid identities | [liquid](../experiments/liquid/) |
| S11 | Same-supervision AIMNet2 control | [training](../experiments/training/) |
| S12–S13 | Nested size results | [shell_size](../experiments/shell_size/) |
| S14 | Controlled global-branch ablation | [global_branch](../experiments/global_branch/) |
| S15 | Held-out-water identities and structures | [heldout_water](../experiments/heldout_water/) |
| S16 | W4 coupling | [heldout_water](../experiments/heldout_water/) |
| S17 | Exact multipole control | [heldout_water](../experiments/heldout_water/) |
| S18 | Frozen ensemble-member robustness | [panel_1](../experiments/panel_1/) |
| S19 | Reported artifact hashes | [training](../experiments/training/) |
| S20 | Distance–orientation sweep | [distance_sweep](../experiments/distance_sweep/) |

Unique historical derivation tables, including detailed development-control and ensemble-member analyses, are retained in [provenance](../provenance/previous_release_tables/). They are supporting analysis records rather than a second set of authoritative experiment data. The [filename map](../provenance/figure_filename_map.json) explains old plotting identifiers.
