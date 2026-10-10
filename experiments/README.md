# Choose an experiment

Start from the scientific question, then open the corresponding guide. Raw array filenames are stable identifiers resolved by their registries. The [data dictionary](../docs/data_format.md) explains units and schemas.

| Experiment | Configurations | Purpose |
|:--|--:|:--|
| [Training](training/) | 48 | The only response labels used to fit the published checkpoint. |
| [Panel I](panel_1/) | 48 | First prospective solute test and retained joint-gate failure. |
| [Panel II](panel_2/) | 96 | A separately selected 24-solute test. |
| [Panel III](panel_3/) | 80 | Twenty FreeSolv identities, without hydration values. |
| [Water-contact audit](water_contact_audit/) | 224 screened, 9 flagged | Eight short solute–water contacts and one separate water–water contact; original panel scores remain primary. |
| [Non-water geometry challenge](nonwater/) | 72 | Frozen NH₃, CH₃OH and CH₃CN replacements, with a post hoc QM contact audit. |
| [Non-water contact follow-up](nonwater_contact/) | 36 | Separately frozen, fixed-solute MMFF94s contact geometries with new QM references. |
| [Liquid](liquid/) | 24 | Four explicit-water configurations per solute. |
| [Shell size](shell_size/) | 36 | Nine parents evaluated with 1, 3, 6 and 12 waters. |
| [ARROW comparison](arrow_comparison/) | 68 shared | Post hoc, atom-type-limited independent polarizable-model comparison on five source sets. |
| [Held-out water](heldout_water/) | 10 bases, 90 outer waters | Frozen response coupling and directional observables. |
| [Pilot](heldout_water_pilot/) | 6 bases | Original downstream protocol-development test. |
| [Distance sweep](distance_sweep/) | 10 bases, many placements | Spatial range of frozen-field reuse. |
| [Global branch](global_branch/) | 224 evaluation cases | A matched post hoc training intervention. |
| [Dissociation](dissociation/) | 224 | Move the inducing fragment itself to 3–20 Å. |
| [Extended dissociation](dissociation_extended/) | 56 | Follow-up at 50 and 100 Å. |
| [Dissociation QM check](dissociation_qm/) | 20 | Direct complex-minus-fragments QM references for one solute across the same 3–100 Å trajectories. |

## Prospective solute transfer

Panels I–III contain 56 different solutes and 224 three-water geometries. Their original GLIDER predictions preceded their response QM. They all evaluate the same frozen checkpoint. Panel I was subsequently reused for a different model family, which is clearly separated in [data lineage](../docs/data_lineage.md). The original checkpoint was not refit.

All reported results remain, including the unsuccessful joint gate and the later model that failed to improve.
