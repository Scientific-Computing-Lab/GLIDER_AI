# Water-panel short-contact audit

This **post hoc sensitivity check** retains the three frozen water panels as
the primary analysis. It asks whether their response-field comparison depends
on unusually short solute–water contacts. It does not alter a geometry,
reference, prediction, or model checkpoint.

Across the 224 archived configurations, eight have a nearest solute–water
atomic distance below 1.5 Å. A separate all-pairs check finds one additional
Panel-I case with a 1.279 Å H···H contact between two waters. All nine are in
the orientation regime. The [archived orientation generator](../../provenance/sampling/manybody_completion/src/build_clusters.py)
rotates each water about its oxygen without rechecking atom–atom clearance
after the rotation.
The threshold flags geometries for inspection; distance alone does not classify
an interaction as energetically repulsive. The [solute–water cases](short_contacts.csv)
and [water–water case](short_water_water_contacts.csv) list exact IDs, atom types,
distances, and response-ESP errors.

An energy-only counterpoise rerun at the response-reference level finds seven
of the eight flagged solute–water clusters net repulsive (interaction energies
+7.5 to +95.6 kcal/mol) and one net attractive (−11.1 kcal/mol). The short
water–water pair is repulsive when evaluated as an isolated frozen dimer
(+20.3 kcal/mol). The [eight cluster energies](cp_interaction_energies_short.csv)
and [water-pair calculation](short_water_pair_cp.json) include component
energies and calculation provenance. The isolated water-pair value does not
include its solute or third water. These nine targeted calculations do not
describe the interaction-energy distribution of all 224 configurations.

| Panel | Primary cases | Excluding short solute–water cases | Reduction, primary | Reduction, screened |
|---|---:|---:|---:|---:|
| I | 48 | 45 | 38.6% | 40.1% |
| II | 96 | 94 | 40.1% | 40.1% |
| III | 80 | 77 | 32.3% | 32.3% |

Excluding the ninth, water–water case as well leaves 44 Panel-I cases and a
41.3% reduction. Each original panel remains the primary result.

The [full result table](panel_sensitivity.csv) records both models' equal-solute
NRMSE, paired differences, 100,000-draw solute-blocked bootstrap intervals,
and solute wins. Every screened interval remains below zero. An affected solute
contributes its mean over the remaining three configurations, so each solute
still has equal weight. The original 48/96/80-case scores remain the primary
reported endpoints.

Reproduce with the repository environment:

```sh
python scripts/diagnostics/audit_water_contacts.py
python scripts/diagnostics/compute_water_cp_energies.py
python scripts/diagnostics/compute_short_water_pair_cp.py
```

Inputs are the `geometries/configuration_registry.csv`,
`geometries/configurations.extxyz`, and per-configuration results files under
[`panel_1`](../panel_1), [`panel_2`](../panel_2), and [`panel_3`](../panel_3).
The script independently recomputes each closest interfragment distance from
the archived coordinates and verifies it against the registry before scoring.

The released water-panel response-reference NPZ files do not contain the three
QM component **total energies** needed to derive counterpoise interaction
energies. The separate energy-only script
[`compute_water_cp_energies.py`](../../scripts/diagnostics/compute_water_cp_energies.py)
reruns the reference-level Hamiltonian on the archived geometries and records
converged component energies. It defaults to the eight flagged solute–water
clusters; `--all` can calculate the full 224-case distribution. The pairwise
energy uses [`compute_short_water_pair_cp.py`](../../scripts/diagnostics/compute_short_water_pair_cp.py).
The full distribution has not yet been computed, and no original response
reference or prediction was changed.
