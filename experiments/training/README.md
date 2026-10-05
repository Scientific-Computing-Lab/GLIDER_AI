# The 48 response-training configurations

This folder contains exactly the examples used by the published response head. It does not contain the later Panel-I development pool.

| Solute group | Solutes | Arrangements per solute | Waters | Configurations |
|:--|--:|:--|--:|--:|
| Stress chemistries | 6 | compact, cooperative, orientation, compressed | 3 | 24 |
| Development chemistries | 8 | compact, cooperative, orientation | 4 | 24 |
| **Total** | **14** | | **3 or 4** | **48** |

The stored `equilibrium` regime corresponds to the compact motif. Solute conformers were built with RDKit ETKDG and MMFF relaxation, then waters were placed with deterministic geometry rules. These are constructed clusters rather than MD snapshots.

[Geometries and identities](geometries/) · [QM labels](references/) · [Frozen MACE features and averaged prior](features/) · [Original development scores](original_validation/)

Each configuration supplies a potential on 589–916 retained surface points (35,840 points altogether) and a three-component response dipole. The training loss uses 384 deterministically selected potential points per example. Constrained fitting to the QM field and dipole yields auxiliary site-charge/site-dipole targets. These fitted sites are derived labels, not independent QM measurements. The response head receives no energy or force labels.

The final five-member ensemble uses 64 epochs per seed. [Model and training details](../../docs/model.md). Leave-one-chemistry-out evaluation separated held chemistries during development. The final checkpoint uses all 48 examples.

Earlier archives combined these examples with Panel I under development terminology. The exact training membership is recovered from the frozen feature IDs and checked against the 48 released geometries. [Why the later pool exists](../../docs/data_lineage.md).
