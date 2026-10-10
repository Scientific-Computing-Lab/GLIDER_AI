# ARROW source tables

These five folders preserve the collaborator's case-level tables in the
paper's experiment order. Each corresponds to the geometry, QM reference and
prediction arrays in the sibling experiment directory. Use the
[comparison guide](../README.md) for the common-subset summary and the
[scorer](../../../scripts/reproduce/score_arrow_release.py) to regenerate the
reported response-ESP errors from the arrays.

| Folder | Paper experiment | Cases with ARROW atom types |
|:--|:--|--:|
| [`panel_1/`](panel_1/) | Later-reused Panel I | 4 |
| [`panel_2/`](panel_2/) | Prospective Panel II | 16 |
| [`panel_3/`](panel_3/) | Prospective Panel III | 16 |
| [`liquid/`](liquid/) | Liquid-derived clusters | 8 |
| [`shell_size/`](shell_size/) | Nested 1/3/6/12-water series | 24 |

Within each folder, `arrow_coverage.csv` lists available and unavailable
cases. `arrow_smeared_results.csv` uses the native electron-cloud potential;
`arrow_results.csv` evaluates the same response dipoles as points. The
`covered_subset_*` tables compare only configurations on which ARROW,
GLIDER, MACE-L and QM all have values. These are historical scoring records;
the released arrays and current scorer are the numerical authority.
