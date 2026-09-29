# Figure 1c · the QM field layers

These four transparent PNGs are **unchanged layers from Figure 1c of the camera-ready paper**. For each water arrangement, `field_*.png` contains the QM response-ESP point cloud and `molecule_*.png` contains the matching atom/bond overlay. The pairs share pixel dimensions and alignment.

| Arrangement | Response-potential layer | Molecular overlay | Response-dipole magnitude |
|:--|:--|:--|--:|
| Equilibrium | [`field_equilibrium.png`](field_equilibrium.png) | [`molecule_equilibrium.png`](molecule_equilibrium.png) | 0.464 D |
| Cooperative | [`field_cooperative.png`](field_cooperative.png) | [`molecule_cooperative.png`](molecule_cooperative.png) | 0.468 D |

[`scripts/render_readme_field.py`](../../scripts/render_readme_field.py) overlays each matched pair and adds only labels, frames and a qualitative sign key to make [`spatial-response.png`](../spatial-response.png). It does not recalculate or alter the QM fields. Blue represents negative response potential and red positive response potential. Similar *magnitudes* do not imply equal full dipole vectors. Full colour-scale units and scientific caption are in the paper.
