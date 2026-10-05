# Canonical Panel-III chemical-space contribution

Panel III was selected before prediction and scoring, so this analysis is
descriptive. The reference population is the 543 identities that passed the
frozen FreeSolv eligibility and historical-overlap filters.

## Quantitative placement

| Quantity | Panel III | Eligible FreeSolv population |
|---|---:|---:|
| identities | 20 | 543 |
| molecular weight (Da) | 84.16–260.39 | 54.09–295.34 |
| heavy atoms | 4–15 | 4–15 |
| TPSA (A²) | 0–110.38 | 0–121.38 |
| calculated logP | −3.22–3.75 | −3.59–5.61 |
| rings | 0–3 | 0–3 |
| heteroatoms | 0–8 | 0–8 |
| maximum Morgan similarity to historical GLIDER chemistry | 0.158–0.273 | <0.70 by eligibility |

The median selected maximum historical similarity is 0.182. Four molecules
come from each preregistered stratum. The panel adds ordinary nonpolar aliphatic
and aromatic solutes, oxygen-centred alcohol/ether chemistry, nitrogen-rich and
heteroaromatic structures, four sulfur/phosphorus structures, and four heavily
halogenated structures including Br and I. It includes donor-rich polyols,
multiple carbonyl/nitro motifs and mixed polar/nonpolar surfaces. No
carboxylic-acid identity was selected by the frozen algorithm; this remains a
specific functional-group gap rather than a post-hoc reason to alter the panel.

`source_data/Fig7_chemical_space.csv` contains every descriptor point and PCA
coordinate used in the descriptive figure. `source_data/Fig7_descriptor_percentiles.csv`
contains the eligible-population percentile data. Selection itself used the
preregistered Morgan/descriptor max–min distance, not PCA.
