# Results → paper map

This is the shortest route from a statement in the SIMBIOCHEM paper to its numbers. **Paper figure numbers** refer to the camera-ready workshop version. Some source CSVs retain names from an earlier figure sequence; the mapping below is intentional. Values in the repository's SVG gallery are read from these exact files, while the conceptual artwork is clearly marked.

| Camera-ready location | Question and panels | Numeric source | Repository view |
|:--|:--|:--|:--|
| Figure 1a–d | What is the fixed-geometry interaction-induced response, and why is a dipole insufficient? | [`Fig1_geometry.csv`](../data/source_data/Fig1_geometry.csv), [`Fig1_physical_response.csv`](../data/source_data/Fig1_physical_response.csv), [`Fig1_response_map.csv`](../data/source_data/Fig1_response_map.csv) | [Conceptual model diagram](../assets/architecture.svg); use the paper itself for the QM potential maps. |
| Figure 2a | Frozen features, starting response, learned corrections and constraints. | [Checkpoint manifest](../checkpoints/glider_site_response_ensemble.manifest.json), [frozen source](../evidence/frozen_source/prospective_1/README.md) | [Architecture](../assets/architecture.svg) (conceptual). |
| Figure 2b,c | Development controls and representation error floor. | [`Fig2_controls.csv`](../data/figure_data/Fig2_controls.csv), [`FigS2_means.csv`](../data/figure_data/FigS2_means.csv) | Data only; this repo does not duplicate the manuscript panel. |
| Figure 3a | Three prospective response-ESP panels, 56 solutes and 224 configurations. | [`Fig4_glider.csv`](../data/figure_data/Fig4_glider.csv), [`Fig4_mace_polar_l.csv`](../data/figure_data/Fig4_mace_polar_l.csv), other `Fig4_*` comparator tables | [Prospective transfer](../assets/prospective-transfer.svg). |
| Figure 3b | Paired GLIDER-minus-polar-baseline intervals and solute wins. | [`Fig4_paired.csv`](../data/figure_data/Fig4_paired.csv) | [Prospective transfer](../assets/prospective-transfer.svg); paired intervals remain in the CSV. |
| Figure 3c | Response-dipole performance as a secondary global observable. | `data/figure_data/Fig4_*.csv` | Data only; the paper shows all eligible models. |
| Figure 4a–c | Frozen coupling to W4: energy, torque and orientation. | [`PLOTTING_MANIFEST.json`](../data/figure_data/downstream/PLOTTING_MANIFEST.json), [`energy_pairs.csv`](../data/figure_data/downstream/energy_pairs.csv), [`torque_pairs.csv`](../data/figure_data/downstream/torque_pairs.csv), [`orientation_pairs.csv`](../data/figure_data/downstream/orientation_pairs.csv) | [Coupling overview](../assets/frozen-coupling.svg) shows energy; the paper carries torque and orientation panels. |
| Figure 4d | Same base response, independent waters W4–W12. | [`energy_rank.csv`](../data/figure_data/downstream/energy_rank.csv) | [Coupling overview](../assets/frozen-coupling.svg). |
| Figure 5a,b | Fixed-charge distance sweep and far-field recovery of higher exact moments. | [`FigS15_dense_energy.csv`](../data/figure_data/FigS15_dense_energy.csv), [`FigS15_dense_nrmse.csv`](../data/figure_data/FigS15_dense_nrmse.csv) | [Distance and moments](../assets/distance-and-moments.svg) shows MAE. |
| Table 1 / neighbour extension | Same frozen model on a separate 12-solute, three-species test. | [`Fig5_nonwater_effects.csv`](../data/figure_data/Fig5_nonwater_effects.csv), [`Fig5_nonwater_molecule_pairs.csv`](../data/figure_data/Fig5_nonwater_molecule_pairs.csv), [`Fig5_nonwater_per_case_metrics.csv`](../data/source_data/Fig5_nonwater_per_case_metrics.csv) | Tables. |

Supplementary figure and table values are under [`data/figure_data/`](../data/figure_data/README.md) and [`data/source_data/`](../data/source_data/README.md). The [source-data workbook](../data/source_data/SourceData.xlsx) makes the plotted values easier to inspect outside Python. The directory filenames are frozen source identifiers, not a claim that an earlier manuscript's figure numbering is still current.

## How to read the primary comparison

The three `Fig4_*` prospective tables each report an **equal-solute** aggregate. The per-solute inputs and original configuration-level archives are in [`data/prospective_1`](../data/prospective_1/README.md), [`data/prospective_2`](../data/prospective_2/README.md) and [`canonical_panel_3`](../canonical_panel_3/README.md). Paired 95% intervals resample whole solutes (100,000 replicates), preserving within-solute configurations. `MACE-L` is MACE-POLAR-1-L used as a public polar baseline without response-specific fitting. It is distinct from the QM reference.

The panel reductions are 38.6%, 40.1% and 32.3%. Their sum is **54/56 improved solutes**. These claims concern response potential, not uniform improvement in global dipoles. A 12-solute non-water extension is separately reported; it is not included in the 56-solute count.

## What the physical test can and cannot show

Figure 4 holds the solute-plus-W1–W3 base fixed and adds W4 only as an external electrostatic probe. The W4 error is in kcal mol⁻¹ for the **response-coupling contribution**. The exact-response-dipole control compresses the QM spatial response to a full vector at one predefined origin. It is not the total molecular dipole. The W4 result supports local reuse and near-field directionality; it does not measure total binding, solvation free energy, mutual polarization or molecular-dynamics performance. Figure 5 shows the expected qualification: higher exact moments recover in the far field.

## Frozen numerical outputs

[`prospective_1/`](prospective_1/README.md) and [`prospective_2/`](prospective_2/README.md) retain the original configuration, molecule and regime summary tables. Panel III's corresponding metrics sit with its independent-source [freeze record](../canonical_panel_3/README.md). The [integrity script](../scripts/reproduce/verify_release.py) checks archived predictions against their registries and physical constraints; the [statistics script](../scripts/reproduce/recompute_statistics.py) recomputes the aggregate and paired statistics from frozen per-solute tables.
