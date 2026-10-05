# Preparing an independent polarizable-physics comparison

No ARROW results are claimed in this release. The geometries and QM references needed to run that comparison are now included for Panels I–III and the nested 1/3/6/12-water test.

1. Read a configuration from the experiment ExtXYZ and its `n_solute_atoms` fragment boundary.
2. Evaluate the full complex and both isolated fragments at the same nuclear geometry with the same ARROW settings. Include the response of both solute and environment. The target is their complex-minus-fragments potential, not the total potential or only the solute induced dipole.
3. Sample that response on the exact `points_angstrom` array in the corresponding QM reference NPZ. Preserve units (Eh/e) and retain the full three-component response dipole (D).
4. Write method predictions indexed by the existing `config_id`, then apply the same configuration → solute → panel aggregation. Report any element/parameter coverage gaps explicitly, with a fixed coverage rule before examining errors.
5. If the model supplies site dipoles but no potential routine, its physical decoder, damping and units must be specified before comparing the near field. A far-field point-dipole substitute would test a different representation.

For the nested-size test, all 36 [geometries, references and frozen predictions](../experiments/shell_size/) are available. [Array dictionary](data_format.md). The public polar baseline remains a related pretrained model, so an independent physics baseline would add a different comparison rather than simply repeat it.
