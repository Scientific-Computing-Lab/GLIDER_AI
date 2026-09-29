# Prospective Panel I · 12 solutes / 48 configurations

[`configurations.extxyz`](configurations.extxyz) and [`configuration_registry.csv`](configuration_registry.csv) lock the solute–three-water geometries. Each configuration has a response observable in [`references/`](references/README.md), and archived model outputs in [`predictions/`](predictions/README.md). The `glider` subfolder is the frozen candidate scored in the camera-ready paper. The [panel's molecule-level results](../../results/prospective_1/README.md) and the [paper map](../../results/README.md) explain the aggregate.

Paper Figure 3 reports GLIDER response-ESP NRMSE **0.322** versus **0.524** for the response-unfitted MACE-L polar baseline, with **12/12** solutes improved. Panel I did *not* pass a separate preregistered **joint** gate requiring both ESP and dipole improvement. The panel was retained intact. Its solutes were excluded from GLIDER response fitting, not proven absent from broad MACE pretraining.
