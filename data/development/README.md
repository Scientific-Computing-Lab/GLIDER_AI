# Sparse-supervision development set

The response head was fitted with **48 three-water environments from 14 chemistries**. `combined_development.extxyz` stores the geometries; `observables/` stores the QM response labels keyed by the observable registry; `frozen_training/` stores the already-extracted MACE-M features and averaged MACE-M/L starting response used for the frozen fit.

The camera-ready paper's Figure 2b,c uses leave-one-chemistry-out development controls; the compact plotted values are at [`figure_data/Fig2_controls.csv`](../figure_data/Fig2_controls.csv). These are development diagnostics, not prospective test results. No energy or force labels trained GLIDER. The frozen checkpoint and its exact training manifest are in [`checkpoints/`](../../checkpoints/README.md).
