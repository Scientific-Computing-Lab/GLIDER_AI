# Data formats and units

## Identity and geometry

`config_id` identifies a configuration and `molecule_id` identifies its solute. Multiple environments for one solute remain grouped for uncertainty estimates. ExtXYZ uses Å. Its metadata gives `n_solute_atoms`; those atoms come first, followed by the environment. Configuration registries supply SMILES, regimes, neighbour counts and source records where applicable.

## Response-reference NPZ

| Array | Meaning | Unit |
|:--|:--|:--|
| `points_angstrom` | Exact retained ESP probe positions, shape (N,3) | Å |
| `delta_esp_hartree_per_e` | Complex-minus-fragments QM response potential | Eh/e |
| `delta_dipole_debye` | Three Cartesian components of the response dipole | D |

`observable_registry.csv` maps configuration IDs to NPZ filenames. Some archives include fitted site labels and component results; inspect `numpy.load(path).files` for their fields. Fitted sites are not unique atomic populations.

## Prediction NPZ

| Array | Meaning | Unit |
|:--|:--|:--|
| `predicted_esp_hartree_per_e` | Potential at reference probe coordinates | Eh/e |
| `predicted_dipole_debye` | Molecular response dipole | D |
| `predicted_charges_e` | Atom-centred response charges | e |
| `predicted_dipoles_e_bohr` | Atom-centred response dipoles | e·bohr |

Dipole-only comparators legitimately lack ESP arrays. Registry paths resolve their outputs. A method is not assigned an ESP score when it does not supply that observable.

## QM matrices and physical reuse

Held-out-water `base_densities` stores complex and ghost-fragment component density matrices, their difference and the component observables. `probe_densities` stores `probe_density_matrix` and the isolated water geometry. Rotated variants store `density_matrix` and `positions_angstrom`. The decoder and Coulomb-integral code specify the AO basis and ordering. Energy is kcal/mol, torque is kcal/mol and force is kcal/mol/Å in the physical-result tables.

Dense distance-sweep NPZ files store actual water positions and charges, method potentials, fields, energies, forces and torques. Their adjacent JSON lists direction, clearance, orientation IDs and rejected/accepted placements.

## Aggregation

Response-ESP NRMSE is RMSE divided by the RMS QM response for that configuration. Average configurations within each solute, then weight solutes equally. Dipole RMSE is the square root of the equally aggregated Cartesian-component MSE. Paired uncertainty resamples solutes, not individual surface points.

The original 48 training environments provide 48 molecular observations, not tens of thousands of independent surface-point examples. [Executable metric definitions](../glider/metrics/core.py).
