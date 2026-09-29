# Probe sweep records

This folder records the placement grid and per-molecule values behind the distance/orientation study. [`probe_manifest.json`](probe_manifest.json) defines the accepted probe inventory; [`RAW_INPUT_HASHES.json`](RAW_INPUT_HASHES.json) locks the raw inputs. `w4_dense_positions.csv` and `w4_quantum_positions.csv` retain the two distinct probe models, while `orientation_couplings.csv`, `per_molecule.csv` and `summary.csv` hold measured or evaluated observables. The [paper map](../../../results/README.md) links the camera-ready figure to its plotted aggregates.

The fixed-charge sweep is a diagnostic of spatial representation as clearance grows. The quantum-water check is complementary; it should not be substituted for the fixed-charge panel or interpreted as self-consistent polarization.
