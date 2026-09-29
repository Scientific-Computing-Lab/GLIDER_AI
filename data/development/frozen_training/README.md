# Frozen response-head training inputs

These two archives are the byte-exact cached inputs used to train the published
GLIDER response head. They were recovered from the preserved development tree
and match the hashes recorded before prospective scoring.

- `mace_polar_m_features.npz`: frozen MACE-POLAR-1-M product-basis features for
  the 48 response-supervised environments.
- `mace_polar_ml_average.npz`: the preregistered, unfitted equal-weight
  MACE-POLAR-1-M/L site-response prior for the same environments.

The raw geometries and counterpoise-consistent response targets are in the
parent development directory. Run `make verify-training-inputs` to validate the
complete training tuple and `make retrain` to repeat response-head fitting.

The cached features are derivative numerical outputs of the official
MACE-POLAR-1 checkpoint. The upstream checkpoint is not redistributed; its
license and release terms remain applicable.
