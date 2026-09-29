# Physical reconstruction utilities

[`observables.py`](observables.py) reconstructs a spatial ESP and molecular dipole from atom-centred response charges and dipoles, with explicit Bohr/Å/Debye conversions. [`surface.py`](surface.py) defines the deterministic molecular sample surface used by the archived predictor. The [release verifier](../../scripts/reproduce/verify_release.py) checks net response charge and the reconstructed dipole against frozen per-configuration predictions.
