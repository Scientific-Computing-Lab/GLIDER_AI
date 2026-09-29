# Verify and recompute

Run from the repository root after `pip install -e .`:

```bash
python scripts/reproduce/verify_release.py
python scripts/reproduce/recompute_statistics.py --root . --output /tmp/glider-statistics
python scripts/verify_companion.py
```

`verify_release.py` checks checkpoint and freeze hashes, all three prospective registries, total predicted response charge and dipole reconstruction from atomic sites. `recompute_statistics.py` starts from the frozen molecule-level tables and writes panel aggregates plus paired molecule-bootstrap statistics into the requested output directory. It does **not** rerun QM or revise the frozen model. `fetch_mace_polar.py` downloads official third-party MACE-POLAR-1 weights and checks their publisher-file SHA-256s; new-geometry inference requires M and L. The third-party files are excluded from this Git repository.
