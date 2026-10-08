# Move the probe, keep the responding complex fixed

The ten held-out-water bases are reused with a probe at different directions, orientations and distances. The fixed-charge calculation evaluates the field at water charge sites. A smaller quantum-water calculation integrates against an isolated frozen water density.

**Where this appears:** In the SIMBIOCHEM workshop paper, Fig. 5a,b,
Fig. S13a,b, Table S20 and Appendix O. In the updated journal manuscript,
the same sweep is Supplementary Fig. S15, Table S19 and Section 14. The
journal figure draws its plotted values from
[workshop Fig. 5](../../figures/figure_05/) and
[workshop Fig. S13](../../figures/figure_S13/).

[Raw evaluations](raw/) contain the actual probe positions, potentials, fields, energies and torques. [Analysis](analysis/) contains position-level and solute-level tables. [Configurations](configurations/) supply directions, clearances and orientation identities, with repository-relative paths to the [base geometries, QM densities and predictions](../heldout_water/). Every referenced input is included.

```bash
python -m pip install -e '.[qm]'
python scripts/reproduce/run_distance_sweep.py \
  --config experiments/distance_sweep/configurations/A_base07.json \
  --mode dense --job-id reproduce_base07
```

This reruns the methane fixed-charge sweep on CPU and writes to `build/distance_sweep/results/`. The quantum mode also requires more Coulomb integrals and is slower. Raw archive filenames are direction/radius indices; the adjacent JSON gives their physical coordinates and accepted orientations. The [raw-input hash list](analysis/RAW_INPUT_HASHES.json) preserves the original archive identifiers, which map directly to filenames under `raw/`.

The abscissa is oxygen clearance from the **base van der Waals surface**, not solute–water separation. Moving an external probe does not test whether the inducing response disappears. That separate requirement is assessed in [dissociation](../dissociation/).
