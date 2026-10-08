# A held-out water tests the response field

Ten solutes supply a fixed base consisting of the solute and waters W1–W3. Waters W4–W12 probe that base independently. A probe enters neither the base QM calculation nor the model input. Its frozen charge distribution couples to the response field.

**Paper:** Fig. 4; Figs. S8–S12; Tables S15–S17; Appendix L.

| Open | Contents |
|:--|:--|
| [Geometries](geometries/) | Ten bases, W4 and all 90 outer waters. |
| [Source structures](source_structures/) | Original SDFs, solvated topology files and OpenMM equilibration states, with verified source commits and hashes. |
| [Solute identities](geometries/source_records.csv) | Canonical SMILES and exact source paths. `2OJ9-original` denotes the BMI ligand, not a QM protein complex. |
| [QM density matrices](references/) | Base complex/fragment densities and isolated probe states, including orientation variants. |
| [Predicted response sites](predictions/) | GLIDER and the polar baseline on the ten fixed bases. |
| [Per-probe results](results.csv) | Coupling to W4–W12. |
| [All physical evaluations](evaluated/) | Energy, torque, forces and orientation profiles. |
| [Analysis](analysis/) | Equal-solute aggregation, paired uncertainty, cumulative coupling and directionality fractions. |
| [Exact multipoles](multipole_control/) | Single-origin QM dipole, quadrupole and octupole controls. |

```bash
python -m pip install -e '.[qm]'
python scripts/reproduce/recompute_coupling.py --solute methane --water-rank 4
```

This evaluates Coulomb integrals from released density matrices on CPU. It does not run SCF. Larger solutes require more memory. The archived calculation used density fitting; the small CPU example reproduces its values to better than 5 × 10⁻¹¹ kcal mol⁻¹.

The measured energy is a response contribution to frozen Coulomb coupling. It is not total binding energy. An exact response dipole has the exact global vector, but represents its field at one origin. Near-field disagreement measures the limitations of that representation, not a defective QM dipole.

[Distance sweep](../distance_sweep/) · [Separate protocol pilot](../heldout_water_pilot/) · [Data formats](../../docs/data_format.md)
