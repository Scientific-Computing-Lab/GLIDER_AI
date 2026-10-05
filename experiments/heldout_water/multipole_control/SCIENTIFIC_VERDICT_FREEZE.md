# Frozen scientific verdict: exact global multipole hierarchy control

Frozen after completion of the predeclared numerical validation and all comparative scoring, and before any manuscript edit.

## Verdict

The exact global quadrupole and octupole do not explain away GLIDER's advantage over the exact response dipole for the held-out W4 probe. Under the predeclared panel-wide origin rule, GLIDER has W4 energy MAE 0.0487 kcal/mol, compared with 0.149 for the exact dipole, 0.252 for dipole plus quadrupole, and 0.508 for dipole plus quadrupole plus octupole. The corresponding torque vector RMSEs are 0.0686, 0.961, 2.14, and 5.80 kcal/mol; the mean absolute orientation-profile RMSEs are 0.0496, 0.253, 0.394, and 0.812 kcal/mol. GLIDER's paired error is lower than every tested exact hierarchy with solute-blocked 95% intervals excluding zero for all three primary endpoints, and the result persists under equal-family and leave-one-solute-out summaries.

This result has an essential geometric qualification. A finite global multipole expansion is an exterior, asymptotic representation, whereas the held-out waters sample the response near field of extended solute--water complexes. At W4, the probe oxygen lies outside the radius containing 99% of the absolute response density for only 0/10 solute-COM, 1/10 base-COM, and 1/10 nuclear-centre expansions; the median oxygen-to-`r99` ratio is 0.65--0.67. Across W4--W12, increasing the order does not approach the exact coupling systematically, and most probes remain inside `r99`. The growing errors are therefore evidence that the tested low-order *global truncations* do not recover the finite-distance information in these configurations, not evidence that multipoles in general fail or that sufficiently high-order, distributed, or alternative compact representations could not reproduce the field.

## Numerical validation

- The physical response charge reconstructed from the positive electron-number response density is neutral to `6.67e-14 e`; the required electronic charge sign is applied explicitly.
- The density-derived dipole agrees with the stored exact response dipole to `4.76e-8 D`.
- The order-1 implementation reproduces every existing W4--W12 energy, W4 torque component, and W4 orientation energy to `4.44e-16 kcal/mol`.
- Direct moments at all three origins obey exact translation identities through rank 3 (maximum absolute discrepancies: `4.49e-13`, `5.14e-12`, and `8.26e-11` in the corresponding atomic units).
- Quadrupole and octupole traces vanish to `2.66e-15 e bohr^2` and `1.99e-13 e bohr^3`.
- Independently recomputed rotated densities reproduce tensor covariance to relative errors of `1.32e-4` for the dipole, `1.36e-4`--`9.51e-4` for the quadrupole, and `1.95e-4`--`6.22e-4` for the octupole.
- Changing the derivative stencil step from the primary 0.01 bohr to 0.0075 or 0.015 bohr changes quadrupole couplings by at most `2.96e-8 kcal/mol` and octupole couplings by at most `4.08e-5 kcal/mol` on the predefined smallest/largest validation cases.

## Frozen interpretation

Demonstrated: the exact first global response moment omits finite-distance energy and directional information that GLIDER recovers, and adding the exact global quadrupole or octupole in the tested near-field truncations does not recover that information.

Suggested: direct spatial prediction is more useful than low-order global compression for these neighbouring-water configurations.

Not demonstrated: that a converged multipole expansion, a distributed multipole model, or any other compact spatial representation must fail; that the full spatial field is uniquely necessary; or that these response couplings constitute a complete intermolecular Hamiltonian.

The Abstract should retain its exact-dipole statement. The main Results and Discussion may add one concise hierarchy result with its exterior-convergence qualification; the complete convention, validation, statistics, and figure belong in the Supplement.
