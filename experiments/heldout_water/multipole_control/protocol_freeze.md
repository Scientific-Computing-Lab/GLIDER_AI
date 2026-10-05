# Frozen protocol: exact-QM global-multipole hierarchy

Frozen on 20 August 2026 before evaluating any quadrupole or octupole
downstream score. This is a post hoc control nested inside the already frozen
ten-solute confirmatory held-out-water experiment. It changes no configuration,
water ordering, QM density, model output, endpoint, or prospective result. The
machine-authoritative specification is `protocol_freeze.json`.

## Question and immutable inputs

The control asks whether the downstream difference between GLIDER and the
exact response dipole is recovered by adding the next two **global** moments of
the same exact QM base-response density. The base remains solute plus W1--W3,
and W4--W12 remain frozen external probes excluded from the base and model
input. All ten confirmatory solutes, 90 initial probe densities, 60 W4 torque
variants and 240 W4 orientations are required. No new base or probe QM will be
run, and no model will be retrained.

## Multipole convention

The stored response density matrix is a positive electron-number density
difference. Nuclear terms cancel in the counterpoise-consistent subtraction,
so the physical response charge density is its negative,
`rho_charge = -rho_electron`. About origin `O`, primitive moments are

`M^(n)_(i1...in) = integral rho_charge(r) product_k (r_ik-O_ik) dr`.

The dipole is `p_i=M^(1)_i`. The traceless Cartesian quadrupole and octupole
are

`Q_ij = 3 M^(2)_ij - delta_ij Tr(M^(2))`,

`O_ijk = 5 M^(3)_ijk - delta_ij T_k - delta_ik T_j - delta_jk T_i`,

with `T_i=sum_a M^(3)_iaa`. Their exterior potential is

`V(R) = q/R + p_i R_i/R^3 + Q_ij R_i R_j/(2 R^5)
        + O_ijk R_i R_j R_k/(2 R^7) + ...`.

The interaction with the complete frozen isolated-water QM potential `Phi_W`
is evaluated as

`E = q Phi_W(O) + p_i d_i Phi_W(O)
     + Q_ij d_ij Phi_W(O)/6 + O_ijk d_ijk Phi_W(O)/30`.

Moments come from analytic AO overlap, `r`, `rr` and `rrr` integrals. The
published exact-dipole coupling routine is called unchanged for the order-1
term. Second and third derivatives use a tensor-product seven-point central
stencil with primary spacing 0.01 bohr; 0.0075 and 0.015 bohr are fixed
convergence checks on methane and 2OJ9-original, selected as the smallest and
largest bases before outcomes.

## Origins and reporting

The three existing origins are unchanged: base mass-weighted COM, solute
mass-weighted COM, and base centre of nuclear charge. Every origin is retained.
There is no casewise origin selection. For each hierarchy order, a displayed
W4 energy control may use only the single panel-wide origin with lowest W4
MAE; torque analogously uses vector RMSE, and orientation uses mean absolute-
profile RMSE. The selected orientation origin is also used for the mean-centred
profile diagnostic. W4--W12 energy retains the origin chosen once by W4 energy
and cannot be reselected by rank or distance.

## Validation before interpretation

1. Physical charge and dipole signs are checked directly: `-Tr(D S)` is the
   physical response charge and `-Tr(D r)` must reproduce the stored response
   dipole within `1e-6 D`.
2. The order-1 field must reproduce every existing exact-dipole coupling used
   here within `1e-10 kcal/mol` because it invokes the same routine.
3. Direct moments at all three origins are checked against the exact primitive-
   moment translation identities through rank 3.
4. The existing independently recomputed rigidly rotated methane density checks
   `p'=R p`, `Q'=R Q R^T` and the rank-3 tensor rotation rule.
5. Primary/check finite-difference spacings are reported, never chosen from
   comparative performance.
6. All fixed cases and rigid water geometries remain accounted for.

## Exterior-convergence diagnostic

A finite global multipole expansion is formally an exterior/asymptotic
representation. Gaussian densities are noncompact, so no finite sphere
strictly encloses the response. For each origin, a level-4 Becke grid will
measure the radii containing 90%, 95% and 99% of the integrated absolute
response charge, together with the maximum base-nuclear radius and each probe's
oxygen/minimum-nuclear distance. Distance-to-extent ratios will be reported.
Poor performance for probes inside these practical extents will be interpreted
as a limitation of a low-order near-field truncation, never as a failure of
multipoles in general.

## Endpoints and statistics

The primary endpoints are W4 response-coupling energy, W4 response torque about
oxygen, and the complete 24-rotation absolute orientation profile. Mean-centred
profile shape and W4--W12 energy versus rank/distance are secondary. Existing
QM, GLIDER, MACE-POLAR-1-L and zero results are reused without recalculation.

Solute is the independent unit. Point estimates weight all ten solutes equally;
95% percentile intervals use the existing 100,000 PCG64 blocked resamples and
seed `2026081902`. Equal-family and leave-one-solute-out results are required.
No pass/fail threshold will be added after outcomes are seen.

The conclusion may say only how much information the tested low-order global
hierarchy recovers. It cannot rule out higher-order, distributed or alternative
compact representations, and the frozen response coupling remains one Coulomb
component rather than a complete intermolecular model.
