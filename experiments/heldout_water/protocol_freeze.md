# Frozen confirmatory protocol: downstream response breadth, range, and directionality

Frozen 19 August 2026 before any new confirmatory QM reference, GLIDER score,
or comparator score was evaluated. This is a post hoc confirmatory extension;
it does not alter the chronology or contents of the six-solute pilot or any
prospective GLIDER result. The machine-authoritative protocol is
`protocol_freeze.json`.

## Fixed panel and source rule

The full-workspace and Git-history audit found exactly ten new eligible neutral
solutes with completed explicit-water equilibration states: 2OJ9-original,
cdk2-1h1q, cdk2-28, ethane, ethanol, jnk1-17124, jnk1-18631, methane,
tyk2-ejm_42, and tyk2-ejm_45. Every eligible distinct identity is retained;
there is one pre-existing state per solute. Repeated test/alias copies and the
alternate 2OJ9 tautomer are identity-collapsed by the target-free rule stated
in the audit. The exact Git commits, paths, source hashes, state times, and
coordinates are in `configurations/configurations.csv` (SHA-256
`47a94ea95ce8ba3259a676c47e68643446a756d32fc80f8b2ed770468672f10e`).

These ten solutes arise from five broader source families. Solute is the
primary independent unit, while equal-family and leave-one-family-out analyses
will expose the related cdk2, jnk1, and tyk2 pairs. The original six pilot
solutes are excluded from confirmatory statistics and may appear only as a
separate pilot or explicitly exploratory combined summary.

Within each periodic parent state, waters are ordered by minimum-image oxygen
distance to the nearest solute atom, with zero-based source-water index as the
tie-break. The base is the solute plus W1--W3. W4--W12 are independent frozen
external probes and never enter the base response or either model input.
Parent salt ions and force-field lone-pair sites are excluded before ranking.
All ten bases and all 90 probes must remain accounted for.

## Frozen physics and models

The response target and electronic structure are unchanged: three
counterpoise-consistent restricted DF-RKS omegaB97X-D3(BJ)/def2-TZVPD
calculations at grid level 4, SCF tolerance `1e-10`, and 160-cycle maximum
produce `D_full - D_solute-ghost - D_waters-ghost`. Each external water uses
its isolated frozen QM density at the same level. Couplings are evaluated from
density-fitted AO Coulomb integrals, never by interpolation from the published
surface grid.

GLIDER uses the unchanged checkpoint SHA-256
`a288afa285128e13cb7a05ba9459f1dce07e79b8638c6e272c2509ae68d7a288`.
MACE-POLAR-1-L uses the official unchanged checkpoint SHA-256
`9f65f8dc6ddaff1d631e299cb531376a7da5e68d1bef04f34a2d5073d5ef114b`
and its normalized Gaussian `l<=1` density. No model is retrained. The exact
QM response dipole is evaluated at the base centre of mass, solute centre of
mass, and base centre of nuclear charge; the first is primary. All three are
shown, and the conservative best-origin summary selects one origin at panel
level, never a different origin case by case. Zero response is also retained.

For probe nuclei A and positive electron density `rho_j`, the measured scalar
is

`E_resp(B,Wj) = sum_A Z_A DeltaV_resp^B(R_A) - integral rho_j(r) DeltaV_resp^B(r) dr`.

This is only the frozen external-water Coulomb coupling to the electronic
response of the fixed base. It is not a full interaction energy, total force,
solvation energy, Hamiltonian, dynamics, or thermodynamic observable.

## Experiment A: independent-chemistry breadth

The primary W4 endpoints are response interaction energy, response torque
about W4 oxygen, and the complete deterministic orientation-energy profile.
Rigid-translation force is secondary. Energy is reported with equal-solute
MAE/RMSE, bias, correlations, every raw error, molecule-blocked uncertainty,
the fraction of solutes improved over MACE-L and the conservative fixed-origin
dipole comparison, and leave-one-solute-out/family sensitivity. Torque uses
vector and magnitude errors and directional metrics only above the fixed
`0.01 kcal/mol` reference threshold. Force uses the identical pilot central
difference (`0.002 Angstrom`, checks at `0.001/0.004 Angstrom`) and a
`0.05 kcal/mol/Angstrom` directional threshold.

W4 orientation uses the unchanged 24 proper rotations of SciPy's octahedral
group in canonical-quaternion lexicographic order, with oxygen fixed. Primary
profile errors cover the entire absolute and mean-centred profiles; correlation,
modulation error, and discrete extremum errors are secondary. Extremum angles
are interpreted only for QM modulation of at least `0.05 kcal/mol`.

## Experiment B: W4--W12 spatial range

The same fixed base response is coupled independently to W4, W5, ..., W12.
Energy error and paired advantage are plotted and tabulated against both rank
and actual nearest-solute oxygen distance. Fixed 0.5-A bins supplement the raw
continuous-distance display; binning will not be redesigned after outcomes.
Torque by rank/distance is secondary. The cumulative observable is the direct
sum `sum_{j=4}^k E_resp(B,Wj)` for every `k=4,...,12`; it remains a frozen
outer-environment response coupling, not a shell or solvation energy.

No unmeasured crossover will be extrapolated. A loss of clear full-field
advantage is described only within the observed distance/rank support and only
when effect size and blocked uncertainty support it.

## Experiment C: energy offset versus directionality

The same QM machinery also evaluates the full frozen base--water electrostatic
coupling. For W4, the primary ratios are response/full absolute energy,
response/full orientation modulation, and response/full torque magnitude.
W4--W12 energy and torque fractions are secondary. Raw numerators and
denominators are always retained. Ratios are undefined rather than inflated
when the full denominator is below `0.10 kcal/mol` (energy), `0.05 kcal/mol`
(orientation modulation), or `0.01 kcal/mol` (torque).

## Statistics, checks, and change control

Point estimates weight ten solutes equally. Percentile 95% intervals use
100,000 solute-blocked NumPy PCG64 resamples with seed `2026081902`.
Equal-family and leave-one-family-out results are sensitivity analyses.

The pilot's validated linearity, rigid-geometry, covariance, finite-difference,
quadrature, unit, and accounting checks are retained. New checks reconstruct
the exact dipole directly from the response density and verify near-zero net
response charge. All failures remain in the registry; no difficult or poorly
performing case may be removed.

There is no numerical editorial pass/fail gate. The verdict will use effect
size, cross-solute consistency, molecule-blocked uncertainty, family
sensitivity, numerical stability, and physical coherence. This protocol is
immutable after its freeze commit. A numerical bug may be corrected only in a
dated amendment that preserves the old rule and records the reason and all
affected outputs.

