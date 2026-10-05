# Frozen protocol: downstream response coupling to a held-out water

Frozen 19 August 2026, before any downstream interaction-energy, force,
torque, orientation, or model-comparison score was generated. This is a post
hoc diagnostic. It does not alter or extend the chronology of any existing
GLIDER result.

## Scientific question and boundary

For each of the complete 24-case liquid bridge, the base is the solute plus
the three nearest waters and the external probe is the fourth-nearest water
from the same parent snapshot. The probe water is excluded from the GLIDER
input, the comparator input, and the base response definition. The measured
quantity is only the Coulomb coupling between the frozen probe-water charge
distribution and the electronic redistribution of the fixed base complex. It
is not a four-water interaction energy, a total force, or a solvation energy.

No model is retrained. No case may be removed or replaced on the basis of a
prediction, reference value, convergence difficulty, or error. A failed case
remains in the accounting and is reported as failed.

## Fixed configurations

The authoritative case list is
`configurations/configurations.csv` (SHA-256
`1afebf4e1d4b74b49a323630adb41bc4e7b1e6a95fe8184c43670ad705c2040b`).
It contains all six liquid-bridge solutes and all four pre-existing parent
snapshots per solute. The base and probe coordinates are in
`configurations/base_3water.extxyz` and
`configurations/heldout_w4.extxyz`; the combined file is for coupling
bookkeeping only. Water ordering is minimum-image oxygen distance to the
nearest solute atom, with source water index as the tie-break. W1--W3 form the
base; W4 is held out. Reconstruction of the pre-existing base geometries has a
maximum absolute discrepancy of `4.999e-9` Angstrom.

## Frozen models

- GLIDER: unchanged checkpoint
  `checkpoints/glider_site_response_ensemble.pt`, SHA-256
  `a288afa285128e13cb7a05ba9459f1dce07e79b8638c6e272c2509ae68d7a288`.
  The existing 24 liquid-bridge site predictions are used. Its response is the
  released atom-centred point monopole/dipole representation.
- MACE-POLAR-1-L: unchanged official checkpoint, SHA-256
  `9f65f8dc6ddaff1d631e299cb531376a7da5e68d1bef04f34a2d5073d5ef114b`.
  The existing frozen 24-case multipoles are used with the official normalized
  Gaussian `l<=1` density, smearing width 1.5 Angstrom. Point decoding is not
  substituted.
- Exact-dipole-only: the exact QM base response dipole for that case, converted
  to a single point dipole. The primary expansion origin is the mass-weighted
  centre of mass of the complete base. Predefined sensitivity origins are the
  solute centre of mass and the base centre of nuclear charge. No origin is
  selected by performance.
- Zero response: identically zero energy, force, and torque.

MACE-POLAR-1-M/S are not included: they are secondary models and are not
needed to answer the predeclared GLIDER-versus-best-L comparison.

## Reference electronic structure and frozen probe

The base response uses the manuscript Hamiltonian: restricted density-fitted
RKS omegaB97X-D3(BJ)/def2-TZVPD, PySCF grid level 4, SCF tolerance `1e-10`,
maximum 160 cycles, and the three counterpoise-consistent calculations in
identical coordinates and the complete base basis:

1. full solute + W1 + W2 + W3;
2. solute with W1--W3 ghost basis;
3. W1--W3 with solute ghost basis.

The response density matrix is `D_full - D_solute_ghost - D_waters_ghost`.
The isolated W4 density is computed at the same functional, basis, grid, and
SCF tolerance at its frozen trajectory geometry, without base ghost centres.
It never polarizes the base. The historical surface arrays are used only to
verify the recomputed response dipole and field, never to evaluate the new
coupling.

For probe nuclei `A` and positive electron number density `rho_W4`,

`E_resp = sum_A Z_A DeltaV_resp(R_A) - integral rho_W4(r) DeltaV_resp(r) dr`.

The primary QM value is evaluated from the response and W4 AO density matrices
with density-fitted Coulomb integrals, consistent with the reference SCF; it
is not interpolated from a surface grid. The same machinery evaluates the
frozen electrostatic interaction with the full base density as a secondary
scale comparison.

GLIDER point multipoles are coupled to W4 through the isolated-W4 potential
and its gradient at every GLIDER site. The gradient uses a symmetric
`1e-4`-bohr displacement and is checked at `5e-5` and `2e-4` bohr. As an
independent check, the resulting coupling is compared with direct level-5
quadrature of the GLIDER field over W4.

MACE-POLAR's smooth Gaussian response potential is integrated directly over
the isolated-W4 electron density using a PySCF level-5 atom-centred quadrature,
plus the exact nuclear term. Levels 4 and 6 are the predeclared convergence
check on the six parent-00 cases. The level is not chosen by model error.

Units are Hartree internally, converted with
`627.5094740631 kcal mol^-1 Hartree^-1`; coordinates use
`0.529177210903 Angstrom bohr^-1`. Room-temperature context uses
`kBT = 0.592484949 kcal/mol` at 298.15 K.

## Experiment 1: energy

Every case is evaluated for QM, GLIDER, MACE-POLAR-1-L, each of the three
exact-dipole origins, and zero response. Primary metrics are MAE, RMSE, signed
bias, Pearson correlation, Spearman correlation, and all per-case absolute
errors in kcal/mol. Signal summaries are RMS `E_resp^QM`, median
`|E_resp^QM|`, range, and the secondary response-to-full-frozen-electrostatic
scale ratio.

The independent statistical unit is the solute. Point estimates give equal
weight to the six solutes (four snapshots each). Percentile 95% intervals use
100,000 molecule-blocked bootstrap resamples with NumPy PCG64 seed
`2026081901`. Paired GLIDER-minus-baseline MAE and RMSE are computed inside
each resample. All six leave-one-solute-out estimates are reported. Undefined
zero-baseline correlations remain undefined.

Planned figures are (1) predicted-versus-QM energy with the zero line and
identity line and (2) paired per-case absolute errors for GLIDER, MACE-L, and
the primary exact-dipole origin. Planned tables are the complete per-case
values, per-solute metrics, aggregate metrics, blocked intervals, and
leave-one-solute-out sensitivity.

## Experiment 2: response force and torque

The response force is the negative rigid-translation derivative of this
coupling with the base fixed. It is evaluated by three-axis central finite
differences with primary step `0.002 Angstrom`. Convergence steps are `0.001`
and `0.004 Angstrom` on the six fixed parent-00 cases.

The response torque is about the held-out water oxygen. It is evaluated by
central rigid rotations about the fixed laboratory x, y, and z axes with
primary step `0.1 degree`; checks use `0.05` and `0.2 degree` on the same six
parent-00 cases. Torque is reported in kcal/mol/radian (dimensionally
kcal/mol). Water internal geometry is invariant under every transformation.

Force metrics are vector RMSE (`sqrt(mean ||F_pred-F_QM||^2)`), magnitude MAE
and RMSE, cosine similarity, and angular error. Angular force metrics are
reported only when `||F_QM|| >= 0.05 kcal/mol/Angstrom` and the prediction is
nonzero. Torque uses the analogous vector/magnitude metrics; angular metrics
require `||tau_QM|| >= 0.01 kcal/mol` and nonzero prediction. Thresholded case
counts are always shown. Zero response has well-defined vector and magnitude
errors but no direction.

## Deterministic orientation test

Orientation profiles use exactly the parent-00 case from each of the six
solutes, selected without reading a response or model value. The W4 oxygen is
fixed. Its rigid geometry is acted on by the 24 proper rotations of SciPy's
octahedral group `Rotation.create_group("O")`, with rotations placed in stable
lexicographic order of canonicalized quaternions. There is no random seed.

The primary profile metric is RMSE over all 24 absolute orientation energies.
Secondary metrics are Pearson/Spearman profile correlation, modulation
amplitude error, and discrete minimum/maximum orientation geodesic error.
Extremum orientation errors are interpreted only if the QM modulation is at
least `0.05 kcal/mol`; otherwise they are reported as numerically
uninformative. The planned figure contains six small-multiple profiles and a
summary panel of profile RMSE and modulation amplitude.

## Numerical checks fixed before scoring

1. Coupling is linear in the three component density matrices.
2. A common translation by `[1.234, -0.731, 0.419]` Angstrom leaves scalar
   energy invariant. For the first lexicographic case, a common 37-degree
   rotation about normalized `[1, 2, 3]` is also checked by recomputing the
   rotated isolated/base densities.
3. Direct GLIDER site evaluation reproduces every existing stored GLIDER
   surface field to `1e-8 Hartree/e` or better.
4. Force and torque step convergence is reported before derivative results are
   interpreted.
5. All rigid transformations preserve both O-H lengths and the H-O-H angle to
   `1e-10 Angstrom`/`1e-10 degree`.
6. Recomputed exact response dipoles must agree with the stored liquid-bridge
   values within `1e-5 D`; otherwise scoring stops for diagnosis.
7. Every model site count and coordinate is checked against the frozen base;
   W4 coordinates may appear only in the coupling evaluator.
8. The registry must contain exactly 24 cases, six solutes, four snapshots per
   solute, with no omissions.
9. Electron counts, net response charge, Hartree conversions, and force units
   are checked independently.
10. Near-zero reference force/torque follows the fixed angular thresholds and
    is not forced into a directional statistic.

## Frozen interpretation and action criteria

The signal is called physically nontrivial only if RMS `|E_resp^QM|` is at
least `0.10 kcal/mol` and at least ten times the largest estimated numerical
uncertainty. A material energy improvement means GLIDER equal-solute MAE is at
least 20% lower than both MACE-L and the best of the three predeclared exact-
dipole origins, with GLIDER lower per-solute MAE in at least four of six
solutes. Broad intervals remain visible and are never converted into a
significance claim.

- Recommend **MAIN TEXT** only if the signal and material energy criteria hold
  and GLIDER also lowers either force vector RMSE or orientation-profile RMSE
  by at least 15% against both MACE-L and the best exact-dipole origin.
- Recommend **SI ONLY** if the signal is nontrivial and GLIDER lowers energy
  MAE by at least 15% against both baselines, but derivative/orientation
  evidence does not meet the main-text criterion.
- Recommend **DO NOT USE** if the response signal is below the fixed physical/
  numerical threshold, GLIDER does not improve materially over the relevant
  baselines, or the exact dipole is within 10% of GLIDER's energy MAE without
  a derivative/orientation advantage.

These thresholds govern interpretation only. They do not govern inclusion,
settings, case selection, or whether an unfavorable result is retained.

## Change control

The exact-QM-density implementation is mandatory for the primary experiment.
If it proves technically invalid, scoring stops and the blocker is documented
before any fixed-charge fallback is considered. This file and
`protocol_freeze.json` are immutable after the freeze commit. A numerical bug
may be corrected only through a separately dated amendment that records the
old rule, new rule, affected outputs, and reason.

