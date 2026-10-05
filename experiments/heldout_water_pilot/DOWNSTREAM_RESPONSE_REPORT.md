# Downstream response coupling to a held-out water

## Verdict

Yes, with an important statistical qualification. The base-cluster response
has a measurable Coulomb coupling to the held-out fourth water (QM RMS 0.287
kcal/mol), and GLIDER predicts its aggregate energy, torque, and
orientation-dependent profile more accurately than MACE-POLAR-1-L and than a
field built from the exact QM response dipole alone. GLIDER reduces aggregate
energy MAE from 0.199 to 0.0965 kcal/mol relative to MACE and from 0.173 to
0.0965 kcal/mol relative to the best of the three predeclared dipole origins.
Its energy advantage over MACE is supported by the solute-blocked interval,
whereas the interval relative to the exact-dipole controls includes zero and
the per-solute energy advantage occurs in only three of six solutes. The
torque and deterministic orientation results provide stronger evidence for
useful higher spatial structure, but the independent sample remains six
solutes. Under the frozen decision rules, the recommended manuscript action
is **SI ONLY**, not main text.

This is a post hoc diagnostic. It does not alter the chronology or artifacts
of any existing GLIDER result, and no manuscript file was edited.

## Frozen protocol

The protocol was committed as `db534143233ba05bf8cb76410aaf126309594a23`
before a downstream score was generated. The immutable human- and
machine-readable specifications are `protocol_freeze.md` and
`protocol_freeze.json`.

- Cases: all 24 pre-existing liquid-bridge parents, comprising six solutes and
  four fixed snapshots per solute.
- Base: solute plus W1--W3, the three nearest waters under the repository's
  deterministic minimum-image oxygen-to-nearest-solute-atom ordering.
- Probe: W4, the fourth-nearest water from the same parent. W4 is absent from
  both model inputs and all three base-response calculations.
- Models: unchanged GLIDER and MACE-POLAR-1-L predictions; neither was
  retrained. Exact-dipole-only and zero-response controls were evaluated.
- Statistical unit: solute. Intervals use 100,000 solute-blocked bootstrap
  resamples with NumPy PCG64 seed 2026081901. All 24 paired points and all six
  leave-one-solute-out estimates are retained.
- Exact-dipole origins: base centre of mass (primary), solute centre of mass,
  and base centre of nuclear charge. None was selected after seeing errors.

The exact case identities, trajectory/topology hashes, source frame indices,
W1--W4 source indices, and distances are in
`configurations/configurations.csv`.

## Physical definition

For a frozen base complex B = solute + W1 + W2 + W3 and the isolated quantum
charge distribution of W4,

\[
E_{\rm resp}(B,W4) = \sum_A Z_A\,\Delta V_{\rm resp}^B(R_A)
- \int \rho_{e,W4}(r)\,\Delta V_{\rm resp}^B(r)\,dr .
\]

The base response is the counterpoise-consistent density difference from
three DF-RKS omegaB97X-D3(BJ)/def2-TZVPD calculations at fixed coordinates
and complete base basis. W4 is computed in isolation at the same functional,
basis, grid, and SCF tolerance. The coupling is evaluated from AO density
matrices with density-fitted Coulomb integrals, not from interpolation on the
published surface grid. W4 does not polarize the base.

This quantity is only the response-electrostatic contribution experienced by
the frozen external water. It is not a four-water interaction energy, a total
intermolecular energy, a solvation energy, a total force or torque, or a
dynamics result.

## Experiment 1: response interaction energy

### Signal magnitude

Across the 24 cases, QM `E_resp` has:

- RMS magnitude: 0.287 kcal/mol;
- median absolute magnitude: 0.0645 kcal/mol;
- range: -0.862 to +0.516 kcal/mol;
- RMS magnitude: 0.485 `kBT` at 298.15 K;
- median `|E_resp|/|E_frozen electrostatic|`: 3.93%.

The signal passes the predeclared physical-scale threshold (RMS at least 0.10
kcal/mol and at least ten times the numerical uncertainty). It is
heterogeneous: the median case is much smaller than the RMS, and the ratio to
the frozen electrostatic interaction is only contextual because that
denominator is not a total interaction energy.

### Aggregate prediction

| Method | MAE (kcal/mol) | RMSE (kcal/mol) | Bias (kcal/mol) | Pearson r | Spearman rho |
|---|---:|---:|---:|---:|---:|
| GLIDER | 0.0965 | 0.157 | +0.0610 | 0.862 | 0.725 |
| MACE-POLAR-1-L | 0.199 | 0.321 | -0.0738 | 0.325 | 0.324 |
| Exact dipole, base COM | 0.178 | 0.266 | +0.00214 | 0.537 | 0.563 |
| Exact dipole, solute COM | 0.173 | 0.276 | +0.0230 | 0.478 | 0.401 |
| Exact dipole, nuclear centre | 0.177 | 0.263 | +0.000068 | 0.548 | 0.574 |
| Zero response | 0.185 | 0.287 | +0.0294 | undefined | undefined |

Thus GLIDER's aggregate MAE is 51.5% lower than MACE and 44.2% lower than
the best aggregate exact-dipole MAE. The conclusion does not depend on a
favourable dipole origin: all three predeclared origins give MAE 0.173--0.178
kcal/mol.

The solute-blocked 95% interval for paired GLIDER-minus-MACE MAE is -0.210 to
-0.0174 kcal/mol. The corresponding interval relative to the best-MAE dipole
origin is -0.197 to +0.0350 kcal/mol and therefore includes zero. GLIDER has
lower per-solute MAE than MACE in five of six solutes, but lower MAE than the
globally best dipole origin in only three of six. The GLIDER leave-one-solute-
out MAE spans 0.0847--0.111 kcal/mol; the corresponding ranges are
0.152--0.226 for MACE and 0.126--0.199 for the best-MAE dipole origin.

The energy result therefore supports an aggregate downstream advantage and a
robust advantage over MACE, but it does not establish a solute-general energy
advantage over an exact dipole with only six independent solutes.

## Experiment 2: response force, torque, and orientation

### Rigid-translation response force

The QM response-force RMS magnitude is 0.420 kcal/mol/Angstrom and its median
magnitude is 0.157 kcal/mol/Angstrom. GLIDER's vector RMSE is 0.317
kcal/mol/Angstrom, compared with 0.341 for MACE and 0.498 for the best
exact-dipole origin. This is only a 7.1% reduction relative to MACE, below the
predeclared 15% derivative threshold. Magnitude MAE is essentially tied
(0.112 for GLIDER, 0.113 for MACE). For the 17 cases above the frozen force
threshold, median angular error is 30.7 degrees for GLIDER, 64.4 degrees for
MACE, and 45.3 degrees for the primary exact-dipole field.

The force result is therefore mixed: GLIDER improves direction and modestly
improves vector error, but it does not provide a strong force-magnitude result
against MACE.

### Response torque about W4 oxygen

The QM response-torque RMS magnitude is 0.577 kcal/mol per radian
(dimensionally kcal/mol) and its median magnitude is 0.272. GLIDER's vector
RMSE is 0.151 kcal/mol, compared with 0.317 for MACE and 0.739 for the best
exact-dipole origin. These are reductions of 52.3% and 79.5%, respectively.
Median angular error over all 24 meaningful cases is 15.0 degrees for GLIDER,
22.1 degrees for MACE, and 29.5 degrees for the primary exact dipole.

The solute-blocked paired GLIDER-minus-MACE torque vector-RMSE interval is
-0.315 to -0.0382 kcal/mol. GLIDER is better than MACE in five of six solutes;
the exception is the low-torque heteroaromatic case. The torque result is the
clearest derivative evidence that the learned spatial field changes what the
held-out water feels beyond a global molecular dipole.

### Deterministic orientation profiles

The orientation test uses the fixed parent-00 case from every solute and all
24 proper octahedral rotations with the W4 oxygen fixed. The primary metric
is error in the complete absolute energy profile.

| Method | Mean profile RMSE (kcal/mol) | RMS profile error (kcal/mol) | Mean absolute amplitude error (kcal/mol) | Median Pearson r |
|---|---:|---:|---:|---:|
| GLIDER | 0.0611 | 0.0728 | 0.0437 | 0.974 |
| MACE-POLAR-1-L | 0.180 | 0.235 | 0.202 | 0.869 |
| Exact dipole, best mean RMSE origin | 0.183 | 0.225 | 0.273 | 0.777 |
| Zero response | 0.178 | 0.207 | 0.612 | undefined |

GLIDER has lower profile RMSE than MACE in all six fixed cases. The
solute-blocked GLIDER-minus-MACE mean-profile-RMSE interval is -0.252 to
-0.0448 kcal/mol; the interval against the best dipole origin is -0.232 to
-0.0252 kcal/mol. Discrete extremum orientations are not emphasized: the
24-member grid, near-degenerate extrema, and one low-modulation case make
those labels unstable even when the whole profile agrees well.

## Numerical and physical checks

All predeclared checks passed. Full values are in
`sanity/numerical_checks.json`.

- All 24 cases are present; no case was removed. Each case has four distinct,
  monotonically ordered water indices, and no selection used model or QM
  performance.
- All 96 base-component/initial-probe SCFs and all 354 newly computed rigid-
  probe SCFs converged. Six identity orientations reused their already
  computed isolated-W4 density.
- Maximum component electron-count error is `2.09e-12 e`; maximum response
  net electron count is `4.54e-14 e`; maximum W4 electron-count error is
  `3.91e-14 e`.
- Maximum linearity error in the response coupling is `1.30e-15 Eh`.
- Direct GLIDER site evaluation reproduces stored GLIDER surface values to
  `1.39e-17 Eh/e`.
- GLIDER site-gradient convergence changes energy by at most `1.03e-9`
  kcal/mol; MACE quadrature levels 4--6 change it by at most `3.84e-8`
  kcal/mol.
- A common translation changes any scalar coupling by at most `3.16e-9`
  kcal/mol. A 37-degree common rotation with independently recomputed QM
  densities changes the QM scalar by `1.28e-5` kcal/mol.
- An independent repeat of the least stable SCF case changes `E_resp` by
  `1.50e-5` kcal/mol. This is treated as the practical scalar numerical floor
  and is more than four orders of magnitude below the QM signal RMS.
- Maximum force-vector change over 0.001/0.002/0.004-Angstrom steps is
  `2.75e-6` kcal/mol/Angstrom; maximum torque-vector change over
  0.05/0.1/0.2-degree steps is `8.28e-6` kcal/mol.
- Rigid rotations preserve bond lengths to `3.22e-15` Angstrom and the water
  angle to `2.13e-13` degrees.

## Limitations

The inference unit is only six solutes, each represented by four correlated
snapshots. The W4 probe is frozen and quantum-mechanical but does not polarize
the base or respond self-consistently. The experiment isolates one Coulomb
component and omits exchange, dispersion, nuclear relaxation, mutual
polarization, and all other terms of a full interaction. Forces and torques
are derivatives only of this response-electrostatic component. Orientation
profiles cover one predeclared snapshot per solute and a discrete 24-member
rotation grid. These boundaries prevent interpreting the results as total
solvent forces, stable dynamics, or thermodynamic accuracy.

## Frozen decision and recommended action

The signal criterion passes. GLIDER exceeds the 20% aggregate energy-
improvement threshold against MACE and every exact-dipole origin, and the
torque/orientation criteria pass. However, the frozen main-text energy rule
also required a lower per-solute energy MAE in at least four of six solutes
against the best exact-dipole control; the observed count is three. The
predeclared **MAIN TEXT** rule therefore does not pass. The **SI ONLY** rule
does pass because aggregate energy MAE improves by more than 15% against both
baselines and the signal is nontrivial. This recommendation is driven by the
protocol, not by editorial preference after seeing the result.

## Manuscript proposals (not implemented)

**Strongest justified one-sentence claim**

> Across 24 liquid-derived complexes, GLIDER's spatial response reduced the
> aggregate frozen-W4 response-coupling energy MAE by 52% relative to
> MACE-POLAR-1-L and reproduced the response torque and orientation profile
> more accurately than either MACE-POLAR-1-L or a field built from the exact
> QM response dipole alone.

**Possible abstract sentence**

> In a post hoc liquid-derived diagnostic, coupling the predicted response to
> the isolated quantum charge density of the next water lowered aggregate
> energy error by 52% relative to MACE-POLAR-1-L and recovered response torque
> and orientation structure more accurately than an exact-dipole-only field.

**Possible main-text paragraph**

> We next asked whether the field improvement changes a physical quantity
> experienced by a neighbouring molecule. In a post hoc diagnostic over all
> 24 liquid-derived parents, we held out the fourth-nearest water and coupled
> its frozen isolated-QM charge density to the response field of the solute
> plus three-water base. The QM coupling had an RMS magnitude of 0.287
> kcal/mol. GLIDER gave an energy MAE of 0.0965 kcal/mol, compared with 0.199
> for MACE-POLAR-1-L and 0.173 for the best of three predeclared fields built
> from the exact response dipole. It also reduced response-torque vector RMSE
> from 0.317 to 0.151 kcal/mol relative to MACE and reproduced six
> deterministic orientation profiles with a mean RMSE of 0.0611 kcal/mol.
> The energy advantage over the exact-dipole control was not uniform across
> the six solutes, but the torque and orientation results show that GLIDER's
> spatial response contains downstream information that is lost in a single
> global moment.

**Possible figure layout**

One four-panel SI figure: (a) predicted versus QM `E_resp` with identity and
zero lines; (b) paired per-case energy absolute errors; (c) force and torque
vector RMSE with the angular summaries in the caption; (d) the six
orientation profiles as small multiples, with a compact inset showing mean
profile RMSE. The full per-case points, blocked intervals, and origin
sensitivity remain in accompanying SI tables.

