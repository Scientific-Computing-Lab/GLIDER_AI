# GLIDER canonical-source prospective validation: preregistration

Frozen on 14 August 2026, before extraction or selection of any Panel-III
identity. This document specifies a one-shot external-validity test. It does
not authorize model fitting, calibration, geometry replacement after the
prediction freeze, or access to FreeSolv hydration properties.

## Scientific question

Does the frozen GLIDER response model retain its prospective advantage for
environment-induced response electrostatic potential (ESP) when solute
identities are selected label-blind from an independently established public
solute collection rather than from an author-designed chemistry pool?

## Immutable source and information boundary

The identity source is FreeSolv database v0.52, MobleyLab/FreeSolv commit
`6c7d19b4b565537365ffd22006aa2cd4643200c6`, file `database.txt`, expected
SHA-256 `2d13f095713bc39b85f85dd7b4e5483fbb12fc694bf253bb1d92a4c4d484f260`.
The extraction program reads only the first two semicolon-delimited fields:
FreeSolv ID and SMILES. It writes an identity-only intermediate table. No
hydration free energy, calculated free energy, uncertainty, compound name, or
other property field may be loaded or copied. The selector consumes only this
sanitized table.

## Eligibility and historical exclusion

Candidates must sanitize in RDKit, contain one covalent fragment, have formal
charge zero, contain 4--15 heavy atoms, contain only H, C, N, O, F, P, S, Cl,
Br or I, contain no isotope labels or radical electrons, and contain no
undefined bond. Exact canonical-isomeric-SMILES or connectivity-InChIKey
matches to any of the 50 historical GLIDER chemistries are excluded. A
candidate is also excluded if its radius-2, 2048-bit Morgan fingerprint has
Tanimoto similarity at least 0.70 to any historical chemistry. Construction
failure under the frozen conformer/environment program excludes the candidate
and is logged; the next deterministically ranked candidate is considered.
No rule may be relaxed.

## Strata and deterministic selection

Exactly 20 identities are selected, four in each mutually exclusive stratum.
Assignment precedence is: (4) any S or P; (5) otherwise any F, Cl, Br or I;
(3) otherwise any N or aromatic n/o/s atom; (2) otherwise any O; (1) all
remaining eligible structures. Stratum 1 has a hard quota of two structures
with nonzero aromatic fraction and two with zero aromatic fraction.

Descriptors are molecular weight, heavy-atom count, H-bond donors, H-bond
acceptors, TPSA, Crippen logP, ring count, rotatable-bond count, aromatic atom
fraction and heteroatom count. Each descriptor is robustly standardized by the
median and interquartile range of the complete eligible pool; a zero IQR is
replaced by one. Descriptor distance is
`min(||z_i-z_j||_2/sqrt(10), 1)`, so both distance components lie in [0,1].
The frozen pair distance is

`D(i,j) = 0.70 * (1 - Tanimoto(i,j)) + 0.30 * d_descriptor(i,j)`.

The score of a candidate is its minimum D to all 50 historical identities and
all already accepted Panel-III identities. Selection proceeds round-robin
through strata 1,2,3,4,5 for four rounds. The stratum-1 requested aromaticity
is aliphatic, aromatic, aliphatic, aromatic in rounds 1--4. Within the active
eligible subset the highest minimum-distance score wins. Exact numerical ties
are ordered lexicographically by SHA-256 of
`FreeSolv_ID|GLIDER_CANONICAL_PANEL_3_FREESOLV_V052_20260814`.

The selector is forbidden from using GLIDER predictions, comparator
predictions, model disagreement, response labels, hydration properties, or
human preference. Its source is frozen in this preregistration commit.

## Physical panel

For each identity, the frozen ETKDG/MMFF conformer construction and existing
three-water placement functions generate exactly four regimes: equilibrium,
cooperative chain, orientation perturbation and compressed contact. Seeds are
deterministic functions of the protocol-file SHA-256, source ID and canonical
SMILES. The task therefore changes identity provenance but not the physical
regimes. The complete panel is 20 molecule blocks and 80 configurations. All
identities and geometries are committed and pushed before prediction. No
geometry can be changed after this freeze.

## Frozen models and prediction chronology

GLIDER is the original response-supervised ensemble checkpoint
`a288afa285128e13cb7a05ba9459f1dce07e79b8638c6e272c2509ae68d7a288`.
Its MACE-POLAR-1-M feature encoder and M/L response prior, site charge
projection, site dipoles, ESP surface and dipole reconstruction remain
unchanged. Family C is not GLIDER.

The full comparator suite is MACE-POLAR-1-S/M/L with the audited official
Gaussian density decoder, AIMNet2's official default ensemble member,
MACE-MDP's official direct dipole, and GFN2-xTB through tblite 0.7.0, plus the
frozen zero response and relevant historical physical/same-data controls where
their frozen implementation and element coverage permit. Public checkpoint
selection is fixed by the preceding comparator audit, not Panel-III results.
Models lacking an observable or element receive no full-panel rank for that
endpoint; any supported-subset diagnostic is explicitly secondary.

All possible GLIDER and comparator predictions are generated, hashed,
committed and pushed before a reference directory is created. Prediction code
receives geometry, identity and fixed model parameters only. The final
pre-reference manifest must state `reference_response_labels_exist=false`.

## Reference response

Only after the remote prediction freeze, each fixed configuration is evaluated
at counterpoise-consistent omegaB97X-D3BJ/def2-TZVPD, grid level 4, SCF
convergence 1e-10. Full complex, ghost-basis solute and ghost-basis three-water
components use identical geometry and complete complex basis. The mutual
response is full minus solute minus waters. Primary references are response ESP
on the same deterministic probe surface and the induced full-complex dipole
vector in Debye.

The recovery ladder changes numerical solution only: (1) ordinary frozen SCF;
(2) restart/continuation with a converged same-system or same-chemistry density
and a larger iteration limit; (3) previously qualified zero-level-shift
second-order SCF continuation initialized from fragment densities. Functional,
basis, grid, tolerance, geometry and final zero-level-shift Hamiltonian cannot
change. Every action is logged. A case unresolved after this ladder remains in
the registry and makes the panel incomplete; it is never replaced.

## Metrics and inference

The primary endpoint is molecule-mean response-ESP NRMSE: per configuration,
RMSE over its frozen probe points divided by the reference response-ESP RMS;
configuration NRMSE values are averaged within each molecule and then equally
over molecules. The secondary endpoint is induced-dipole vector RMSE in Debye:
the square root of the mean squared Euclidean vector error over all 80
configurations. No endpoint, normalization or aggregation may change.

Paired inference uses 100,000 molecule-blocked bootstrap resamples with seed
20260814. For each method pair and endpoint, resampled molecule blocks are
aggregated exactly as above; reported intervals are the 2.5 and 97.5
percentiles of GLIDER minus comparator. Configuration resampling is not a
primary analysis. Molecule wins, regime values, medians and tails are
descriptive and cannot replace the preregistered aggregate.

## Interpretation

For ESP: **strong canonical replication** requires the lowest full-panel ESP
NRMSE among fair full-coverage executable comparators and a paired 95% interval
against the strongest comparator entirely below zero. **Nominal canonical
replication** requires the lowest observed ESP NRMSE with an interval crossing
zero. **Mixed result** applies when rank, molecule or regime evidence is
inconsistent with a clear replication. **Failed canonical replication**
applies when a fair comparator is clearly better or GLIDER exhibits a
catastrophic transferable failure. Dipole is ranked and its paired intervals
reported without an artificial threshold.

No experimental hydration or solvation value may be used for selection,
training, inference or scoring. Transfer language is limited to chemistries
absent from GLIDER response supervision; broad foundation pretraining may have
included exact identities or close analogues.
