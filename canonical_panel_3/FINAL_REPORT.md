# Canonical-source prospective Panel III: final report

## Terminal verdict

**A STRONG CANONICAL REPLICATION.**

GLIDER response-ESP NRMSE is **0.426581** (rank 1). The
strongest eligible comparator is MACE-POLAR-1-L at
**0.629706**. The GLIDER-minus-comparator difference is
-0.203126 (95% molecule-blocked bootstrap CI
[-0.247574, -0.153531]); GLIDER improves
19/20 molecules.

GLIDER induced-dipole RMSE is **0.226697 D** (rank
1). The strongest comparator is
MACE-MDP at **0.246295 D**. The
difference is -0.019598 D (95% CI
[-0.206149, 0.163604]); GLIDER improves
6/20 molecules.

## Integrity

The 20 identities were selected deterministically from an identity-only copy
of pinned FreeSolv v0.52. FreeSolv hydration values were not accessed. All 80
geometries, GLIDER predictions, and eligible comparator predictions were
hash-locked and pushed before reference QM began. No molecule or configuration
was removed or replaced. Four iodine-containing cases used the declared
same-Hamiltonian SCF recovery route; all completed at the frozen reference
theory.

## Scope

Panel III tests mutual electronic response of neutral organic solute--three-water
environments. FreeSolv supplies molecular identities only; this is not a
hydration-free-energy benchmark. Chemical identities are absent from GLIDER
response supervision, but exact identity or analogue exposure during broad
foundation pretraining cannot be excluded.
