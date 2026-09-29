# Pre-selection implementation erratum

Recorded at 2026-08-14T17:53:00Z, before the identity-only FreeSolv file was
extracted and before the selector was executed.

Static review after the first external preregistration timestamp found that a
Python list-membership expression compared dictionaries containing NumPy
arrays. Such comparison is undefined when more than one candidate has already
been accepted and would stop execution. The expression was replaced by exact
membership of the public FreeSolv ID string. This is a purely mechanical fix:
it does not change the preregistered eligibility rules, distances, strata,
round order, tie break, geometry preflight, inputs or outputs. No selected
identity, prediction or response label existed when the correction was made.

At 2026-08-14T17:55:00Z, the first selector invocation stopped before candidate
construction because the prospective 4--15-heavy-atom eligibility constraint
had inadvertently also been applied while loading the historical reference
set. One historical molecule has three heavy atoms. The loader was corrected
so that the size constraint applies only to Panel-III candidates; all 50
historical identities still participate in exact-overlap, similarity and
max--min calculations. The sanitized identity-only file had been created, but
no selected identity, eligible-pool file, geometry, prediction or response
label existed. This repair implements the written protocol and does not alter
any selection rule.
