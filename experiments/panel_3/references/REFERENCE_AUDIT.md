# Panel-III reference audit

**Result: PASS.** All 80 preregistered configurations have complete, finite
counterpoise-consistent response references. Probe coordinates match the
prediction-frozen GLIDER files exactly; no molecule or configuration was
excluded, replaced, or altered.

## Chronology and integrity

- Prediction-freeze commit: `f9a87b1d5a91a45f595ccd1891f59ec8d9445fea` (pushed before QM).
- Pre-QM integrity commit: `5347a97c8f01a998e736f9b453138d9e73419a1f`.
- Frozen prediction-manifest SHA-256: `3dc5641a289f48118e0e372c341286573d4a381a2cf2ea9b7ae2b2b80ec8c73e`.
- All 80 reference files postdate the remote prediction freeze.
- FreeSolv hydration quantities were not accessed.

## Reference protocol

All components use DF-RKS omegaB97X-D3(BJ)/def2-TZVPD, grid level 4 and
SCF convergence 1e-10 at the frozen geometry and complete complex basis.
The mutual response is the full complex minus ghost-basis solute minus
ghost-basis three-water cluster.

Seventy-six cases converged by the ordinary route. The four iodine-containing
environments for `mobley_2727678` used the preregistered same-Hamiltonian
recovery ladder. Final solutions use zero-level-shift second-order/DIIS SCF;
functional, basis, grid, geometry, convergence threshold and Hamiltonian were
unchanged. All three components converged in every case. The exact routes and
hashes are in `QM_RECOVERY_LOG.csv`.

## Superseded raw metadata

`references/manifest.json` is retained byte-for-byte because it was emitted by
the pre-existing QM generator. Its inherited development-era scope and
`prospective_labels_accessed=false` field are stale for Panel III. The audited
`REFERENCE_MANIFEST.json`, `REFERENCE_FILE_AUDIT.csv`, and this report are the
authoritative Panel-III metadata; no numerical reference file was changed.
