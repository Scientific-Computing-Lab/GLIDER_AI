# Diagnostics and geometry follow-ups

The published checkpoint is fixed. These scripts ask where its learned response succeeds or fails; they do not silently revise the training set or the original test panels. See each [experiment folder](../../experiments/) for the released inputs and results.

## Non-water contact follow-up

The original [72-case non-water panel](../../experiments/nonwater/) used frozen replacements at water-oxygen anchors. Its response scores are reproducible, but a later QM energy audit found most contacts repulsive. The scripts here make a **separate** contact-geometry follow-up from the original starting poses. This follow-up was designed after the audit and must not be described as an independent preregistered replication.

The order matters:

1. `generate_nonwater_contact_panel.py` fixes each solute and optimizes its neighbouring NH₃, CH₃OH or CH₃CN with MMFF94s. For each solute/species pair, it chooses the lower-energy converged pose from the two original starting arrangements. Both the source list and the two geometric acceptance rules are fixed before new QM work.
2. `predict_geometry.py` computes frozen GLIDER predictions, and `predict_mace_polar_cp_response.py` plus `materialize_polar_sites.py` compute the unfitted MACE-POLAR-1-L baseline on the same probe points.
3. `freeze_nonwater_contact_panel.py` checks exactly 36 geometry and prediction files, hashes them, and refuses to run after a reference directory exists. Do not overwrite its manifest.
4. `acquire_nonwater_contact_qm_cpu.py` computes the same counterpoise-consistent response target as the paper using ωB97X-D3(BJ)/def2-TZVPD. It checks the frozen geometry and prediction hashes before each calculation. The CPU calculation is expensive. For iodine-containing solute fragments, an isolated atomic-density solution is projected into the full ghost basis as an initial guess; the final SCF, energy and observables still use the unchanged counterpoise Hamiltonian. A same-Hamiltonian second-order stage remains available if DIIS fails. The official `dftd3` library provides the D3(BJ) correction when PySCF's optional dispersion package is unavailable; D3 does not change the SCF potential.
5. `score_nonwater_contact_panel.py` requires all 36 references, checks hashes and shared probe coordinates, and reports equal-solute errors and solute-blocked paired uncertainty without selecting cases by outcome.

Install the model dependencies as in the [inference guide](../benchmark/README.md), plus `pip install -e '.[contact]'`. The release contains the generated geometries, frozen predictions, QM references and score tables, so readers can **re-score** without rerunning the expensive SCF calculations. The generation workflow is for an independent replay.

The original panel and its contact audit remain under `experiments/nonwater/`. The completed follow-up is stored separately under `experiments/nonwater_contact/`. A negative QM interaction energy is a fixed-geometry contact check, not a decomposition of the response into polarization, exchange or charge transfer.

## Fragment-separation QM check

`run_cyclic_carbamate_separation_qm.py` calculates complex-minus-ghost-fragments references for the archived `dev_cyclic_carbamate` trajectories at ten separations, with one water or the intact four-water environment moved. `summarize_cyclic_carbamate_separation_qm.py` validates the 20 case records and regenerates the figure and result note. The complete [released package](../../experiments/dissociation_qm/) can be checked without rerunning SCF using `scripts/verify_companion.py`. See [the command recipe](../../docs/reproduce.md#direct-qm-check-of-one-separation-trajectory).
