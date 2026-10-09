# Fragment separation at 3–20 Å

For each of the 14 original training solutes, construct either a solute plus its first water **with the other waters omitted**, or a solute plus its intact original environment. Translate the water or full environment away from the stationary solute. Keep internal fragment geometry and solute-centred probe points fixed. Recompute the polar features, averaged prior and frozen GLIDER head at every geometry.

[Geometries](configurations.extxyz) · [Distance registry](configuration_registry.csv) · [Protocol](protocol.json) · [Probe points](solute_probe_points.npz) · [Predicted fields and sites](predicted_probe_potentials.npz) · [Per-case results](results.csv) · [Mean curves](summary.csv)

![Fragment-separation diagnostic](../../assets/separation-control.svg)

**The response does not vanish.** At 20 Å, the GLIDER mean surface-potential RMS is 0.564 mEh/e for a single water and 0.925 mEh/e for the intact environment. The averaged prior is lower at 0.195 and 0.301 mEh/e, respectively. The plot above also includes [the 50–100 Å extension](../dissociation_extended/), which shows that GLIDER's amplitude subsequently grows while the averaged prior declines. This folder contains the initial 3–20 Å data. These are predicted response amplitudes, not finite-distance errors against newly computed QM.

Together, the two prediction sweeps identify a failure of fragment-separation consistency. A [later one-solute QM check](../dissociation_qm/) confirms the far-separation error directly. Charge conservation and moment reconciliation do not impose that condition. **Paper:** workshop Fig. S7, journal Fig. 8c,d and Supplementary Section 11.4. [Reproduction](../../docs/reproduce.md#author-review-controls).

The intact-environment mean at 5 Å is dominated by `dev_thio_amide`, whose response RMS reaches about 1.91 Eh/e in both the averaged prior and GLIDER. This one case raises the 14-solute mean to about 138 mEh/e. The raw isolated-environment MACE-L output is already anomalous; its numerical cause remains unresolved. The case is retained in both plots and the [per-case data](results.csv). [Probe-to-environment clearance](probe_clearance.csv) and raw site values are included for diagnosis.
