# Held-out-water plotting values

This archive backs camera-ready Figure 4 and the related supplementary coupling analysis. The base is always a solute plus W1–W3. W4 (or a more distant water) is excluded from the base calculation and from the GLIDER input, then coupled **one-way** to the frozen response.

| Start here | Meaning |
|:--|:--|
| [`PLOTTING_MANIFEST.json`](PLOTTING_MANIFEST.json) | SHA-256s, ten-solute panel size, 90 range probes and panel-wide origins for exact-dipole controls. |
| [`energy_pairs.csv`](energy_pairs.csv), [`torque_pairs.csv`](torque_pairs.csv), [`orientation_pairs.csv`](orientation_pairs.csv) | Paired W4 comparisons across methods and solutes. |
| [`energy_rank.csv`](energy_rank.csv), [`torque_rank.csv`](torque_rank.csv) | Independent external-water ranks W4–W12. |
| [`multipole_*.csv`](multipole_energy.csv) | Exact response-moment hierarchy. |
| [`distance_errors.csv`](distance_errors.csv), [`fractions.csv`](fractions.csv) | Distance and response-fraction diagnostics. |

The numerically distinct energy MAE values used in the root graphic are W4: GLIDER 0.0487, polar baseline 0.114, exact response dipole 0.149 and zero response 0.0853 kcal mol⁻¹. The exact dipole's origin is selected once for the whole endpoint panel, not separately for each solute. These are not total interaction energies.
