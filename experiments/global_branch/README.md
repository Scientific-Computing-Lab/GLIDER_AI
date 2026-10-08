# Does the global dipole branch help?

This controlled post hoc ablation removes the global dipole readout and moment reconciliation while keeping the local mapping and the same averaged frozen prior. Both variants use the original 48 configurations, five seeds, 64 epochs, optimizer and loss weights. Test labels are used only for final evaluation.

[Protocol](protocol.json) · [Matched full checkpoint](matched_full.pt) · [Ablated checkpoint](without_global.pt) · [Predictions](predictions/) · [Per-configuration results](results.csv) · [Panel summary](summary.csv) · [Paired evaluation](evaluation.json)

The matched full model reproduces the published head exactly on the common newly extracted features. Removing the global branch raises mean ESP NRMSE by 0.0114 across the 56 solutes, with a solute-bootstrap 95% interval [0.0008, 0.0203]. Fourteen solutes improve without the branch. This supports a modest average contribution, not uniform benefit.

The new common feature extraction produces small rounding differences from the old primary-panel tables. Those original tables and predictions are unchanged.

[Reproduction procedure](../../docs/reproduce.md#author-review-controls) · **Paper:** Appendix K.3, Table S14.
