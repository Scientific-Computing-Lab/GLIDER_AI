# Molecule-aware metrics

[`core.py`](core.py) contains equal-solute aggregation, paired molecule-blocked bootstrap and leave-one-solute-out influence utilities. The [statistics script](../../scripts/reproduce/recompute_statistics.py) reads archived per-solute results and uses these functions to recompute panel summaries. A configuration is not an independent bootstrap block when several configurations share one solute.
