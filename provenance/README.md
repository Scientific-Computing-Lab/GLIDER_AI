# Historical records and calculation sources

This archive supports the experiment folders. Start with [data lineage](../docs/data_lineage.md) for a readable account of the original training set, Panel-I reuse and later controls.

- [Frozen implementation](frozen_source/) preserves the scientific definitions used by the released checkpoint.
- [Sampling sources](sampling/) document cluster construction and liquid-snapshot selection.
- [Downstream calculation sources](calculation_source/) retain the QM-density and coupling implementations.
- [Panel-III records](panel_3/) preserve its identity selection, source attribution and freeze audits.
- [Later development pool](later_development_pool/) is deliberately separate from the original 48 training examples.
- [Retained failures](negative_results/) record the Panel-I joint gate and later unsuccessful model.
- [Migration map](release_path_map.json) links old archive paths to release locations.
- [Previous table archive](previous_release_tables/) retains unique old derivation tables. It is historical material; use the current [paper map](../docs/paper_map.md) and [figure folders](../figures/) for the final numbering.

Original machine paths in archived manifests document past runs. They are not active download or execution dependencies. Current runnable inputs are repository-relative and checked by the release verifier. Some archived scripts require their original project layout; supported entry points are under `scripts/`.
