# Panel-III scoring chronology

The first invocation of the frozen scoring program completed the identical
metric calculations for every method but stopped while serializing metadata:
the selected-identity table column is named `freesolv_id`, whereas the report
writer requested `source_id`. No model-specific path, prediction, reference,
metric, aggregation, rank or bootstrap logic was involved in the exception.

The objective one-token metadata correction (`source_id` to `freesolv_id`) was
made before any complete score report existed. The full program was then rerun
from all frozen prediction and reference files for every method. No prediction,
reference, benchmark identity, geometry, metric or aggregation rule changed.
