# Panel I · the retained preregistered gate result

Panel I improved the spatial response-ESP endpoint, but it **did not pass** a separate, broader preregistered joint gate requiring both ESP and dipole gains. The 0.3218 ESP NRMSE exceeded that gate's 0.30 cutoff, and the dipole improvement over its gate comparator was too small and statistically unresolved. The frozen model, all 48 environments and this failed decision were retained.

- [`prospective_1_gate_result.json`](prospective_1_gate_result.json) is the machine-readable gate outcome.
- [`PROSPECTIVE_1_PREREGISTERED_GATE_FAILURE.md`](PROSPECTIVE_1_PREREGISTERED_GATE_FAILURE.md) is the original detailed report, copied without editing.
- The exact [frozen scoring script](../../evidence/frozen_source/prospective_1/score_prospective_response.py) and [freeze manifest](../../evidence/manifests/prospective_freeze_manifest.json) preserve the original protocol.

The gate used its original comparator set and conditions. The paper's later response-ESP comparison against the unfitted MACE-POLAR-1-L **polar baseline** is a different, explicitly reported analysis. [Return to the paper-to-data map](../README.md).
