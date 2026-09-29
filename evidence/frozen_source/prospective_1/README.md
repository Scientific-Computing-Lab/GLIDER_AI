# Frozen Panel-I implementation

`response_learning.py` defines the response head, site-basis design and original fitting/inference logic. `extract_mace_features.py` and `predict_mace_polar_cp_response.py` call official MACE-POLAR checkpoints to create geometry features and complex-minus-fragments site predictions; `average_qu_predictions.py` makes the M/L starting response. These files were copied **byte-for-byte** from source commit `a774f6a`; their hashes are in [`PROVENANCE.json`](../../../PROVENANCE.json).

The model head is trained from 48 QM-labelled environments; the MACE encoders remain frozen. The [checkpoint manifest](../../../checkpoints/glider_site_response_ensemble.manifest.json) records training seeds and the response-only supervision boundary. The [clean inference wrapper](../../../scripts/benchmark/README.md) is easier to use on new geometries but still imports this frozen `ResponseHead`.
