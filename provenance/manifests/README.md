# Prospective freeze manifests

These JSON files retain content hashes and timing boundaries for Panel I, Panel II and comparator outputs. [`prospective_freeze_manifest.json`](prospective_freeze_manifest.json) and [`NEXTGEN_FREEZE_MANIFEST.json`](NEXTGEN_FREEZE_MANIFEST.json) are the first two panel anchors; the remaining files document prediction and reference freezes, comparator coverage and related integrity checks.

The [release verifier](../../scripts/reproduce/verify_release.py) checks the archived checkpoint, registries and predicted physical constraints. Hashes and `reference_response_labels_exist` fields support the chronology, but do not imply that third-party pretrained models lacked access to a chemical identity in their much larger pretraining corpora.
