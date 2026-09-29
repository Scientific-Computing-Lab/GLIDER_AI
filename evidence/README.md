# Integrity and provenance

The [top-level provenance file](../PROVENANCE.json) records the exact source commit and SHA-256s of the curated code, checkpoint and camera-ready numeric tables. [`frozen_source/`](frozen_source/README.md) contains the byte-preserved historical response implementation; [`manifests/`](manifests/README.md) records freezes of the prospective experiments. [`CLEAN_IMPLEMENTATION_EQUIVALENCE.json`](CLEAN_IMPLEMENTATION_EQUIVALENCE.json) compares the release-facing predictor against the historical outputs. [`FINAL_STATISTICS.json`](FINAL_STATISTICS.json) retains frozen panel summary statistics.

The archive distinguishes three stages: (1) sparse response-head development; (2) frozen model predictions on new chemistry; (3) QM references and scoring. That chronology supports a prospective *response-supervision* split. It does not establish what molecular identities a broad third-party foundation model encountered during its own pretraining.
