# Lightweight GLIDER utilities

This importable Python package contains small pieces needed to **inspect and score** the response archive: registry validation ([`data/`](data/README.md)), potential/dipole reconstruction and molecular probe surfaces ([`inference/`](inference/README.md)), metric aggregation and paired bootstrap ([`metrics/`](metrics/README.md)), and a pointer to the frozen trained head ([`model/`](model/README.md)).

The source of the scientific response model itself is the byte-preserved [`evidence/frozen_source/prospective_1/response_learning.py`](../evidence/frozen_source/prospective_1/response_learning.py). Keeping that file distinct from convenience utilities prevents a later refactor from being mistaken for the code that generated the frozen predictions.
