# Registry validation

[`registry.py`](registry.py) validates the expected counts and uniqueness of frozen configuration registries. The [release verifier](../../scripts/reproduce/verify_release.py) uses it for the 48, 96 and 80 prospective configurations. It checks record identity, not model quality by itself.
