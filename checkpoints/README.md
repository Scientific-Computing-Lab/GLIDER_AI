# Frozen GLIDER response head

[`glider_site_response_ensemble.pt`](glider_site_response_ensemble.pt) is the **five-member response-head ensemble** selected from the original 14-chemistry, 48-environment development set. Its SHA-256 is `a288afa285128e13cb7a05ba9459f1dce07e79b8638c6e272c2509ae68d7a288`; [`glider_site_response_ensemble.manifest.json`](glider_site_response_ensemble.manifest.json) records the seeds, selection inputs and the absence of energy/force supervision. The [release verifier](../scripts/reproduce/verify_release.py) checks the hash and frozen predictions.

This file is **not** a stand-alone geometry-to-field model. It consumes features from the frozen MACE-POLAR-1-M encoder and a starting response averaged from frozen MACE-POLAR-1-M/-1-L predictions. The [geometry inference command](../scripts/benchmark/README.md) assembles those pieces. Official MACE weights remain with their publisher and can be fetched with the verified [download script](../scripts/reproduce/fetch_mace_polar.py).

The checkpoint is a PyTorch archive. Load it only from a source you trust; the inference path in this repository uses `torch.load(..., weights_only=False)` to restore the historical trained ensemble format.
