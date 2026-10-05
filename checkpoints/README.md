# The frozen five-member response head

[glider_site_response_ensemble.pt](glider_site_response_ensemble.pt) is the original model used throughout the paper. The [manifest](glider_site_response_ensemble.manifest.json) records its training inputs. Its content hash is checked by the release verifier.

It requires frozen MACE-POLAR-1-M features and the averaged M/L response prior. [Complete inference command](../scripts/benchmark/README.md). The new control checkpoints are kept separately under [global_branch](../experiments/global_branch/).

[Training data](../experiments/training/) · [Model definition](../docs/model.md) · [Known dissociation failure](../experiments/dissociation/)
