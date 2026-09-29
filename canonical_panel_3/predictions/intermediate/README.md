# Frozen Panel-III intermediate arrays

These feature/prior caches are the inputs used to apply the response head on Panel III. They come from frozen MACE-POLAR checkpoints and are not QM response labels. Keeping them makes the prediction step auditable without redownloading foundation weights. The [prediction freeze manifest](../ALL_PREDICTIONS_FREEZE_MANIFEST.json) records the model boundary.
