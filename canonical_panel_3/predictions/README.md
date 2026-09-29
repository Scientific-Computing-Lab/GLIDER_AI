# Panel-III predictions, frozen before QM

[`ALL_PREDICTIONS_FREEZE_MANIFEST.json`](ALL_PREDICTIONS_FREEZE_MANIFEST.json) records the prediction-time lock and checkpoint. [`glider/`](glider/README.md) holds the fixed GLIDER outputs; [`comparators/`](comparators/README.md) holds eligible public baselines; [`intermediate/`](intermediate/README.md) retains frozen MACE feature/prior caches. The later QM targets are isolated in [`references/`](../references/README.md). Align arrays by the parent [configuration registry](../data/configuration_registry.csv), not by filename order.
