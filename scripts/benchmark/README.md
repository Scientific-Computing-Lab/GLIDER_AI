# Geometry-to-response inference

[`predict_geometry.py`](predict_geometry.py) is the entry point. It extracts MACE-M features, obtains MACE-M and MACE-L complex-minus-fragments predictions, averages the starting response, then applies the frozen five-member GLIDER head. [`predict_from_features.py`](predict_from_features.py) performs the final cached-feature stage using the historical `ResponseHead` definition.

```bash
python scripts/benchmark/predict_geometry.py \
  --configurations examples/one_response_geometry.extxyz \
  --mace-polar-m third_party/checkpoints/MACE-POLAR-1-M.model \
  --mace-polar-l third_party/checkpoints/MACE-POLAR-1-L.model \
  --output example_prediction --device cpu
```

Install this package with `pip install -e '.[model]'` and the matching [MACE source](https://github.com/ACEsuit/mace) with `pip install 'git+https://github.com/ACEsuit/mace.git@91df5a2032b24ff9e23e0dc7b9407dde0da6fb31'`. Fetch the official foundation weights with [`fetch_mace_polar.py`](../reproduce/fetch_mace_polar.py). CPU is supported but may be slow; CUDA is the script's default device. The output path must not already exist. Each output NPZ contains surface points, response ESP, response dipole and atom-centred charges/dipoles, indexed by `prediction_registry.csv`. Predictions are a frozen electronic-response observable, not a total force field.
