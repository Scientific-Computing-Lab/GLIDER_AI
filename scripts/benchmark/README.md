# Predict a response from geometry

Install the model dependencies and paper-matched MACE implementation:

```bash
python -m pip install -e '.[model]'
python -m pip install 'git+https://github.com/ACEsuit/mace.git@91df5a2032b24ff9e23e0dc7b9407dde0da6fb31'
python scripts/reproduce/fetch_mace_polar.py
python scripts/benchmark/predict_geometry.py \
  --configurations examples/one_response_geometry.extxyz \
  --mace-polar-m third_party/checkpoints/MACE-POLAR-1-M.model \
  --mace-polar-l third_party/checkpoints/MACE-POLAR-1-L.model \
  --output example_prediction --device cpu
```

The wrapper extracts M features, predicts M/L complex-minus-fragment sites, averages the prior and runs the five-member GLIDER head. CUDA is also supported. The output path must not already exist.

[Cached-feature entry point](predict_from_features.py) · [Example geometry](../../examples/one_response_geometry.extxyz) · [Output arrays and units](../../docs/data_format.md)

**Model boundary:** this is a frozen electronic-response predictor for the tested compact neutral clusters. The current checkpoint fails the distant-fragment limit described in the [separation diagnostic](../../experiments/dissociation/).
