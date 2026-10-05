# From frozen polar features to response sites

GLIDER adapts information from MACE-POLAR, whose pretraining includes energies and forces. Its response head learns the interaction-induced potential and dipole directly.

1. Extract geometry-dependent MACE-POLAR-1-M features.
2. Predict complex and isolated-fragment sites with frozen M and L checkpoints. Subtract fragments and average the two resulting responses.
3. Learn local corrections to the atom-centred charges and dipoles, plus a global response-dipole correction.
4. Remove net response charge and distribute the mismatch between site moment and global prediction over atomic dipoles.
5. Evaluate the resulting potential at any requested off-site points.

The global constraint is exact for the predicted moment. It does not make either prediction exact relative to QM, and it does not enforce dissociation.

## Training

The original 48 configurations are listed in [training](../experiments/training/). The final ensemble uses seeds 2026081300–2026081304, 64 epochs, AdamW at 0.002 with weight decay 0.00002, cosine scheduling and gradient clipping at 5. The objective combines response ESP with dipole weight 2 and auxiliary fitted-site weight 0.05. Architecture and transformations are retained in the [frozen implementation](../provenance/frozen_source/prospective_1/response_learning.py).

The response head was selected through nested leave-one-chemistry-out development. All 14 chemistries enter final fitting. Panels I–III never enter that final fitting operation. Panel I later trained a separate unsuccessful model family, as explained in [data lineage](data_lineage.md).

## Inference and physical limits

The [geometry wrapper](../scripts/benchmark/README.md) assembles the official frozen M/L models and the released five-member head. The model supports the compact neutral systems evaluated here. The [dissociation test](../experiments/dissociation/) demonstrates a residual distant response. The current checkpoint should not be treated as a general molecular-dynamics force field or as a validated fragment-separable response model.
