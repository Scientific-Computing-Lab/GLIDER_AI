# GLIDER model card

## Intended use

Prediction of counterpoise-consistent mutual response ESP and induced dipole for complexes containing a neutral organic solute and structured multi-water clusters.

## Inputs and outputs

Input is atomic identity and geometry. The frozen broad MACE-POLAR features and unfitted M/L response prior feed a site-response head. Outputs are net-neutral induced site charges and equivariant site dipoles, from which response ESP and molecular dipole are reconstructed.

## Training

Fourteen development chemistries, 48 environments, five seeds. Direct response ESP and dipole supervision; no energy, force, hydration or prospective labels. Chemistry-blocked model selection. Exact settings are in the checkpoint manifest and frozen source.

## Validation

Prediction-before-label panels of 12, 24 and 20 chemistries absent from GLIDER response supervision. Panel III used a deterministic label-blind selection from pinned FreeSolv v0.52 after excluding historical exact matches and Morgan similarity at least 0.70. FreeSolv hydration values were not accessed.

## Limitations

Neutral organics, three-water environments and four designed regimes only. Not a conservative potential, not a force field, not validated for MD or free energy. MACE-POLAR checkpoint licenses apply to the prerequisite encoders. The broad encoder may have seen exact identities or analogues during foundation pretraining; the transfer claim concerns absence from GLIDER response supervision, not pretraining-naive chemistry.
