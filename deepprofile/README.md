# Deep neural profiling of head and neck cancer

VAE training, Integrated Gradients attribution, enrichment, classification and survival analyses used in
`zenodo_manuscript/`.

## Attribution

The modelling code in this folder is adapted from **DeepProfile**:

> Qiu W, Dincer AB, Janizek JD, Celik S, Pittet MJ, Naxerova K, Lee SI. Deep profiling of gene expression across
> 18 human cancers. *Nature Biomedical Engineering* 2025;9(3):333–355. https://doi.org/10.1038/s41551-024-01290-8

| File | Origin |
|---|---|
| `Create_PCs_for_DeepLearning_Models.py` | Adapted from DeepProfile |
| `MODEL_TRAININGS/VAE_3Layers_Model.py`, `Run_VAE_Models.py` | Adapted from DeepProfile (updated to TensorFlow 2 / Keras 3) |
| `MODEL_TRAININGS/IntegratedGradients.py`, `Get_VAE_IG_Attributions.py` | Adapted from DeepProfile |
| `MODEL_TRAININGS/Create_Ensemble_Labels.py`, `Create_DeepProfile_Training_Embeddings.py`, `Create_DeepProfile_Ensemble_Weights.py`, `Example_Run_All.py`, `gmeans.py` | Adapted from DeepProfile |
| `MODEL_TRAININGS/Select_Latent_Dimension_with_Gmeans.ipynb` | DeepProfile notebook; the saved outputs are from the original DeepProfile run across 18 cancer types, not from this project |
| `RESULT/`, `Single_Gene_Model/`, `SURVIVAL_ANALYSIS/` | Written for this project |

The head and neck cancer expression compendium (643 samples, 26 GEO datasets, 11,020 genes) is the one released with
DeepProfile.

## Configuration used in this project

- Input: 500 principal components of the compendium
- VAE: encoder 250 → 100 units (batch norm + ReLU), 50 latent variables, mirrored decoder; Adam (learning rate
  0.0005), batch size 50, 50 epochs; two runs with different seeds
- Integrated Gradients: zero baseline, 50 steps, projected to genes through the PCA loadings

DeepProfile itself trains 100 models at each of six latent sizes and combines them into a 150-node ensemble; this
project used the reduced configuration above because of computational limits.
