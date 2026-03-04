# scOPE — single-cell Oncological Prediction Explorer

[![PyPI version](https://badge.fury.io/py/scope-bio.svg)](https://badge.fury.io/py/scope-bio)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)

<picture>
  <source media="(prefers-color-scheme: dark)"
          srcset="https://raw.githubusercontent.com/Ashford-A/scOPE/main/assets/figures/scOPE_overview_dark.png">
  <img src="https://raw.githubusercontent.com/Ashford-A/scOPE/main/assets/figures/scOPE_overview_light.png"
       alt="scOPE transfer-learning method overview"
       width="100%">
</picture>

---

## scOPE workflow (transfer learning from bulk → single-cell)

**scOPE is a transfer-learning framework that uses *bulk RNA-seq cohorts with known mutation status* to learn a compact, biologically meaningful latent space, then projects *single-cell RNA-seq* into that same space to predict the probability that **specific cancer-associated gene mutations** are present in individual cells. This provides mutation-informed subclonal structure that can complement CNV methods and increase subclonal granularity.

### Overview

scOPE proceeds in two phases:

### 1) Learn latent factors from bulk RNA-seq and train mutation classifiers (Panels a–c)

- **i. Bulk expression matrix**  
  Construct a bulk cohort expression matrix **A_bulk** (rows = patient samples, columns = genes).  
  The bulk matrix is **normalized / centered / scaled** to ensure comparable gene-wise signal.

- **ii. Latent feature mapping via SVD**  
  Decompose the normalized bulk matrix using SVD:

  <p align="center">
    <picture>
      <source media="(prefers-color-scheme: dark)"
              srcset="https://latex.codecogs.com/svg.image?\color{white}{A_{\text{bulk}}=U_{\text{bulk}}\Sigma_{\text{bulk}}V^{\top}}">
      <img src="https://latex.codecogs.com/svg.image?A_{\text{bulk}}=U_{\text{bulk}}\Sigma_{\text{bulk}}V^{\top}"
           alt="A_bulk = U_bulk Σ_bulk V^T">
    </picture>
  </p>

  - **U_bulk**: sample scores (rows = patients, columns = latent factors)  
  - **Σ_bulk**: diagonal matrix of singular values  
  - **V**: gene loadings (rows = genes, columns = latent factors)

  Define the bulk latent representation (patient-by-factor embedding):

  <p align="center">
    <picture>
      <source media="(prefers-color-scheme: dark)"
              srcset="https://latex.codecogs.com/svg.image?\color{white}{Z_{\text{bulk}}=U_{\text{bulk}}\Sigma_{\text{bulk}}}">
      <img src="https://latex.codecogs.com/svg.image?Z_{\text{bulk}}=U_{\text{bulk}}\Sigma_{\text{bulk}}"
           alt="Z_bulk = U_bulk Σ_bulk">
    </picture>
  </p>

- **iii. Train mutation-prediction models in latent space**  
  For each mutation / gene-of-interest, train a supervised ML model to predict mutation presence **Y** from **Z_bulk**.  
  This yields one (or multiple) mutation-specific classifiers operating on the learned latent factors.

---

### 2) Project scRNA-seq into the bulk-derived latent space and predict mutations per cell (Panels d–f)

- **i. Single-cell expression matrix**  
  Construct a single-cell expression matrix **A_sc** (rows = single cells, columns = genes).

- **ii. Normalize scRNA-seq using bulk-derived parameters**  
  Apply the *same* gene-wise normalization / centering / scaling learned from the bulk cohort to obtain **A′_sc**.  
  (This alignment step makes the projection comparable across bulk and single-cell.)

- **iii. Project cells into latent space and infer mutation probabilities**  
  Use the bulk-derived gene loadings **V** to compute the single-cell latent representation:

  <p align="center">
    <picture>
      <source media="(prefers-color-scheme: dark)"
              srcset="https://latex.codecogs.com/svg.image?\color{white}{Z_{\text{sc}}=A'_{\text{sc}}V}">
      <img src="https://latex.codecogs.com/svg.image?Z_{\text{sc}}=A'_{\text{sc}}V"
           alt="Z_sc = A'_sc V">
    </picture>
  </p>

  Then apply the trained bulk models to **Z_sc** to predict **per-cell mutation probabilities**, producing mutation-informed cellular maps that can be analyzed alongside expression programs, clusters, and CNV signals.

---

## Installation

### From PyPI
```bash
pip install scope-bio
```

### With optional dependencies (UMAP, XGBoost, LightGBM, SHAP)
```bash
pip install scope-bio[full]
```

### From conda-forge
```bash
conda install -c conda-forge scope-bio
```

### Development install
```bash
git clone https://github.com/Ashford-A/scOPE.git
cd scOPE
conda env create -f environments/scope-dev.yml
conda activate scope-dev
pip install -e ".[dev]"
```

---

## Quick start

```python
import anndata as ad
import pandas as pd
from scope import BulkPipeline, SingleCellPipeline
from scope.io import load_mutation_labels

# --- Phase 1: Bulk --------------------------------------------------------
adata_bulk = ad.read_h5ad("bulk_cohort.h5ad")
mutation_labels = load_mutation_labels("mutations.csv", sample_col="sample_id")

bulk_pipe = BulkPipeline(
    norm_method="cpm",
    decomposition="svd",   # "svd" | "nmf" | "ica" | "pca"
    n_components=50,
    classifier="logistic", # "logistic" | "random_forest" | "gbm" | "xgboost" | "lightgbm" | "svm" | "mlp"
)
bulk_pipe.fit(adata_bulk, mutation_labels, cv=5)
bulk_pipe.save("models/bulk_pipeline.pkl")

# --- Phase 2: Single cell --------------------------------------------------
adata_sc = ad.read_h5ad("sc_tumor.h5ad")

# Prepare a preprocessed bulk reference for moment matching
adata_bulk_pp = bulk_pipe.preprocessor_.transform(adata_bulk)

sc_pipe = SingleCellPipeline(
    bulk_pipeline=bulk_pipe,
    alignment_method="z_score_bulk",  # "z_score_bulk" | "moment_matching" | "quantile" | "none"
)
sc_pipe.fit(adata_bulk_pp, adata_sc)
adata_sc = sc_pipe.transform(adata_sc)

# adata_sc.obs now contains columns: mutation_prob_KRAS, mutation_prob_TP53, ...

# --- Visualise -------------------------------------------------------------
from scope.visualization import compute_umap, plot_mutation_probabilities

adata_sc = compute_umap(adata_sc, obsm_key="X_svd")
fig = plot_mutation_probabilities(adata_sc, mutations=["KRAS", "TP53"])
fig.savefig("mutation_probs.pdf", bbox_inches="tight")
```

---

## API reference

### Preprocessing
| Class | Description |
|---|---|
| `BulkNormalizer` | CPM / TPM / median-ratio / TMM normalisation |
| `BulkScaler` | Gene-wise centering and scaling |
| `BulkPreprocessor` | Combined normalise + scale |
| `SingleCellPreprocessor` | QC filter + normalise + optional scale |
| `BulkSCAligner` | z-score / moment-matching / quantile alignment |

### Decomposition
| Class | Description |
|---|---|
| `SVDDecomposition` | Truncated SVD (randomized / ARPACK / full) |
| `NMFDecomposition` | Non-negative matrix factorization |
| `ICADecomposition` | FastICA |
| `PCADecomposition` | PCA via sklearn |
| `get_decomposition(name)` | Factory function |

### Classification
| Class | Description |
|---|---|
| `LogisticMutationClassifier` | L1/L2/ElasticNet logistic regression |
| `RandomForestMutationClassifier` | Random forest |
| `GBMMutationClassifier` | Gradient boosting (sklearn) |
| `XGBMutationClassifier` | XGBoost |
| `LGBMMutationClassifier` | LightGBM |
| `SVMMutationClassifier` | SVM + Platt calibration |
| `MLPMutationClassifier` | Multi-layer perceptron |
| `PerMutationClassifierSet` | Trains/stores one classifier per mutation |
| `get_classifier(name)` | Factory function |

### Evaluation
| Function | Description |
|---|---|
| `evaluate_classifier` | AUROC, AUPRC, Brier score |
| `evaluate_all` | Evaluate all mutations at once |
| `cross_validate_classifiers` | Stratified k-fold CV |
| `roc_curve_data` / `pr_curve_data` | Curve arrays for plotting |

### Visualization
| Function | Description |
|---|---|
| `compute_umap` | UMAP on latent embedding |
| `compute_tsne` | t-SNE on latent embedding |
| `plot_embedding` | Scatter by categorical or continuous |
| `plot_mutation_probabilities` | Grid of per-mutation probability overlays |
| `plot_scree` | Singular value / EVR scree plot |
| `plot_mutation_heatmap` | Mean probability per cluster heatmap |

---

## Running tests

```bash
pytest tests/ -v --cov=scope --cov-report=term-missing
```

---

## Citation

If you use scOPE in your research, please cite:

> Ashford, A. et al. (2024). scOPE: transfer-learning from bulk RNA-seq to infer
> per-cell mutation probabilities in single-cell transcriptomics.
> *[Journal]* doi: ...

---

## License

MIT — see [LICENSE](LICENSE).
