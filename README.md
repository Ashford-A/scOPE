# scOPE — single-cell Oncological Prediction Explorer

[![PyPI version](https://badge.fury.io/py/scope-bio.svg)](https://badge.fury.io/py/scope-bio)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)

**scOPE** is a transfer-learning framework that learns mutation-linked
expression programs from bulk RNA-seq cohorts and projects them into
single-cell RNA-seq data to infer per-cell mutation probabilities —
revealing subclonal evolutionary structure that complements CNV analyses.

---

## Workflow overview

```
Bulk RNA-seq  ──► Normalise ──► SVD / NMF / ICA / PCA ──► Train mutation classifiers
                                        │
                                        ▼  gene loadings V
scRNA-seq     ──► Normalise ──► Moment-match to bulk ──► Project (Z_sc = A'·V) ──► Predict P(mut | cell)
```

**Phase 1** (bulk):
1. Normalise bulk counts (CPM / TPM / median-ratio / TMM).
2. Centre and scale gene-wise.
3. Decompose with SVD (default), NMF, ICA, or PCA.
4. Train per-mutation classifiers (logistic, random forest, GBM, XGBoost, LightGBM, SVM, MLP) on latent coordinates.

**Phase 2** (single cell):
1. QC-filter and normalise sc counts.
2. Align sc distribution to bulk (z-score / moment-matching / quantile).
3. Project cells into the bulk latent space using stored gene loadings V.
4. Apply trained classifiers → per-cell mutation probability columns in `adata.obs`.

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
