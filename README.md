# scOPE — single-cell Oncological Prediction Explorer

[![PyPI version](https://img.shields.io/pypi/v/scope-bio?v=0.1.8)](https://pypi.org/project/scope-bio/)
[![pypi downloads](https://img.shields.io/pepy/dt/scope-bio?label=pypi%20downloads)](https://pepy.tech/project/scope-bio)
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

scOPE is a transfer-learning framework that uses *bulk RNA-seq cohorts with known mutation status* to learn a compact, biologically meaningful latent space, then projects *single-cell RNA-seq* into that same space to predict the probability that **specific cancer-associated gene mutations** are present in individual cells. This provides mutation-informed subclonal structure that can complement CNV methods and increase subclonal granularity.

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
    decomposition="svd",   # "svd" | "nmf" | "ica" | "pca" | "fa" | "cnmf"
    n_components=50,
    classifier="logistic", # "logistic" | "random_forest" | "gbm" | "xgboost" | "lightgbm" | "svm" | "mlp"
)
bulk_pipe.fit(adata_bulk, mutation_labels, cv=5)
bulk_pipe.save("models/bulk_pipeline.pkl")

# --- Phase 2: Single cell --------------------------------------------------
adata_sc = ad.read_h5ad("sc_tumor.h5ad")

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

## Preprocessing options

`BulkPreprocessor` and `SingleCellPreprocessor` are designed to handle the full range of input states — from raw counts to already-normalized matrices — without requiring you to re-implement preprocessing outside the pipeline.

### Bulk: handling already-processed inputs

```python
bulk_pipe = BulkPipeline(
    norm_method="none",               # skip library-size normalization
    log1p=False,                      # data is already log-transformed
    decomposition="svd",
    n_components=50,
)

# Or pass flags directly to BulkPreprocessor:
from scope.preprocessing import BulkPreprocessor

preprocessor = BulkPreprocessor(
    norm_method="cpm",
    log1p=True,
    already_log_transformed=False,    # set True if input is e.g. log2-TPM from GEO
    min_samples_expressed=5,          # remove genes expressed in fewer than 5 samples
    min_expression=0.5,               # expression threshold for the above filter
    gene_blacklist=["MALAT1", "NEAT1"],
    auto_remove_mito=True,            # remove MT- genes
    auto_remove_ribo=True,            # remove RPS/RPL genes
    run_hvg=True,                     # select highly variable genes (requires scanpy)
    n_hvg=3000,
    hvg_flavor="seurat_v3",
    batch_key="cohort",               # batch correction via obs column
    batch_method="combat",            # "combat" | "harmony"
)
```

### Single-cell: QC, mito filtering, and doublet removal

```python
from scope.preprocessing import SingleCellPreprocessor

sc_prep = SingleCellPreprocessor(
    filter_strategy="both",           # "min_counts" | "min_genes" | "both" | "none"
    min_counts=500,
    min_genes=300,
    max_counts=25000,
    max_genes=6000,                   # upper bound (doublet proxy)
    max_mito_pct=20.0,                # remove cells with >20% mitochondrial reads
    auto_flag_mito=True,              # annotate pct_mito in adata.obs regardless
    run_doublet_detection=True,       # Scrublet-based doublet removal (pip install scrublet)
    doublet_threshold=None,           # None = automatic Scrublet threshold
    already_qc_filtered=False,        # True = skip all cell-level filters
    already_normalized=False,         # True = skip library-size normalization
    already_log_transformed=False,    # True = skip log1p
)
```

---

## Decomposition methods

scOPE supports several latent-space methods, all sharing the same `fit` / `transform` / `components_` interface and usable as a drop-in via the `decomposition=` argument in `BulkPipeline`.

| Key | Class | Notes |
|---|---|---|
| `"svd"` | `SVDDecomposition` | Default. Linear, interpretable. Gene loadings V enable direct sc projection. |
| `"nmf"` | `NMFDecomposition` | Non-negative, additive gene programs (metagenes). Requires non-negative input. |
| `"ica"` | `ICADecomposition` | Independent components. Useful for finding non-Gaussian expression sources. |
| `"pca"` | `PCADecomposition` | Standard PCA. Equivalent to SVD on centred data. |
| `"fa"` | `FactorAnalysisDecomposition` | Probabilistic FA. Accounts for gene-specific noise variance (heteroscedasticity). |
| `"cnmf"` | `ConsensusNMFDecomposition` | Consensus NMF (Kotliar et al., eLife 2019). Runs NMF *n* times, clusters components for stability. Recommended over single-run NMF for publication. |

```python
# Consensus NMF example
bulk_pipe = BulkPipeline(
    decomposition="cnmf",
    decomposition_kwargs={"n_iter": 50, "n_components": 20},
    classifier="logistic",
)

# Factor Analysis example
bulk_pipe = BulkPipeline(
    decomposition="fa",
    n_components=30,
    classifier="logistic",
)
```

---

## SVD evaluation

When using `decomposition="svd"`, `SVDEvaluator` produces a comprehensive suite of plots and a gene program table that characterize which components drive classification and which genes define them. This is useful for reviewers, manuscript figures, and biological interpretation.

```python
from scope.evaluation import SVDEvaluator

# After fitting:
adata_pp = bulk_pipe.transform_bulk(adata_bulk)
Z_bulk = adata_pp.obsm[bulk_pipe.obsm_key_]

ev = SVDEvaluator(bulk_pipe, Z_bulk, mutation="KRAS")
ev.run_all(output_dir="figures/svd_eval_KRAS")
```

`run_all()` saves the following to `output_dir/`:

| Output | Description |
|---|---|
| `weighted_scree.png` | Scree plot with bars colour-coded by `\|coef\|×σ` classifier importance |
| `component_importance.png` | Ranked bar chart of importance; decomposed `\|coef\|` vs σ |
| `gene_loading_heatmap.png` | Hierarchically clustered heatmap of top gene loadings × top components |
| `top_genes_per_component.png` | Signed gene bar charts for each top-important component |
| `latent_scatter.png` | Pairwise scatter of top components coloured by mutation label |
| `separation_violins.png` | Component score distributions by mutation status + Mann-Whitney p |
| `component_label_correlation.png` | Spearman ρ heatmap between components and mutation label (FDR-corrected) |
| `roc_ablation.png` | Cross-validated ROC curve + AUC vs. n_components retained |
| `permutation_test.png` | Observed AUROC vs. permutation null distribution |
| `gene_biplot.png` | Sample scores + gene loading arrows for top-2 components |
| `shap_summary_dot.png` | SHAP dot summary (requires `shap`) |
| `shap_summary_bar.png` | SHAP bar summary (requires `shap`) |
| `umap_zbulk.png` | UMAP of Z_bulk coloured by mutation label (requires `umap-learn`) |
| `calibration_curve.png` | Reliability diagram for predicted mutation probabilities |
| `component_crosscorr.png` | Pearson correlation among SVD components (flags batch bleed) |
| `gene_program_table.csv` | Tidy table: component rank, σ, `\|coef\|`, importance, top genes + loadings |

Individual plots can also be called directly:

```python
ev.plot_weighted_scree(output_dir=Path("figures/"))
ev.plot_separation_violins(output_dir=Path("figures/"), top_components=12)
ev.export_gene_program_table(output_dir=Path("figures/"), top_components=10)
```

---

## API reference

### Preprocessing
| Class | Description |
|---|---|
| `BulkNormalizer` | CPM / TPM / median-ratio / TMM normalisation |
| `BulkScaler` | Gene-wise centering and scaling |
| `BulkPreprocessor` | Combined normalise + scale, with gene filtering, HVG, and batch correction |
| `SingleCellPreprocessor` | QC filter + mito filter + doublet removal + normalise + optional scale |
| `BulkSCAligner` | z-score / moment-matching / quantile alignment |

### Decomposition
| Class | Description |
|---|---|
| `SVDDecomposition` | Truncated SVD (randomized / ARPACK / full) |
| `NMFDecomposition` | Non-negative matrix factorization |
| `ICADecomposition` | FastICA |
| `PCADecomposition` | PCA via sklearn |
| `FactorAnalysisDecomposition` | Probabilistic factor analysis (heteroscedastic noise) |
| `ConsensusNMFDecomposition` | Consensus NMF for stable gene program discovery |
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
| Class / Function | Description |
|---|---|
| `SVDEvaluator` | Full SVD component interpretation suite (15 plots + gene program CSV) |
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

## Optional dependencies

| Package | Purpose | Install |
|---|---|---|
| `umap-learn` | UMAP embeddings and SVDEvaluator UMAP plot | `pip install scope-bio[full]` |
| `shap` | SHAP component importance in SVDEvaluator | `pip install scope-bio[full]` |
| `xgboost` | XGBoost classifier | `pip install scope-bio[full]` |
| `lightgbm` | LightGBM classifier | `pip install scope-bio[full]` |
| `scrublet` | Doublet detection in SingleCellPreprocessor | `pip install scrublet` |
| `combat` | ComBat batch correction in BulkPreprocessor | `pip install combat` |
| `harmonypy` | Harmony batch correction in BulkPreprocessor | `pip install harmonypy` |
| `statsmodels` | FDR correction in SVDEvaluator correlation heatmap | `pip install statsmodels` |

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

(Work currently in progress - will be filled out later once preprint is up on bioRxiv)

---

## License

MIT — see [LICENSE](LICENSE).
