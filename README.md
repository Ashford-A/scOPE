# scOPE

<picture>
  <source media="(prefers-color-scheme: dark)"
          srcset="https://raw.githubusercontent.com/Ashford-A/scOPE/main/assets/figures/scOPE_overview_dark.png">
  <img src="https://raw.githubusercontent.com/Ashford-A/scOPE/main/assets/figures/scOPE_overview_light.png"
       alt="scOPE transfer-learning method overview"
       width="100%">
</picture>

---

## scOPE workflow (transfer learning from bulk → single-cell)

**scOPE (single-cell Oncological Prediction Explorer)** is a transfer-learning framework that uses *bulk RNA-seq cohorts with known mutation status* to learn a compact, biologically meaningful latent space, then projects *single-cell RNA-seq* into that same space to predict the probability that **specific cancer-associated gene mutations** are present in individual cells. This provides mutation-informed subclonal structure that can complement CNV methods and increase subclonal granularity.

### Overview

scOPE proceeds in two phases:

### 1) Learn latent factors from bulk RNA-seq and train mutation classifiers (Panels a–c)

1. **Bulk expression matrix**  
   Construct a bulk cohort expression matrix **A_bulk** (rows = patient samples, columns = genes).  
   The bulk matrix is **normalized / centered / scaled** to ensure comparable gene-wise signal.

2. **Latent feature mapping via SVD**  
   Decompose the normalized bulk matrix using SVD:

   <p align="center">
     <picture>
       <source media="(prefers-color-scheme: dark)"
               srcset="https://latex.codecogs.com/svg.image?\color{white}{A_{\text{bulk}}=U_{\text{bulk}}\Sigma_{\text{bulk}}V^{\top}}">
       <img src="https://latex.codecogs.com/svg.image?A_{\text{bulk}}=U_{\text{bulk}}\Sigma_{\text{bulk}}V^{\top}"
            alt="A_bulk = U_bulk Σ_bulk V^T">
     </picture>
   </p>

   * **U_bulk**: sample scores (rows = patients, columns = latent factors)  
   * **Σ_bulk**: diagonal matrix of singular values  
   * **V**: gene loadings (rows = genes, columns = latent factors)

   Define the bulk latent representation (patient-by-factor embedding):

   <p align="center">
     <picture>
       <source media="(prefers-color-scheme: dark)"
               srcset="https://latex.codecogs.com/svg.image?\color{white}{Z_{\text{bulk}}=U_{\text{bulk}}\Sigma_{\text{bulk}}}">
       <img src="https://latex.codecogs.com/svg.image?Z_{\text{bulk}}=U_{\text{bulk}}\Sigma_{\text{bulk}}"
            alt="Z_bulk = U_bulk Σ_bulk">
     </picture>
   </p>

3. **Train mutation-prediction models in latent space**  
   For each mutation / gene-of-interest, train a supervised ML model to predict mutation presence **Y** from **Z_bulk**.  
   This yields one (or multiple) mutation-specific classifiers operating on the learned latent factors.

---

### 2) Project scRNA-seq into the bulk-derived latent space and predict mutations per cell (Panels d–f)

4. **Single-cell expression matrix**  
   Construct a single-cell expression matrix **A_sc** (rows = single cells, columns = genes).

5. **Normalize scRNA-seq using bulk-derived parameters**  
   Apply the *same* gene-wise normalization / centering / scaling learned from the bulk cohort to obtain **A′_sc**.  
   (This alignment step makes the projection comparable across bulk and single-cell.)

6. **Project cells into latent space and infer mutation probabilities**  
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

### Inputs / Outputs (at a glance)

**Inputs**
- Bulk RNA-seq expression + mutation labels (per sample)
- scRNA-seq expression (same gene universe or mapped to it)

**Outputs**
- Bulk latent space (**Z_bulk**) and learned gene loadings (**V**)
- Trained mutation classifiers (bulk-supervised)
- Single-cell latent space (**Z_sc**) and **per-cell mutation probability estimates**
