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

<ol>
  <li><b>Learn latent factors from bulk RNA-seq and train mutation classifiers (Panels a–c)</b>
    <ol type="i">
      <li><b>Bulk expression matrix</b><br>
        Construct a bulk cohort expression matrix <b>A_bulk</b> (rows = patient samples, columns = genes).<br>
        The bulk matrix is <b>normalized / centered / scaled</b> to ensure comparable gene-wise signal.
      </li>

      <li><b>Latent feature mapping via SVD</b><br>
        Decompose the normalized bulk matrix using SVD:
        <p align="center">…your picture/img equation…</p>
      </li>

      <li><b>Train mutation-prediction models in latent space</b><br>
        For each mutation / gene-of-interest, train a supervised ML model to predict mutation presence <b>Y</b> from <b>Z_bulk</b>.
      </li>
    </ol>
  </li>

  <li><b>Project scRNA-seq into the bulk-derived latent space and predict mutations per cell (Panels d–f)</b>
    <ol type="i">
      <li><b>Single-cell expression matrix</b><br>
        Construct a single-cell expression matrix <b>A_sc</b> (rows = single cells, columns = genes).
      </li>

      <li><b>Normalize scRNA-seq using bulk-derived parameters</b><br>
        Apply the same gene-wise normalization / centering / scaling learned from the bulk cohort to obtain <b>A′_sc</b>.
      </li>

      <li><b>Project cells into latent space and infer mutation probabilities</b><br>
        Use the bulk-derived gene loadings <b>V</b> to compute:
        <p align="center">…your picture/img equation…</p>
      </li>
    </ol>
  </li>
</ol>

---

### 2) Project scRNA-seq into the bulk-derived latent space and predict mutations per cell (Panels d–f)

 i. **Single-cell expression matrix**  
   Construct a single-cell expression matrix **A_sc** (rows = single cells, columns = genes).

 ii. **Normalize scRNA-seq using bulk-derived parameters**  
   Apply the *same* gene-wise normalization / centering / scaling learned from the bulk cohort to obtain **A′_sc**.  
   (This alignment step makes the projection comparable across bulk and single-cell.)

 iii. **Project cells into latent space and infer mutation probabilities**  
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
