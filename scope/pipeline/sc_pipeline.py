"""SingleCellPipeline: orchestrates Phase 2 of scOPE.

Phase 2 takes the fitted :class:`~scope.pipeline.bulk_pipeline.BulkPipeline`
and projects single-cell data into the bulk-derived latent space to infer
per-cell mutation probabilities.

Typical usage
-------------
>>> sc_pipe = SingleCellPipeline(bulk_pipeline=bulk_pipe)
>>> adata_sc = sc_pipe.transform(adata_sc)
>>> # adata_sc.obs now contains mutation_prob_* columns
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np
from anndata import AnnData
from sklearn.base import BaseEstimator

from scope.pipeline.bulk_pipeline import BulkPipeline
from scope.preprocessing.alignment import BulkSCAligner
from scope.preprocessing.single_cell import SingleCellPreprocessor
from scope.utils.gene_utils import subset_to_shared_genes
from scope.utils.logging import get_logger
from scope.visualization.embeddings import compute_umap, plot_mutation_probabilities

log = get_logger(__name__)

PathLike = Union[str, Path]


class SingleCellPipeline(BaseEstimator):
    """End-to-end single-cell projection and mutation inference phase.

    Parameters
    ----------
    bulk_pipeline:
        A fitted :class:`~scope.pipeline.bulk_pipeline.BulkPipeline` instance.
    sc_filter_strategy:
        QC filter for single-cell data: ``"both"``, ``"min_counts"``,
        ``"min_genes"``, or ``"none"``.
    sc_min_counts:
        Minimum total counts per cell.
    sc_min_genes:
        Minimum genes detected per cell.
    sc_target_sum:
        Library-size normalisation target.
    alignment_method:
        Bulk–sc expression alignment: ``"z_score_bulk"`` (default),
        ``"moment_matching"``, ``"quantile"``, or ``"none"``.
    clip_percentile:
        Clip aligned sc values to this bulk percentile.
    layer:
        AnnData layer to use; ``None`` → ``adata.X``.
    add_to_obs:
        If ``True``, write mutation probabilities directly into
        ``adata_sc.obs`` (keys: ``mutation_prob_<gene>``).
    """

    def __init__(
        self,
        bulk_pipeline: BulkPipeline,
        sc_filter_strategy: str = "both",
        sc_min_counts: int = 200,
        sc_min_genes: int = 200,
        sc_target_sum: float = 1e4,
        alignment_method: str = "z_score_bulk",
        clip_percentile: float | None = 99.0,
        layer: str | None = None,
        add_to_obs: bool = True,
    ):
        self.bulk_pipeline = bulk_pipeline
        self.sc_filter_strategy = sc_filter_strategy
        self.sc_min_counts = sc_min_counts
        self.sc_min_genes = sc_min_genes
        self.sc_target_sum = sc_target_sum
        self.alignment_method = alignment_method
        self.clip_percentile = clip_percentile
        self.layer = layer
        self.add_to_obs = add_to_obs

    # ------------------------------------------------------------------
    def fit(self, adata_bulk_pp: AnnData, adata_sc: AnnData) -> SingleCellPipeline:
        """Fit sc preprocessing and bulk–sc aligner.

        Parameters
        ----------
        adata_bulk_pp:
            Preprocessed (normalised, centred, scaled) bulk AnnData — the
            same object produced during ``BulkPipeline.fit``.
        adata_sc:
            Raw (or pre-filtered) single-cell AnnData.

        Returns
        -------
        self
        """
        log.info("=== SingleCellPipeline.fit ===")

        # ── SC QC / normalisation ──────────────────────────────────────
        self.sc_preprocessor_ = SingleCellPreprocessor(
            filter_strategy=self.sc_filter_strategy,
            min_counts=self.sc_min_counts,
            min_genes=self.sc_min_genes,
            target_sum=self.sc_target_sum,
            log1p=True,
            scale=False,  # Alignment handles scaling
            layer_in=self.layer,
            layer_out=self.layer,
        )
        adata_sc_pp = self.sc_preprocessor_.fit_transform(adata_sc)

        # ── Gene universe alignment ────────────────────────────────────
        adata_bulk_sub, adata_sc_sub = subset_to_shared_genes(adata_bulk_pp, adata_sc_pp)

        # ── Moment matching / alignment ────────────────────────────────
        self.aligner_ = BulkSCAligner(
            method=self.alignment_method,
            clip_percentile=self.clip_percentile,
            layer_bulk=self.layer,
            layer_sc_in=self.layer,
            layer_sc_out=self.layer,
        )
        self.aligner_.fit(adata_bulk_sub)
        self.is_fitted_ = True
        log.info("SingleCellPipeline.fit complete.")
        return self

    # ------------------------------------------------------------------
    def transform(self, adata_sc: AnnData) -> AnnData:
        """Project single-cell data into bulk latent space and infer mutations.

        This is the main inference method. It:

        1. QC-filters cells.
        2. Normalises and log-transforms sc expression.
        3. Subsets to the bulk gene universe.
        4. Aligns sc distribution to bulk.
        5. Projects cells into the bulk latent space.
        6. Applies trained mutation classifiers to get per-cell probabilities.
        7. Writes probabilities to ``adata_sc.obs``.

        Parameters
        ----------
        adata_sc:
            Raw (or minimally preprocessed) single-cell AnnData.

        Returns
        -------
        AnnData
            Copy with:

            * ``obsm[bulk_pipeline.obsm_key_]`` — latent embedding.
            * ``obs[mutation_prob_<gene>]`` — per-cell mutation probabilities.
        """
        log.info("=== SingleCellPipeline.transform ===")
        bp = self.bulk_pipeline

        # ── 1. SC QC / normalisation ───────────────────────────────────
        if hasattr(self, "sc_preprocessor_"):
            adata_pp = self.sc_preprocessor_.transform(adata_sc)
        else:
            # Fast path: fit+transform in one shot
            sc_prep = SingleCellPreprocessor(
                filter_strategy=self.sc_filter_strategy,
                min_counts=self.sc_min_counts,
                min_genes=self.sc_min_genes,
                target_sum=self.sc_target_sum,
                log1p=True,
                scale=False,
                layer_in=self.layer,
                layer_out=self.layer,
            )
            adata_pp = sc_prep.fit_transform(adata_sc)

        # ── 2. Subset to bulk gene universe ───────────────────────────
        bulk_genes = bp.gene_names_
        sc_genes = list(adata_pp.var_names)
        shared = [g for g in bulk_genes if g in set(sc_genes)]
        missing = len(bulk_genes) - len(shared)
        if missing > 0:
            log.warning(
                "%d / %d bulk genes absent from sc data (will be zero-padded).",
                missing, len(bulk_genes),
            )

        # Re-index sc to exact bulk gene order, zero-filling gaps
        import scipy.sparse as sp
        X_sc = adata_pp.X
        if sp.issparse(X_sc):
            X_sc = X_sc.toarray()
        gene_idx = {g: i for i, g in enumerate(sc_genes)}
        X_aligned = np.zeros((adata_pp.n_obs, len(bulk_genes)), dtype=np.float32)
        for j, g in enumerate(bulk_genes):
            if g in gene_idx:
                X_aligned[:, j] = X_sc[:, gene_idx[g]]

        import anndata as ad
        adata_aligned = ad.AnnData(X=X_aligned)
        adata_aligned.obs_names = list(adata_pp.obs_names)
        adata_aligned.var_names = bulk_genes
        # Copy obs metadata
        adata_aligned.obs = adata_pp.obs.copy()

        # ── 3. Bulk–SC alignment ───────────────────────────────────────
        if hasattr(self, "aligner_"):
            adata_aligned = self.aligner_.transform(adata_aligned)

        # ── 4. Project into bulk latent space ─────────────────────────
        adata_emb = bp.decomposer_.transform(adata_aligned)
        Z_sc = adata_emb.obsm[bp.obsm_key_]
        log.info(
            "Projected %d cells into %d-D latent space.", Z_sc.shape[0], Z_sc.shape[1]
        )

        # ── 5. Mutation probabilities ──────────────────────────────────
        prob_df = bp.classifier_set_.predict_proba(Z_sc)
        log.info(
            "Inferred probabilities for %d mutations.", len(prob_df.columns)
        )

        # ── 6. Write back to AnnData ───────────────────────────────────
        result = adata_emb.copy()
        if self.add_to_obs:
            prob_df.index = list(result.obs_names)
            for col in prob_df.columns:
                result.obs[col] = prob_df[col].values
        result.uns["scope_mutations"] = list(bp.classifier_set_.classifiers_.keys())
        result.uns["scope_decomposition"] = bp.decomposition
        result.uns["scope_n_components"] = bp.n_components
        return result

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------

    def run_umap(
        self,
        adata_sc: AnnData,
        n_neighbors: int = 15,
        min_dist: float = 0.3,
    ) -> AnnData:
        """Run UMAP on the latent embedding and return updated AnnData."""
        obsm_key = self.bulk_pipeline.obsm_key_
        return compute_umap(
            adata_sc,
            obsm_key=obsm_key,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
        )

    def plot_mutations(self, adata_sc: AnnData, **kwargs):
        """Shortcut to :func:`~scope.visualization.embeddings.plot_mutation_probabilities`."""
        return plot_mutation_probabilities(adata_sc, **kwargs)
