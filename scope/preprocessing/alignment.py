"""Bulk–single-cell distribution alignment (moment matching).

Before projecting single-cell data into the bulk-derived latent space the two
datasets must share a comparable feature space.  This module provides several
strategies for aligning marginal distributions gene by gene.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
import scipy.sparse as sp
from anndata import AnnData
from sklearn.base import BaseEstimator, TransformerMixin

from scope.utils.logging import get_logger

log = get_logger(__name__)

AlignMethod = Literal[
    "moment_matching",     # match mean and std of each gene
    "quantile",            # quantile normalisation to bulk reference
    "z_score_bulk",        # subtract bulk mean, divide by bulk std
    "none",
]


class BulkSCAligner(BaseEstimator, TransformerMixin):
    """Align single-cell expression to bulk gene-wise statistics.

    Learns gene-wise statistics from bulk data (via ``fit``) then applies the
    alignment to single-cell data (via ``transform``).  This ensures the sc
    projection into the bulk latent space is on the same scale.

    Parameters
    ----------
    method:
        Alignment strategy:

        * ``"moment_matching"`` — shift and scale sc gene distributions to
          match bulk mean and variance.
        * ``"quantile"``        — quantile normalise sc genes to the bulk
          reference distribution (gene-wise).
        * ``"z_score_bulk"``    — subtract bulk mean, divide by bulk std
          (i.e., apply the bulk z-score transformation to sc data).
        * ``"none"``            — no alignment (pass-through).
    n_quantiles:
        Number of quantile breakpoints for ``method="quantile"``.
    clip_percentile:
        After alignment, clip sc values to ``[bulk_pX_lo, bulk_pX_hi]`` where
        X is this percentile. ``None`` → no clipping.
    layer_bulk:
        Layer of bulk AnnData to learn from. ``None`` → ``adata.X``.
    layer_sc_in:
        Layer of sc AnnData to read. ``None`` → ``adata.X``.
    layer_sc_out:
        Layer of sc AnnData to write into. ``None`` → ``adata.X``.
    """

    def __init__(
        self,
        method: AlignMethod = "z_score_bulk",
        n_quantiles: int = 1000,
        clip_percentile: float | None = 99.0,
        layer_bulk: str | None = None,
        layer_sc_in: str | None = None,
        layer_sc_out: str | None = None,
    ):
        self.method = method
        self.n_quantiles = n_quantiles
        self.clip_percentile = clip_percentile
        self.layer_bulk = layer_bulk
        self.layer_sc_in = layer_sc_in
        self.layer_sc_out = layer_sc_out

    # ------------------------------------------------------------------
    def fit(self, adata_bulk: AnnData, y=None) -> BulkSCAligner:
        """Learn gene-wise statistics from *adata_bulk*.

        Parameters
        ----------
        adata_bulk:
            Preprocessed bulk AnnData (samples × genes).
        """
        X_bulk = self._get_X(adata_bulk, self.layer_bulk)
        self.bulk_mean_ = X_bulk.mean(axis=0)           # (n_genes,)
        self.bulk_std_ = X_bulk.std(axis=0)
        self.bulk_std_ = np.where(self.bulk_std_ == 0, 1.0, self.bulk_std_)

        if self.method == "quantile":
            quantile_probs = np.linspace(0, 1, self.n_quantiles)
            self._bulk_quantiles_ = np.quantile(X_bulk, quantile_probs, axis=0)

        if self.clip_percentile is not None:
            lo = (100.0 - self.clip_percentile) / 2.0
            hi = 100.0 - lo
            self._clip_lo_ = np.percentile(X_bulk, lo, axis=0)
            self._clip_hi_ = np.percentile(X_bulk, hi, axis=0)

        self.n_genes_ = X_bulk.shape[1]
        log.info("BulkSCAligner fitted (method=%s, n_genes=%d).", self.method, self.n_genes_)
        return self

    def transform(self, adata_sc: AnnData, y=None) -> AnnData:
        """Apply bulk alignment to *adata_sc*.

        Parameters
        ----------
        adata_sc:
            Single-cell AnnData with the *same* gene ordering as the bulk data
            used in ``fit``.  Use :func:`~scope.utils.gene_utils.subset_to_shared_genes`
            to ensure this.

        Returns
        -------
        AnnData
            Modified copy with aligned expression in *layer_sc_out*.
        """
        X_sc = self._get_X(adata_sc, self.layer_sc_in).copy()

        if self.method == "moment_matching":
            sc_mean = X_sc.mean(axis=0)
            sc_std = X_sc.std(axis=0)
            sc_std = np.where(sc_std == 0, 1.0, sc_std)
            X_sc = (X_sc - sc_mean) / sc_std * self.bulk_std_ + self.bulk_mean_

        elif self.method == "z_score_bulk":
            X_sc = (X_sc - self.bulk_mean_) / self.bulk_std_

        elif self.method == "quantile":
            X_sc = self._quantile_align(X_sc)

        elif self.method == "none":
            pass

        else:
            raise ValueError(f"Unknown alignment method: {self.method!r}")

        if self.clip_percentile is not None:
            X_sc = np.clip(X_sc, self._clip_lo_, self._clip_hi_)

        adata_sc = adata_sc.copy()
        layer_out = self.layer_sc_out
        if layer_out:
            adata_sc.layers[layer_out] = X_sc.astype(np.float32)
        else:
            adata_sc.X = X_sc.astype(np.float32)
        return adata_sc

    # ------------------------------------------------------------------
    def _quantile_align(self, X_sc: np.ndarray) -> np.ndarray:
        """Gene-wise quantile alignment of sc to bulk distribution."""
        n_q = self._bulk_quantiles_.shape[0]
        q_probs = np.linspace(0, 1, n_q)
        out = np.empty_like(X_sc)
        for j in range(X_sc.shape[1]):
            sc_col = X_sc[:, j]
            # Compute quantiles of sc column
            sc_quantiles = np.quantile(sc_col, q_probs)
            # Map sc values → uniform → bulk quantile (interpolation)
            out[:, j] = np.interp(sc_col, sc_quantiles, self._bulk_quantiles_[:, j])
        return out

    @staticmethod
    def _get_X(adata: AnnData, layer: str | None) -> np.ndarray:
        X = adata.layers[layer] if layer else adata.X
        if sp.issparse(X):
            X = X.toarray()
        return X.astype(np.float64)
