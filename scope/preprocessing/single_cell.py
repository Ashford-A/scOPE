"""Single-cell RNA-seq preprocessing utilities.

Provides AnnData-aware wrappers around standard scRNA-seq normalisation steps
and integrates with scOPE's bulk-derived parameter store so that the same
transformations can be applied consistently across datasets.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
import scipy.sparse as sp
from anndata import AnnData
from sklearn.base import BaseEstimator, TransformerMixin

from scope.utils.logging import get_logger

log = get_logger(__name__)

FilterStrategy = Literal["min_counts", "min_genes", "both", "none"]


class SingleCellPreprocessor(BaseEstimator, TransformerMixin):
    """Standard scRNA-seq preprocessing: filter → normalise → log1p → scale.

    Parameters
    ----------
    filter_strategy:
        Which quality-control filter to apply:

        * ``"min_counts"``  — keep cells with at least *min_counts* total UMIs.
        * ``"min_genes"``   — keep cells expressing at least *min_genes* genes.
        * ``"both"``        — apply both filters.
        * ``"none"``        — skip filtering.
    min_counts:
        Minimum total counts per cell (used when ``filter_strategy`` includes
        ``"min_counts"``).
    min_genes:
        Minimum genes detected per cell (used when ``filter_strategy`` includes
        ``"min_genes"``).
    max_counts:
        Maximum total counts per cell (doublet proxy). ``None`` → no upper
        bound.
    target_sum:
        Library size target for per-cell normalisation (e.g. 1e4 or 1e6).
    log1p:
        Apply log(x + 1) after library-size normalisation.
    scale:
        Apply gene-wise z-scoring after log-normalisation.
    max_value:
        Clip scaled values to ``[-max_value, max_value]``. ``None`` → no
        clipping. A value of 10 is commonly used.
    layer_in:
        AnnData layer to read; ``None`` → ``adata.X``.
    layer_out:
        AnnData layer to store result; ``None`` → ``adata.X``.
    """

    def __init__(
        self,
        filter_strategy: FilterStrategy = "both",
        min_counts: int = 200,
        min_genes: int = 200,
        max_counts: int | None = None,
        target_sum: float = 1e4,
        log1p: bool = True,
        scale: bool = False,
        max_value: float | None = 10.0,
        layer_in: str | None = None,
        layer_out: str | None = None,
    ):
        self.filter_strategy = filter_strategy
        self.min_counts = min_counts
        self.min_genes = min_genes
        self.max_counts = max_counts
        self.target_sum = target_sum
        self.log1p = log1p
        self.scale = scale
        self.max_value = max_value
        self.layer_in = layer_in
        self.layer_out = layer_out

    # ------------------------------------------------------------------
    def fit(self, adata: AnnData, y=None) -> SingleCellPreprocessor:
        """Learn gene-wise statistics (mean, std) if scaling is requested."""
        if self.scale:
            filtered = self._apply_filter(adata)
            X = self._get_X(filtered)
            X = self._normalise(X)
            self._scale_mean_ = X.mean(axis=0)
            _std = X.std(axis=0)
            self._scale_std_ = np.where(_std == 0, 1.0, _std)
        log.info("SingleCellPreprocessor fitted.")
        return self

    def transform(self, adata: AnnData, y=None) -> AnnData:
        """Apply preprocessing pipeline, returning a new AnnData."""
        adata = self._apply_filter(adata)
        X = self._get_X(adata)
        X = self._normalise(X)
        if self.scale:
            X = (X - self._scale_mean_) / self._scale_std_
            if self.max_value is not None:
                X = np.clip(X, -self.max_value, self.max_value)
        adata = adata.copy()
        self._set_X(adata, X.astype(np.float32))
        return adata

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _apply_filter(self, adata: AnnData) -> AnnData:
        if self.filter_strategy == "none":
            return adata
        X = self._get_X(adata)
        cell_counts = X.sum(axis=1)
        gene_counts = (X > 0).sum(axis=1)
        mask = np.ones(adata.n_obs, dtype=bool)
        if self.filter_strategy in ("min_counts", "both"):
            mask &= cell_counts >= self.min_counts
        if self.filter_strategy in ("min_genes", "both"):
            mask &= gene_counts >= self.min_genes
        if self.max_counts is not None:
            mask &= cell_counts <= self.max_counts
        n_removed = (~mask).sum()
        if n_removed > 0:
            log.info("QC filter: removed %d / %d cells.", n_removed, adata.n_obs)
        return adata[mask].copy()

    def _normalise(self, X: np.ndarray) -> np.ndarray:
        lib = X.sum(axis=1, keepdims=True)
        lib = np.where(lib == 0, 1.0, lib)
        X = X / lib * self.target_sum
        if self.log1p:
            X = np.log1p(X)
        return X

    def _get_X(self, adata: AnnData) -> np.ndarray:
        X = adata.layers[self.layer_in] if self.layer_in else adata.X
        if sp.issparse(X):
            X = X.toarray()
        return X.astype(np.float64)

    def _set_X(self, adata: AnnData, X: np.ndarray) -> None:
        if self.layer_out:
            adata.layers[self.layer_out] = X
        else:
            adata.X = X
