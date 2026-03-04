"""Bulk RNA-seq preprocessing: normalization, centering, and scaling.

All transformers follow the scikit-learn fit / transform / fit_transform API
and store learned parameters in their attributes so they can be applied
identically to single-cell data later.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
import scipy.sparse as sp
from anndata import AnnData
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import StandardScaler

from scope.utils.logging import get_logger

log = get_logger(__name__)

# ---------------------------------------------------------------------------
# Count normalisation
# ---------------------------------------------------------------------------

NormMethod = Literal["cpm", "tpm", "median_ratio", "tmm", "none"]


class BulkNormalizer(BaseEstimator, TransformerMixin):
    """Normalise a bulk expression matrix (samples × genes).

    Parameters
    ----------
    method:
        Normalisation strategy:

        * ``"cpm"``          — counts per million.
        * ``"tpm"``          — TPM-style (requires *gene_lengths*).
        * ``"median_ratio"`` — DESeq2-style geometric mean / median ratio.
        * ``"tmm"``          — edgeR-inspired TMM (approximate).
        * ``"none"``         — no library-size normalisation.
    log1p:
        If ``True``, apply log(x + 1) after normalisation.
    target_sum:
        Total counts per sample after CPM-style normalisation. Ignored for
        methods other than ``"cpm"``/``"tpm"``.
    gene_lengths:
        Array of gene lengths in bp, shape ``(n_genes,)``. Required for
        ``method="tpm"``.
    layer_in:
        AnnData layer to read from. ``None`` → ``adata.X``.
    layer_out:
        AnnData layer to write normalised values into. ``None`` → ``adata.X``.
    """

    def __init__(
        self,
        method: NormMethod = "cpm",
        log1p: bool = True,
        target_sum: float = 1e6,
        gene_lengths: np.ndarray | None = None,
        layer_in: str | None = None,
        layer_out: str | None = None,
    ):
        self.method = method
        self.log1p = log1p
        self.target_sum = target_sum
        self.gene_lengths = gene_lengths
        self.layer_in = layer_in
        self.layer_out = layer_out

    # ------------------------------------------------------------------
    def fit(self, adata: AnnData, y=None) -> BulkNormalizer:
        """Learn normalisation parameters from *adata*.

        For ``"median_ratio"`` and ``"tmm"`` the reference is computed here.
        """
        X = self._get_X(adata)
        if self.method == "median_ratio":
            self._size_factors_ = self._median_ratio_sf(X)
        elif self.method == "tmm":
            self._size_factors_ = self._tmm_sf(X)
        else:
            self._size_factors_ = None
        self.n_genes_ = X.shape[1]
        log.info("BulkNormalizer fitted (method=%s).", self.method)
        return self

    def transform(self, adata: AnnData, y=None) -> AnnData:
        """Apply normalisation, returning a modified copy of *adata*."""
        X = self._get_X(adata).copy()
        X = self._normalise(X)
        if self.log1p:
            np.log1p(X, out=X)
        adata = adata.copy()
        self._set_X(adata, X)
        return adata

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _get_X(self, adata: AnnData) -> np.ndarray:
        X = adata.layers[self.layer_in] if self.layer_in else adata.X
        if sp.issparse(X):
            X = X.toarray()
        return X.astype(np.float64)

    def _set_X(self, adata: AnnData, X: np.ndarray) -> None:
        if self.layer_out:
            adata.layers[self.layer_out] = X.astype(np.float32)
        else:
            adata.X = X.astype(np.float32)

    def _normalise(self, X: np.ndarray) -> np.ndarray:
        if self.method == "cpm":
            return self._cpm(X)
        if self.method == "tpm":
            return self._tpm(X)
        if self.method in ("median_ratio", "tmm"):
            return X / self._size_factors_[:, np.newaxis]
        if self.method == "none":
            return X
        raise ValueError(f"Unknown normalisation method: {self.method!r}")

    @staticmethod
    def _cpm(X: np.ndarray, target: float = 1e6) -> np.ndarray:
        lib_sizes = X.sum(axis=1, keepdims=True)
        lib_sizes = np.where(lib_sizes == 0, 1.0, lib_sizes)
        return X / lib_sizes * target

    def _tpm(self, X: np.ndarray) -> np.ndarray:
        if self.gene_lengths is None:
            raise ValueError("gene_lengths must be provided for TPM normalisation.")
        lengths = np.asarray(self.gene_lengths, dtype=np.float64)
        rpk = X / lengths[np.newaxis, :]
        per_sample_sum = rpk.sum(axis=1, keepdims=True)
        per_sample_sum = np.where(per_sample_sum == 0, 1.0, per_sample_sum)
        return rpk / per_sample_sum * self.target_sum

    @staticmethod
    def _median_ratio_sf(X: np.ndarray) -> np.ndarray:
        """DESeq2-style median-of-ratios size factors."""
        log_X = np.log(np.where(X > 0, X, np.nan))
        log_geo_mean = np.nanmean(log_X, axis=0)
        finite_mask = np.isfinite(log_geo_mean)
        if finite_mask.sum() == 0:
            return np.ones(X.shape[0])
        log_ratios = log_X[:, finite_mask] - log_geo_mean[np.newaxis, finite_mask]
        sf = np.exp(np.nanmedian(log_ratios, axis=1))
        sf = np.where(sf == 0, 1.0, sf)
        return sf

    @staticmethod
    def _tmm_sf(X: np.ndarray) -> np.ndarray:
        """Approximate TMM: trim extreme M/A values and compute weighted mean."""
        lib_sizes = X.sum(axis=1)
        lib_sizes = np.where(lib_sizes == 0, 1.0, lib_sizes)
        ref_idx = np.argmin(np.abs(lib_sizes - np.quantile(lib_sizes, 0.75)))
        ref = X[ref_idx]
        sf = np.ones(X.shape[0])
        for i, row in enumerate(X):
            mask = (row > 0) & (ref > 0)
            if mask.sum() < 10:
                continue
            M = np.log2(row[mask] / lib_sizes[i]) - np.log2(ref[mask] / lib_sizes[ref_idx])
            A = 0.5 * (
                np.log2(row[mask] / lib_sizes[i])
                + np.log2(ref[mask] / lib_sizes[ref_idx])
            )
            # Trim top/bottom 30% by M, 5% by A
            m_lo, m_hi = np.quantile(M, [0.30, 0.70])
            a_lo, a_hi = np.quantile(A, [0.05, 0.95])
            keep = (m_lo <= M) & (m_hi >= M) & (a_lo <= A) & (a_hi >= A)
            if keep.sum() < 5:
                continue
            w = 1.0 / (1 / row[mask][keep] + 1 / ref[mask][keep])
            sf[i] = 2 ** (np.sum(w * M[keep]) / np.sum(w))
        return sf


# ---------------------------------------------------------------------------
# Gene-wise centering and scaling
# ---------------------------------------------------------------------------


class BulkScaler(BaseEstimator, TransformerMixin):
    """Gene-wise centering and/or scaling of a bulk expression matrix.

    Wraps :class:`sklearn.preprocessing.StandardScaler` with AnnData-aware
    input/output and stores parameters for downstream projection onto sc data.

    Parameters
    ----------
    center:
        Subtract gene-wise mean.
    scale:
        Divide by gene-wise standard deviation.
    layer_in:
        AnnData layer to read; ``None`` → ``adata.X``.
    layer_out:
        AnnData layer to write; ``None`` → ``adata.X``.
    """

    def __init__(
        self,
        center: bool = True,
        scale: bool = True,
        layer_in: str | None = None,
        layer_out: str | None = None,
    ):
        self.center = center
        self.scale = scale
        self.layer_in = layer_in
        self.layer_out = layer_out

    def fit(self, adata: AnnData, y=None) -> BulkScaler:
        X = self._get_X(adata)
        self._scaler = StandardScaler(with_mean=self.center, with_std=self.scale)
        self._scaler.fit(X)
        log.info(
            "BulkScaler fitted (center=%s, scale=%s).", self.center, self.scale
        )
        return self

    def transform(self, adata: AnnData, y=None) -> AnnData:
        X = self._get_X(adata)
        X_scaled = self._scaler.transform(X).astype(np.float32)
        adata = adata.copy()
        self._set_X(adata, X_scaled)
        return adata

    @property
    def mean_(self) -> np.ndarray:
        """Gene-wise means (shape: n_genes)."""
        return self._scaler.mean_

    @property
    def scale_(self) -> np.ndarray:
        """Gene-wise std (shape: n_genes)."""
        return self._scaler.scale_

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


# ---------------------------------------------------------------------------
# Convenience pipeline wrapper
# ---------------------------------------------------------------------------


class BulkPreprocessor(BaseEstimator, TransformerMixin):
    """Combined normalisation + centering/scaling for bulk RNA-seq.

    Sequentially applies :class:`BulkNormalizer` then :class:`BulkScaler`.

    Parameters
    ----------
    norm_method:
        Passed to :class:`BulkNormalizer`.
    log1p:
        Apply log(x+1) after library-size normalisation.
    center, scale:
        Passed to :class:`BulkScaler`.
    gene_lengths:
        Required for ``norm_method="tpm"``.
    layer_in:
        Input layer. ``None`` → ``adata.X``.
    layer_out:
        Output layer. ``None`` → ``adata.X``.
    """

    def __init__(
        self,
        norm_method: NormMethod = "cpm",
        log1p: bool = True,
        center: bool = True,
        scale: bool = True,
        gene_lengths: np.ndarray | None = None,
        layer_in: str | None = None,
        layer_out: str | None = None,
    ):
        self.norm_method = norm_method
        self.log1p = log1p
        self.center = center
        self.scale = scale
        self.gene_lengths = gene_lengths
        self.layer_in = layer_in
        self.layer_out = layer_out

    def fit(self, adata: AnnData, y=None) -> BulkPreprocessor:
        self.normalizer_ = BulkNormalizer(
            method=self.norm_method,
            log1p=self.log1p,
            gene_lengths=self.gene_lengths,
            layer_in=self.layer_in,
            layer_out=self.layer_out,
        )
        adata_normed = self.normalizer_.fit_transform(adata)
        self.scaler_ = BulkScaler(
            center=self.center,
            scale=self.scale,
            layer_in=self.layer_out,
            layer_out=self.layer_out,
        )
        self.scaler_.fit(adata_normed)
        return self

    def transform(self, adata: AnnData, y=None) -> AnnData:
        adata = self.normalizer_.transform(adata)
        adata = self.scaler_.transform(adata)
        return adata
