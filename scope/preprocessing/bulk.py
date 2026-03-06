"""Bulk RNA-seq preprocessing: normalization, centering, and scaling.

All transformers follow the scikit-learn fit / transform / fit_transform API
and store learned parameters in their attributes so they can be applied
identically to single-cell data later.

Changes from v0.1.1
-------------------
* ``BulkPreprocessor`` now accepts ``already_log_transformed``,
  ``min_samples_expressed``, ``min_expression``, ``gene_blacklist``,
  ``auto_remove_mito``, ``auto_remove_ribo``, ``run_hvg``, ``n_hvg``,
  ``batch_key``, and ``batch_method`` parameters — covering the full range
  of "already preprocessed" to "raw counts" input scenarios.
* The ``transform`` path respects all learned masks so sc data is aligned
  identically to the bulk fit.
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
# Optional deps
# ---------------------------------------------------------------------------
try:
    import importlib.util as _ilu

    _HAS_SCANPY = _ilu.find_spec("scanpy") is not None
except Exception:
    _HAS_SCANPY = False

try:
    from combat.pycombat import pycombat as _pycombat

    _HAS_COMBAT = True
except ImportError:
    _HAS_COMBAT = False

try:
    import harmonypy as _harmonypy

    _HAS_HARMONY = True
except ImportError:
    _HAS_HARMONY = False

# ---------------------------------------------------------------------------
# Count normalisation
# ---------------------------------------------------------------------------

NormMethod = Literal["cpm", "tpm", "median_ratio", "tmm", "none"]
BatchMethod = Literal["combat", "harmony"]


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

    def fit(self, adata: AnnData, y=None) -> BulkNormalizer:
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
        X = self._get_X(adata).copy()
        X = self._normalise(X)
        if self.log1p:
            np.log1p(X, out=X)
        adata = adata.copy()
        self._set_X(adata, X)
        return adata

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
        lib_sizes = X.sum(axis=1)
        lib_sizes = np.where(lib_sizes == 0, 1.0, lib_sizes)
        ref_idx = np.argmin(np.abs(lib_sizes - np.quantile(lib_sizes, 0.75)))
        ref = X[ref_idx]
        sf = np.ones(X.shape[0])
        for i, row in enumerate(X):
            mask = (row > 0) & (ref > 0)
            if mask.sum() < 10:
                continue
            M = np.log2(row[mask] / lib_sizes[i]) - np.log2(
                ref[mask] / lib_sizes[ref_idx]
            )
            A = 0.5 * (
                np.log2(row[mask] / lib_sizes[i])
                + np.log2(ref[mask] / lib_sizes[ref_idx])
            )
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
    """Gene-wise centering and/or scaling of a bulk expression matrix."""

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
        log.info("BulkScaler fitted (center=%s, scale=%s).", self.center, self.scale)
        return self

    def transform(self, adata: AnnData, y=None) -> AnnData:
        X = self._get_X(adata)
        X_scaled = self._scaler.transform(X).astype(np.float32)
        adata = adata.copy()
        self._set_X(adata, X_scaled)
        return adata

    @property
    def mean_(self) -> np.ndarray:
        return self._scaler.mean_

    @property
    def scale_(self) -> np.ndarray:
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
# Convenience pipeline wrapper  (EXTENDED)
# ---------------------------------------------------------------------------


class BulkPreprocessor(BaseEstimator, TransformerMixin):
    """Combined normalisation + centering/scaling for bulk RNA-seq.

    Extends the v0.1.1 implementation with explicit flags for every common
    "already preprocessed" starting point, plus optional HVG selection,
    gene filtering, and batch correction.

    Parameters
    ----------
    norm_method:
        Passed to :class:`BulkNormalizer`.  Set ``"none"`` if data are already
        library-size normalised.
    log1p:
        Apply log(x+1) after normalisation.  Set ``False`` if the matrix is
        already log-transformed (e.g. downloaded as log2-TPM from GEO).
    center, scale:
        Gene-wise centering / scaling (always applied — these are required for
        SVD to be meaningful).
    gene_lengths:
        Required for ``norm_method="tpm"``.
    layer_in / layer_out:
        AnnData layers for I/O.

    Preprocessing skip flags
    ~~~~~~~~~~~~~~~~~~~~~~~~
    already_log_transformed:
        If ``True``, the ``log1p`` step is silently skipped regardless of the
        ``log1p`` parameter.  Useful when ``norm_method="none"`` and the input
        is already log-normalised.
    min_samples_expressed:
        Remove genes expressed (> *min_expression*) in fewer than this many
        samples.  ``0`` → no filter.
    min_expression:
        Expression threshold for the above filter.
    gene_blacklist:
        Gene symbols to unconditionally exclude before fitting.
    auto_remove_mito:
        Remove genes whose symbol starts with ``MT-`` (human) or ``mt-``
        (mouse).
    auto_remove_ribo:
        Remove ribosomal protein genes (``^RPS`` / ``^RPL``).

    HVG selection
    ~~~~~~~~~~~~~
    run_hvg:
        Select highly variable genes before fitting the scaler.  Requires
        ``scanpy``.
    n_hvg:
        Number of HVGs to retain.
    hvg_flavor:
        Scanpy HVG flavor: ``"seurat"``, ``"seurat_v3"``, or
        ``"cell_ranger"``.

    Batch correction
    ~~~~~~~~~~~~~~~~
    batch_key:
        Column in ``adata.obs`` used to identify batch labels.  If ``None``,
        batch correction is skipped.
    batch_method:
        ``"combat"`` (requires ``combat`` package) or ``"harmony"`` (requires
        ``harmonypy``).
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
        # --- skip flags ---
        already_log_transformed: bool = False,
        min_samples_expressed: int = 0,
        min_expression: float = 0.0,
        gene_blacklist: list[str] | None = None,
        auto_remove_mito: bool = False,
        auto_remove_ribo: bool = False,
        # --- HVG ---
        run_hvg: bool = False,
        n_hvg: int = 3000,
        hvg_flavor: str = "seurat_v3",
        # --- batch ---
        batch_key: str | None = None,
        batch_method: BatchMethod = "combat",
    ):
        self.norm_method = norm_method
        self.log1p = log1p
        self.center = center
        self.scale = scale
        self.gene_lengths = gene_lengths
        self.layer_in = layer_in
        self.layer_out = layer_out
        self.already_log_transformed = already_log_transformed
        self.min_samples_expressed = min_samples_expressed
        self.min_expression = min_expression
        self.gene_blacklist = gene_blacklist or []
        self.auto_remove_mito = auto_remove_mito
        self.auto_remove_ribo = auto_remove_ribo
        self.run_hvg = run_hvg
        self.n_hvg = n_hvg
        self.hvg_flavor = hvg_flavor
        self.batch_key = batch_key
        self.batch_method = batch_method

    # ------------------------------------------------------------------
    def fit(self, adata: AnnData, y=None) -> BulkPreprocessor:
        # ── gene filtering (fit stores the kept gene mask) ──────────────
        gene_mask = self._build_gene_mask(adata)
        self._gene_mask_: np.ndarray = gene_mask
        self._kept_genes_: list[str] = list(np.array(adata.var_names)[gene_mask])
        adata = adata[:, gene_mask].copy()

        # ── normalisation ───────────────────────────────────────────────
        effective_log1p = self.log1p and not self.already_log_transformed
        self.normalizer_ = BulkNormalizer(
            method=self.norm_method,
            log1p=effective_log1p,
            gene_lengths=self.gene_lengths,
            layer_in=self.layer_in,
            layer_out=self.layer_out,
        )
        adata_normed = self.normalizer_.fit_transform(adata)

        # ── batch correction ────────────────────────────────────────────
        if self.batch_key is not None:
            adata_normed = self._batch_correct(adata_normed)

        # ── HVG selection ───────────────────────────────────────────────
        hvg_mask: np.ndarray | None = None
        if self.run_hvg:
            hvg_mask = self._select_hvg_mask(adata_normed)
            self._hvg_mask_: np.ndarray = hvg_mask
            adata_normed = adata_normed[:, hvg_mask].copy()
        else:
            self._hvg_mask_ = np.ones(adata_normed.n_vars, dtype=bool)

        # ── scaler ──────────────────────────────────────────────────────
        self.scaler_ = BulkScaler(
            center=self.center,
            scale=self.scale,
            layer_in=self.layer_out,
            layer_out=self.layer_out,
        )
        self.scaler_.fit(adata_normed)
        log.info(
            "BulkPreprocessor fitted: %d genes → %d after filtering%s.",
            adata.n_vars,
            adata_normed.n_vars,
            f" → {hvg_mask.sum()} HVGs" if hvg_mask is not None else "",
        )
        return self

    def transform(self, adata: AnnData, y=None) -> AnnData:
        # 1. Subset to the genes seen during fit (zero-pads missing genes)
        adata = self._apply_gene_mask(adata)
        # 2. Normalise using the already-fitted normalizer
        #    (it was built with the correct effective log1p flag during fit)
        adata = self.normalizer_.transform(adata)
        # 3. Apply HVG mask if used during fit
        if self.run_hvg and hasattr(self, "_hvg_mask_"):
            adata = adata[:, self._hvg_mask_].copy()
        # 4. Scale using the already-fitted scaler
        adata = self.scaler_.transform(adata)
        return adata

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_gene_mask(self, adata: AnnData) -> np.ndarray:
        """Return boolean mask of genes to keep."""
        n_genes = adata.n_vars
        mask = np.ones(n_genes, dtype=bool)
        gene_names = np.array(adata.var_names)

        # Blacklist
        bl = {g.upper() for g in self.gene_blacklist}
        for i, g in enumerate(gene_names):
            gu = g.upper()
            if gu in bl:
                mask[i] = False
            if self.auto_remove_mito and gu.startswith("MT-"):
                mask[i] = False
            if self.auto_remove_ribo and (gu.startswith("RPS") or gu.startswith("RPL")):
                mask[i] = False

        # Expression filter
        if self.min_samples_expressed > 0:
            X = self._raw_X(adata)
            expressed = (self.min_expression < X).sum(axis=0)
            mask &= expressed >= self.min_samples_expressed

        n_removed = (~mask).sum()
        if n_removed:
            log.info("Gene filter: removed %d / %d genes.", n_removed, n_genes)
        return mask

    def _apply_gene_mask(self, adata: AnnData) -> AnnData:
        """Subset adata to kept genes; zero-pad if a gene is absent."""
        import anndata as ad

        kept = self._kept_genes_
        present = set(adata.var_names)
        missing = [g for g in kept if g not in present]
        if missing:
            log.warning(
                "%d / %d bulk genes absent from new data (zero-padded).",
                len(missing),
                len(kept),
            )
            # Build zero-padded matrix
            X_src = self._raw_X(adata)
            src_idx = {g: i for i, g in enumerate(adata.var_names)}
            X_out = np.zeros((adata.n_obs, len(kept)), dtype=np.float32)
            for j, g in enumerate(kept):
                if g in src_idx:
                    X_out[:, j] = X_src[:, src_idx[g]]
            adata_aligned = ad.AnnData(X=X_out, obs=adata.obs.copy())
            adata_aligned.var_names = kept
            return adata_aligned
        return adata[:, kept].copy()

    def _batch_correct(self, adata: AnnData) -> AnnData:
        if self.batch_key not in adata.obs.columns:
            log.warning(
                "batch_key '%s' not in adata.obs — skipping batch correction.",
                self.batch_key,
            )
            return adata

        if self.batch_method == "combat":
            if not _HAS_COMBAT:
                log.warning("pycombat not installed (pip install combat). Skipping.")
                return adata
            import pandas as pd

            X = self._raw_X(adata)
            batch = adata.obs[self.batch_key].values
            df = pd.DataFrame(X.T, index=list(adata.var_names))
            corrected = _pycombat(df, batch)
            adata = adata.copy()
            adata.X = corrected.values.T.astype(np.float32)
            log.info("ComBat batch correction applied.")

        elif self.batch_method == "harmony":
            if not _HAS_HARMONY:
                log.warning("harmonypy not installed. Skipping batch correction.")
                return adata
            X = self._raw_X(adata)
            ho = _harmonypy.run_harmony(X, adata.obs, self.batch_key)
            adata = adata.copy()
            adata.X = ho.Z_corr.T.astype(np.float32)
            log.info("HarmonyPy batch correction applied.")
        return adata

    def _select_hvg_mask(self, adata: AnnData) -> np.ndarray:
        if not _HAS_SCANPY:
            log.warning("scanpy not installed — falling back to top-variance HVG.")
            X = self._raw_X(adata)
            variances = X.var(axis=0)
            idx = np.argsort(variances)[::-1][: self.n_hvg]
            mask = np.zeros(adata.n_vars, dtype=bool)
            mask[idx] = True
            return mask

        import scanpy as sc

        adata_tmp = adata.copy()
        sc.pp.highly_variable_genes(
            adata_tmp, n_top_genes=self.n_hvg, flavor=self.hvg_flavor
        )
        mask = adata_tmp.var["highly_variable"].values
        log.info("HVG selection: %d / %d genes.", mask.sum(), adata.n_vars)
        return mask

    @staticmethod
    def _raw_X(adata: AnnData) -> np.ndarray:
        X = adata.X
        if sp.issparse(X):
            X = X.toarray()
        return X.astype(np.float64)
