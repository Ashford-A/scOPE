"""Single-cell RNA-seq preprocessing utilities.

Changes from v0.1.1
-------------------
* ``max_genes`` filter (upper bound on detected genes per cell — doublet proxy
  complementary to ``max_counts``).
* ``max_mito_pct`` — mitochondrial fraction filter (common scRNA-seq QC step).
* ``auto_flag_mito`` — annotate ``adata.obs['pct_mito']`` even when
  ``filter_strategy='none'``.
* ``run_doublet_detection`` — optional Scrublet-based doublet scoring and
  removal (requires ``scrublet``).
* ``doublet_threshold`` — manual score threshold; ``None`` = automatic.
* ``already_qc_filtered`` — skip all cell-level filtering when ``True``.
* ``already_normalized`` — skip library-size normalisation.
* ``already_log_transformed`` — skip log1p.
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

# Optional deps -----------------------------------------------------------
try:
    import scrublet as _scrublet

    _HAS_SCRUBLET = True
except ImportError:
    _HAS_SCRUBLET = False


class SingleCellPreprocessor(BaseEstimator, TransformerMixin):
    """Standard scRNA-seq preprocessing: filter → normalise → log1p → scale.

    Parameters
    ----------
    filter_strategy:
        QC filter to apply:

        * ``"min_counts"``  — keep cells with ≥ *min_counts* total UMIs.
        * ``"min_genes"``   — keep cells expressing ≥ *min_genes* genes.
        * ``"both"``        — apply both count and gene filters.
        * ``"none"``        — skip cell filtering entirely.
    min_counts:
        Minimum total counts per cell.
    min_genes:
        Minimum detected genes per cell.
    max_counts:
        Maximum total counts per cell (doublet proxy). ``None`` → no upper
        bound.
    max_genes:
        Maximum detected genes per cell (doublet proxy).  ``None`` → no upper
        bound.
    max_mito_pct:
        Maximum mitochondrial fraction (0–100).  Cells exceeding this are
        removed.  ``None`` → no mito filter.  Requires gene symbols to
        have ``MT-`` / ``mt-`` prefix.
    auto_flag_mito:
        If ``True``, compute ``pct_mito`` in ``adata.obs`` even when mito
        filtering is disabled.
    target_sum:
        Library-size normalisation target (1e4 or 1e6 are typical).
    log1p:
        Apply log(x + 1) after library-size normalisation.
    scale:
        Gene-wise z-scoring after log-normalisation.
    max_value:
        Clip scaled values to ``[-max_value, max_value]``. ``None`` → no clip.
    layer_in:
        AnnData layer to read; ``None`` → ``adata.X``.
    layer_out:
        AnnData layer to store result; ``None`` → ``adata.X``.

    Skip flags (for already-preprocessed inputs)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    already_qc_filtered:
        Skip all cell-level QC filters (count, gene, mito, doublet).
    already_normalized:
        Skip library-size normalisation.
    already_log_transformed:
        Skip the log1p step.

    Doublet detection
    ~~~~~~~~~~~~~~~~~
    run_doublet_detection:
        Run Scrublet on raw counts and remove predicted doublets.
        Requires ``scrublet`` (``pip install scrublet``).
    doublet_threshold:
        Scrublet doublet-score threshold.  ``None`` = automatic threshold
        chosen by Scrublet.
    """

    def __init__(
        self,
        filter_strategy: FilterStrategy = "both",
        min_counts: int = 200,
        min_genes: int = 200,
        max_counts: int | None = None,
        max_genes: int | None = None,
        max_mito_pct: float | None = 20.0,
        auto_flag_mito: bool = True,
        target_sum: float = 1e4,
        log1p: bool = True,
        scale: bool = False,
        max_value: float | None = 10.0,
        layer_in: str | None = None,
        layer_out: str | None = None,
        # skip flags
        already_qc_filtered: bool = False,
        already_normalized: bool = False,
        already_log_transformed: bool = False,
        # doublet detection
        run_doublet_detection: bool = False,
        doublet_threshold: float | None = None,
    ):
        self.filter_strategy = filter_strategy
        self.min_counts = min_counts
        self.min_genes = min_genes
        self.max_counts = max_counts
        self.max_genes = max_genes
        self.max_mito_pct = max_mito_pct
        self.auto_flag_mito = auto_flag_mito
        self.target_sum = target_sum
        self.log1p = log1p
        self.scale = scale
        self.max_value = max_value
        self.layer_in = layer_in
        self.layer_out = layer_out
        self.already_qc_filtered = already_qc_filtered
        self.already_normalized = already_normalized
        self.already_log_transformed = already_log_transformed
        self.run_doublet_detection = run_doublet_detection
        self.doublet_threshold = doublet_threshold

    # ------------------------------------------------------------------
    def fit(self, adata: AnnData, y=None) -> SingleCellPreprocessor:
        """Learn gene-wise statistics (mean, std) if scaling is requested."""
        if self.scale:
            filtered = self._run_qc(adata)
            X = self._get_X(filtered)
            if not self.already_normalized:
                X = self._lib_normalise(X)
                if self.log1p and not self.already_log_transformed:
                    X = np.log1p(X)
            elif self.log1p and not self.already_log_transformed:
                X = np.log1p(X)
            self._scale_mean_ = X.mean(axis=0)
            _std = X.std(axis=0)
            self._scale_std_ = np.where(_std == 0, 1.0, _std)
        log.info("SingleCellPreprocessor fitted.")
        return self

    def transform(self, adata: AnnData, y=None) -> AnnData:
        """Apply the preprocessing pipeline, returning a new AnnData."""
        # Mito annotation (non-destructive, always useful)
        if self.auto_flag_mito or self.max_mito_pct is not None:
            adata = self._annotate_mito(adata)

        # QC filtering
        adata = self._run_qc(adata)

        # Doublet detection (on raw counts — before normalisation)
        if self.run_doublet_detection and not self.already_qc_filtered:
            adata = self._run_doublet_detection(adata)

        # Expression matrix
        X = self._get_X(adata)

        # Normalisation + log1p (skipped if data already processed)
        if not self.already_normalized:
            X = self._lib_normalise(X)
            if self.log1p and not self.already_log_transformed:
                X = np.log1p(X)
        elif self.log1p and not self.already_log_transformed:
            # Already normalised but not yet log-transformed
            X = np.log1p(X)

        # Scaling
        if self.scale:
            X = (X - self._scale_mean_) / self._scale_std_
            if self.max_value is not None:
                X = np.clip(X, -self.max_value, self.max_value)

        adata = adata.copy()
        self._set_X(adata, X.astype(np.float32))
        return adata

    # ------------------------------------------------------------------
    # QC helpers
    # ------------------------------------------------------------------

    def _annotate_mito(self, adata: AnnData) -> AnnData:
        """Add pct_mito column to adata.obs (non-destructive copy)."""
        X = self._get_X(adata)
        mito_mask = np.array([g.upper().startswith("MT-") for g in adata.var_names])
        if mito_mask.sum() == 0:
            log.debug("No MT- genes found in var_names — skipping mito annotation.")
            return adata
        total_counts = X.sum(axis=1)
        mito_counts = X[:, mito_mask].sum(axis=1)
        pct = np.where(total_counts > 0, mito_counts / total_counts * 100, 0.0)
        adata = adata.copy()
        adata.obs["pct_mito"] = pct
        return adata

    def _run_qc(self, adata: AnnData) -> AnnData:
        if self.already_qc_filtered or self.filter_strategy == "none":
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
        if self.max_genes is not None:
            mask &= gene_counts <= self.max_genes
        if self.max_mito_pct is not None and "pct_mito" in adata.obs.columns:
            mask &= adata.obs["pct_mito"].values <= self.max_mito_pct

        n_removed = int((~mask).sum())
        if n_removed > 0:
            log.info(
                "QC filter: removed %d / %d cells (%d retained).",
                n_removed,
                adata.n_obs,
                mask.sum(),
            )
        return adata[mask].copy()

    def _run_doublet_detection(self, adata: AnnData) -> AnnData:
        if not _HAS_SCRUBLET:
            log.warning(
                "scrublet not installed — doublet detection skipped "
                "(pip install scrublet)."
            )
            return adata

        X = self._get_X(adata)
        try:
            scrub = _scrublet.Scrublet(X)
            scores, predicted = scrub.scrub_doublets(threshold=self.doublet_threshold)
        except Exception as exc:
            log.warning("Scrublet failed (%s) — skipping doublet removal.", exc)
            return adata

        adata = adata.copy()
        adata.obs["doublet_score"] = scores
        adata.obs["predicted_doublet"] = predicted
        n_doublets = int(predicted.sum())
        if n_doublets > 0:
            log.info("Doublet detection: removing %d predicted doublets.", n_doublets)
            adata = adata[~predicted].copy()
        return adata

    # ------------------------------------------------------------------
    # Normalisation + helpers
    # ------------------------------------------------------------------

    def _lib_normalise(self, X: np.ndarray) -> np.ndarray:
        lib = X.sum(axis=1, keepdims=True)
        lib = np.where(lib == 0, 1.0, lib)
        return X / lib * self.target_sum

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
