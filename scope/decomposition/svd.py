"""Truncated SVD decomposition for scOPE.

This is the primary latent-factor method described in the scOPE manuscript.

Bulk phase
----------
  A_bulk = U_bulk Σ_bulk V^T
  Z_bulk = U_bulk Σ_bulk      (patient × latent factor embedding)

Single-cell projection
----------------------
  Z_sc = A'_sc · V            (cell × latent factor embedding)

  where A'_sc has been bulk-aligned (same gene-wise normalisation).
"""

from __future__ import annotations

from typing import Literal, Optional

import numpy as np
from anndata import AnnData
from sklearn.utils.extmath import randomized_svd

from scope.decomposition.base import BaseDecomposition
from scope.utils.logging import get_logger

log = get_logger(__name__)

SVDAlgorithm = Literal["randomized", "arpack", "full"]


class SVDDecomposition(BaseDecomposition):
    """Truncated SVD for bulk RNA-seq latent space learning.

    Computes the best rank-*k* approximation of the bulk expression matrix and
    stores gene loadings (V) for projecting new single-cell observations.

    Parameters
    ----------
    n_components:
        Number of singular vectors to retain.
    algorithm:
        SVD solver:

        * ``"randomized"`` — fast randomised SVD (sklearn); good default.
        * ``"arpack"``     — ARPACK partial SVD via scipy; more accurate for
          small *n_components*.
        * ``"full"``       — full dense SVD; only practical for small matrices.
    scale_by_singular_values:
        If ``True`` (default), Z = U · Σ (scale latent coords by singular
        values).  If ``False``, Z = U only (unit-variance factors).
    n_iter:
        Number of power-iteration passes for ``algorithm="randomized"``.
    random_state:
        Seed for reproducibility (used by randomised SVD).
    layer:
        AnnData layer to decompose. ``None`` → ``adata.X``.
    obsm_key:
        Key under which the embedding is stored in ``adata.obsm``.
    """

    obsm_key: str = "X_svd"

    def __init__(
        self,
        n_components: int = 50,
        algorithm: SVDAlgorithm = "randomized",
        scale_by_singular_values: bool = True,
        n_iter: int = 4,
        random_state: int = 42,
        layer: Optional[str] = None,
        obsm_key: str = "X_svd",
    ):
        super().__init__(n_components=n_components, layer=layer)
        self.algorithm = algorithm
        self.scale_by_singular_values = scale_by_singular_values
        self.n_iter = n_iter
        self.random_state = random_state
        self.obsm_key = obsm_key

    # ------------------------------------------------------------------
    def fit(self, adata: AnnData, y=None) -> "SVDDecomposition":
        """Decompose bulk expression matrix and store V and Σ.

        Parameters
        ----------
        adata:
            Preprocessed bulk AnnData (samples × genes).
        """
        X = self._get_X(adata)
        k = min(self.n_components, X.shape[0] - 1, X.shape[1] - 1)
        if k != self.n_components:
            log.warning(
                "Requested %d components but matrix allows at most %d; "
                "using %d.",
                self.n_components, k, k,
            )

        if self.algorithm == "randomized":
            U, S, Vt = randomized_svd(
                X,
                n_components=k,
                n_iter=self.n_iter,
                random_state=self.random_state,
            )
        elif self.algorithm == "arpack":
            from scipy.sparse.linalg import svds

            U, S, Vt = svds(X, k=k)
            # svds returns in ascending order — reverse for consistency
            idx = np.argsort(-S)
            U, S, Vt = U[:, idx], S[idx], Vt[idx, :]
        elif self.algorithm == "full":
            U_full, S_full, Vt_full = np.linalg.svd(X, full_matrices=False)
            U, S, Vt = U_full[:, :k], S_full[:k], Vt_full[:k, :]
        else:
            raise ValueError(f"Unknown SVD algorithm: {self.algorithm!r}")

        # Store gene loadings (V) — shape (n_genes, k)
        self.components_ = Vt          # (k, n_genes)
        self.singular_values_ = S      # (k,)
        self.n_components_ = k

        # Compute and store explained variance ratio
        total_var = np.sum(X ** 2)
        self.explained_variance_ = S ** 2 / X.shape[0]
        self._explained_variance_ratio_ = S ** 2 / total_var if total_var > 0 else np.zeros(k)

        log.info(
            "SVD fitted: %d components (cumulative EVR=%.3f).",
            k,
            self._explained_variance_ratio_.sum(),
        )
        return self

    def transform(self, adata: AnnData, y=None) -> AnnData:
        """Project *adata* into the SVD latent space.

        For bulk data this recomputes U·Σ from the stored V.
        For single-cell data the same projection A'·V^T is used.

        Parameters
        ----------
        adata:
            AnnData with the same genes (in the same order) as the bulk data
            used in ``fit``.

        Returns
        -------
        AnnData
            Copy with embedding stored in ``adata.obsm[self.obsm_key]``.
        """
        X = self._get_X(adata)
        # Z = X · V^T  where V^T = self.components_  (k × n_genes)
        # so Z = X · components_.T   →  (n_obs, k)
        Z = X @ self.components_.T

        if self.scale_by_singular_values:
            Z = Z * self.singular_values_[np.newaxis, :]

        adata = adata.copy()
        adata.obsm[self.obsm_key] = Z.astype(np.float32)
        return adata

    def explained_variance_ratio(self) -> Optional[np.ndarray]:
        """Fraction of total variance explained by each component."""
        return self._explained_variance_ratio_

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------

    def select_n_components_by_variance(self, threshold: float = 0.90) -> int:
        """Return the minimum number of components that explain *threshold* variance.

        Parameters
        ----------
        threshold:
            Cumulative explained-variance fraction, e.g. 0.90 for 90%.

        Returns
        -------
        int
        """
        cumvar = np.cumsum(self._explained_variance_ratio_)
        idx = np.searchsorted(cumvar, threshold)
        n = int(idx) + 1
        log.info(
            "%d components explain %.1f%% of variance.", n, threshold * 100
        )
        return n

    def scree_data(self) -> dict:
        """Return a dict suitable for plotting a scree plot."""
        return {
            "component": np.arange(1, self.n_components_ + 1),
            "singular_value": self.singular_values_,
            "explained_variance_ratio": self._explained_variance_ratio_,
            "cumulative_evr": np.cumsum(self._explained_variance_ratio_),
        }
