"""Factor Analysis decomposition for scOPE.

Factor Analysis is a probabilistic linear latent-variable model that
explicitly accounts for gene-specific noise variances (heteroscedasticity).
Unlike SVD/PCA which assume uniform noise, FA fits:

    X_i = W z_i + ε_i,   ε_i ~ N(0, Ψ)

where Ψ is a diagonal matrix of per-gene noise variances.  This makes FA
more appropriate for RNA-seq data where technical noise varies widely across
genes.

Projection onto new data (e.g. single-cell) uses the posterior mean:

    z_new = (W^T Ψ^{-1} W + I)^{-1} W^T Ψ^{-1} x_new
"""

from __future__ import annotations

import numpy as np
from anndata import AnnData
from sklearn.decomposition import FactorAnalysis

from scope.decomposition.base import BaseDecomposition
from scope.utils.logging import get_logger

log = get_logger(__name__)


class FactorAnalysisDecomposition(BaseDecomposition):
    """Factor Analysis for bulk RNA-seq latent space learning.

    Parameters
    ----------
    n_components:
        Number of latent factors.
    max_iter:
        Maximum number of EM iterations.
    tol:
        Convergence tolerance.
    svd_method:
        Algorithm used internally by sklearn FA: ``"lapack"`` or
        ``"randomized"``.
    random_state:
        Seed for reproducibility.
    layer:
        AnnData layer to decompose.  ``None`` → ``adata.X``.
    obsm_key:
        Key for the embedding in ``adata.obsm``.
    """

    obsm_key: str = "X_fa"

    def __init__(
        self,
        n_components: int = 30,
        max_iter: int = 1000,
        tol: float = 1e-3,
        svd_method: str = "randomized",
        random_state: int = 42,
        layer: str | None = None,
        obsm_key: str = "X_fa",
    ):
        super().__init__(n_components=n_components, layer=layer)
        self.max_iter = max_iter
        self.tol = tol
        self.svd_method = svd_method
        self.random_state = random_state
        self.obsm_key = obsm_key

    def fit(self, adata: AnnData, y=None) -> FactorAnalysisDecomposition:
        X = self._get_X(adata)
        self._model = FactorAnalysis(
            n_components=self.n_components,
            max_iter=self.max_iter,
            tol=self.tol,
            svd_method=self.svd_method,
            random_state=self.random_state,
        )
        self._model.fit(X)
        # components_: (n_components, n_genes) — factor loading matrix W^T
        self.components_ = self._model.components_
        self.noise_variance_ = self._model.noise_variance_  # (n_genes,)
        self.n_components_ = self.n_components
        log.info(
            "FactorAnalysis fitted: %d components.  " "Mean noise var = %.4f.",
            self.n_components,
            float(self.noise_variance_.mean()),
        )
        return self

    def transform(self, adata: AnnData, y=None) -> AnnData:
        X = self._get_X(adata)
        Z = self._model.transform(X)  # posterior mean
        adata = adata.copy()
        adata.obsm[self.obsm_key] = Z.astype(np.float32)
        return adata

    @property
    def loadings(self) -> np.ndarray:
        """Factor loading matrix W, shape (n_genes, n_components)."""
        return self.components_.T
