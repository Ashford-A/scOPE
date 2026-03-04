"""PCA decomposition for scOPE (alternative / baseline to SVD)."""

from __future__ import annotations

import numpy as np
from anndata import AnnData
from sklearn.decomposition import PCA

from scope.decomposition.base import BaseDecomposition
from scope.utils.logging import get_logger

log = get_logger(__name__)


class PCADecomposition(BaseDecomposition):
    """PCA via scikit-learn, with AnnData integration.

    Closely related to :class:`~scope.decomposition.svd.SVDDecomposition`
    (they are mathematically equivalent on centred data) but provides
    additional PCA-specific diagnostics such as explained variance ratio
    and loading biplots.

    Parameters
    ----------
    n_components:
        Number of principal components.
    svd_solver:
        Solver passed to :class:`sklearn.decomposition.PCA`.
        ``"auto"`` lets sklearn decide.
    whiten:
        If ``True``, divide by singular values to give unit-variance factors.
    random_state:
        Seed for reproducibility (used for randomised solvers).
    layer:
        AnnData layer; ``None`` → ``adata.X``.
    obsm_key:
        Key for the embedding in ``adata.obsm``.
    """

    obsm_key: str = "X_pca"

    def __init__(
        self,
        n_components: int = 50,
        svd_solver: str = "auto",
        whiten: bool = False,
        random_state: int = 42,
        layer: str | None = None,
        obsm_key: str = "X_pca",
    ):
        super().__init__(n_components=n_components, layer=layer)
        self.svd_solver = svd_solver
        self.whiten = whiten
        self.random_state = random_state
        self.obsm_key = obsm_key

    def fit(self, adata: AnnData, y=None) -> PCADecomposition:
        X = self._get_X(adata)
        self._model = PCA(
            n_components=self.n_components,
            svd_solver=self.svd_solver,
            whiten=self.whiten,
            random_state=self.random_state,
        )
        self._model.fit(X)
        self.components_ = self._model.components_          # (k, n_genes)
        self._evr_ = self._model.explained_variance_ratio_
        self.n_components_ = self.n_components
        log.info(
            "PCA fitted: %d components (cumulative EVR=%.3f).",
            self.n_components,
            self._evr_.sum(),
        )
        return self

    def transform(self, adata: AnnData, y=None) -> AnnData:
        X = self._get_X(adata)
        Z = self._model.transform(X)
        adata = adata.copy()
        adata.obsm[self.obsm_key] = Z.astype(np.float32)
        return adata

    def explained_variance_ratio(self) -> np.ndarray:
        return self._evr_
