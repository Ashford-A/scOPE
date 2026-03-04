"""Independent Component Analysis (ICA) decomposition for scOPE.

ICA identifies statistically independent latent sources, which can reveal
non-Gaussian expression programs that SVD/PCA might miss.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from anndata import AnnData
from sklearn.decomposition import FastICA

from scope.decomposition.base import BaseDecomposition
from scope.utils.logging import get_logger

log = get_logger(__name__)

ICAAlgorithm = Literal["parallel", "deflation"]
ICAFunction = Literal["logcosh", "exp", "cube"]


class ICADecomposition(BaseDecomposition):
    """FastICA decomposition for bulk RNA-seq.

    Parameters
    ----------
    n_components:
        Number of independent components.
    algorithm:
        ``"parallel"`` or ``"deflation"``.
    fun:
        Non-quadratic G-function: ``"logcosh"`` (default), ``"exp"``, or
        ``"cube"``.
    max_iter:
        Maximum iterations for ICA convergence.
    tol:
        Convergence threshold.
    whiten:
        Pre-whiten data before ICA (recommended).
    random_state:
        Seed for reproducibility.
    layer:
        AnnData layer; ``None`` → ``adata.X``.
    obsm_key:
        Key for the embedding in ``adata.obsm``.
    """

    obsm_key: str = "X_ica"

    def __init__(
        self,
        n_components: int = 30,
        algorithm: ICAAlgorithm = "parallel",
        fun: ICAFunction = "logcosh",
        max_iter: int = 500,
        tol: float = 1e-4,
        whiten: bool = True,
        random_state: int = 42,
        layer: str | None = None,
        obsm_key: str = "X_ica",
    ):
        super().__init__(n_components=n_components, layer=layer)
        self.algorithm = algorithm
        self.fun = fun
        self.max_iter = max_iter
        self.tol = tol
        self.whiten = whiten
        self.random_state = random_state
        self.obsm_key = obsm_key

    def fit(self, adata: AnnData, y=None) -> ICADecomposition:
        X = self._get_X(adata)
        self._model = FastICA(
            n_components=self.n_components,
            algorithm=self.algorithm,
            fun=self.fun,
            max_iter=self.max_iter,
            tol=self.tol,
            whiten=self.whiten,
            random_state=self.random_state,
        )
        self._model.fit(X)
        # components_: (n_components, n_features) — mixing direction
        self.components_ = self._model.components_
        self.n_components_ = self.n_components
        self._mixing_pinv_ = np.linalg.pinv(self._model.mixing_)
        log.info("ICA fitted: %d components.", self.n_components)
        return self

    def transform(self, adata: AnnData, y=None) -> AnnData:
        X = self._get_X(adata)
        Z = self._model.transform(X)
        adata = adata.copy()
        adata.obsm[self.obsm_key] = Z.astype(np.float32)
        return adata
