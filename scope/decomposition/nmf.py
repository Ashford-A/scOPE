"""Non-negative Matrix Factorization (NMF) decomposition for scOPE.

NMF is an alternative to SVD that produces non-negative, parts-based gene
programs — often more biologically interpretable as metagenes.

  A_bulk ≈ W_bulk · H        (W: samples × k,  H: k × genes)

Projection: Z_sc = A'_sc · H^+   (Moore-Penrose pseudoinverse of H)
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from anndata import AnnData
from sklearn.decomposition import NMF

from scope.decomposition.base import BaseDecomposition
from scope.utils.logging import get_logger

log = get_logger(__name__)

NMFInit = Literal["nndsvda", "nndsvdar", "random", "nndsvd"]
NMFSolver = Literal["cd", "mu"]
NMFBeta = Literal["frobenius", "kullback-leibler", "itakura-saito"]


class NMFDecomposition(BaseDecomposition):
    """NMF for bulk RNA-seq metagene learning.

    Parameters
    ----------
    n_components:
        Number of metagene programs.
    init:
        Initialisation strategy for NMF (see :class:`sklearn.decomposition.NMF`).
        ``"nndsvda"`` is recommended for sparse data.
    solver:
        ``"cd"`` (coordinate descent) or ``"mu"`` (multiplicative updates).
        Use ``"mu"`` for beta-divergence losses.
    beta_loss:
        Loss function: ``"frobenius"`` (Euclidean), ``"kullback-leibler"``, or
        ``"itakura-saito"``.
    max_iter:
        Maximum number of iterations.
    tol:
        Convergence tolerance.
    random_state:
        Random seed.
    l1_ratio:
        Regularisation mixing parameter (0 = L2 only, 1 = L1 only).
    alpha_W, alpha_H:
        Regularisation strengths for W and H respectively.
    shift_negative:
        If ``True``, shift the input matrix so all values ≥ 0 before fitting.
        Required when the input has been z-scored (negative values).
    layer:
        AnnData layer to use; ``None`` → ``adata.X``.
    obsm_key:
        Key for the embedding in ``adata.obsm``.
    """

    obsm_key: str = "X_nmf"

    def __init__(
        self,
        n_components: int = 30,
        init: NMFInit = "nndsvda",
        solver: NMFSolver = "cd",
        beta_loss: NMFBeta = "frobenius",
        max_iter: int = 500,
        tol: float = 1e-4,
        random_state: int = 42,
        l1_ratio: float = 0.0,
        alpha_W: float = 0.0,
        alpha_H: float = 0.0,
        shift_negative: bool = True,
        layer: str | None = None,
        obsm_key: str = "X_nmf",
    ):
        super().__init__(n_components=n_components, layer=layer)
        self.init = init
        self.solver = solver
        self.beta_loss = beta_loss
        self.max_iter = max_iter
        self.tol = tol
        self.random_state = random_state
        self.l1_ratio = l1_ratio
        self.alpha_W = alpha_W
        self.alpha_H = alpha_H
        self.shift_negative = shift_negative
        self.obsm_key = obsm_key

    def fit(self, adata: AnnData, y=None) -> NMFDecomposition:
        X = self._get_X(adata)
        X = self._ensure_nonneg(X)

        self._model = NMF(
            n_components=self.n_components,
            init=self.init,
            solver=self.solver,
            beta_loss=self.beta_loss,
            max_iter=self.max_iter,
            tol=self.tol,
            random_state=self.random_state,
            l1_ratio=self.l1_ratio,
            alpha_W=self.alpha_W,
            alpha_H=self.alpha_H,
        )
        self._model.fit(X)

        self.components_ = self._model.components_   # (k, n_genes)
        self.n_components_ = self.n_components
        # Pseudoinverse of H for projecting new data
        self._H_pinv_ = np.linalg.pinv(self.components_)  # (n_genes, k)
        log.info(
            "NMF fitted: %d components (reconstruction error=%.4f).",
            self.n_components,
            self._model.reconstruction_err_,
        )
        return self

    def transform(self, adata: AnnData, y=None) -> AnnData:
        X = self._get_X(adata)
        X = self._ensure_nonneg(X)
        Z = X @ self._H_pinv_   # (n_obs, k)
        Z = np.maximum(Z, 0.0)   # keep non-negative
        adata = adata.copy()
        adata.obsm[self.obsm_key] = Z.astype(np.float32)
        return adata

    def _ensure_nonneg(self, X: np.ndarray) -> np.ndarray:
        if self.shift_negative and X.min() < 0:
            self._shift_ = -X.min(axis=0)
            return X + self._shift_
        if hasattr(self, "_shift_"):
            return X + self._shift_
        return X

    @property
    def metagenes(self) -> np.ndarray:
        """Gene weight matrix H, shape (n_components, n_genes)."""
        return self.components_
