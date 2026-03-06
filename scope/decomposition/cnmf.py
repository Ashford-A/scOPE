"""Consensus NMF (cNMF) decomposition for scOPE.

cNMF (Kotliar et al., eLife 2019) addresses the instability of a single NMF
run by:

  1. Running NMF *n_iter* times with random initialisations.
  2. Collecting all component vectors (k × n_iter runs → n_iter·k rows).
  3. K-means clustering the combined set into *n_components* consensus
     programs.
  4. Fitting a final NMF run with the consensus H matrix fixed, so only the
     sample activations W are optimised.

The resulting gene programs are substantially more stable and reproducible
than single-run NMF, and are widely used in cancer transcriptomics to
identify robust expression programs.

Reference
---------
Kotliar, D. et al. (2019). Identifying gene expression programs of cell-type
identity and cellular activity with single-cell RNA-seq. eLife, 8, e43803.
"""

from __future__ import annotations

import warnings

import numpy as np
from anndata import AnnData
from sklearn.cluster import KMeans
from sklearn.decomposition import NMF

from scope.decomposition.base import BaseDecomposition
from scope.utils.logging import get_logger

log = get_logger(__name__)


class ConsensusNMFDecomposition(BaseDecomposition):
    """Consensus NMF for robust gene program discovery.

    Parameters
    ----------
    n_components:
        Number of consensus gene programs.
    n_iter:
        Number of independent NMF runs used to generate the component pool.
        More iterations → more stable programs; 50–100 is typical.
    max_iter_nmf:
        Maximum NMF iterations per run.
    tol:
        NMF convergence tolerance.
    random_state:
        Base seed; individual run seeds are derived from this.
    shift_negative:
        Shift input to non-negative range before NMF (required if data have
        been z-scored or centred).
    layer:
        AnnData layer to decompose.  ``None`` → ``adata.X``.
    obsm_key:
        Key for the embedding in ``adata.obsm``.
    """

    obsm_key: str = "X_cnmf"

    def __init__(
        self,
        n_components: int = 20,
        n_iter: int = 50,
        max_iter_nmf: int = 300,
        tol: float = 1e-4,
        random_state: int = 42,
        shift_negative: bool = True,
        layer: str | None = None,
        obsm_key: str = "X_cnmf",
    ):
        super().__init__(n_components=n_components, layer=layer)
        self.n_iter = n_iter
        self.max_iter_nmf = max_iter_nmf
        self.tol = tol
        self.random_state = random_state
        self.shift_negative = shift_negative
        self.obsm_key = obsm_key

    # ------------------------------------------------------------------
    def fit(self, adata: AnnData, y=None) -> ConsensusNMFDecomposition:
        X = self._get_X(adata)
        X = self._ensure_nonneg(X, fit=True)

        rng = np.random.default_rng(self.random_state)
        all_H: list[np.ndarray] = []

        log.info(
            "cNMF: running %d NMF iterations (k=%d)...", self.n_iter, self.n_components
        )
        for _i in range(self.n_iter):
            seed = int(rng.integers(0, 2**31))
            nmf = NMF(
                n_components=self.n_components,
                init="random",
                max_iter=self.max_iter_nmf,
                tol=self.tol,
                random_state=seed,
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                nmf.fit(X)
            all_H.append(nmf.components_)  # (k, n_genes)

        # Stack → (n_iter * k, n_genes) and cluster into k consensus programs
        H_stack = np.vstack(all_H)
        log.info(
            "cNMF: clustering %d component vectors → %d consensus programs.",
            H_stack.shape[0],
            self.n_components,
        )
        km = KMeans(
            n_clusters=self.n_components,
            random_state=self.random_state,
            n_init=10,
        )
        km.fit(H_stack)
        H_consensus = np.clip(km.cluster_centers_, 0, None)  # (k, n_genes)

        # Final NMF pass: fix H = consensus, optimise only W
        self._final_nmf = NMF(
            n_components=self.n_components,
            init="custom",
            max_iter=self.max_iter_nmf * 2,
            tol=self.tol,
            random_state=self.random_state,
        )
        W_init = rng.random((X.shape[0], self.n_components)) + 1e-6
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self._final_nmf.fit_transform(X, W=W_init, H=H_consensus)

        # Overwrite H with consensus (fit_transform may update H slightly)
        self._final_nmf.components_ = H_consensus
        self.components_ = H_consensus  # (k, n_genes)
        self.n_components_ = self.n_components

        log.info(
            "cNMF fitted: %d consensus programs from %d iterations.",
            self.n_components,
            self.n_iter,
        )
        return self

    def transform(self, adata: AnnData, y=None) -> AnnData:
        X = self._get_X(adata)
        X = self._ensure_nonneg(X, fit=False)
        Z = self._final_nmf.transform(X)
        Z = np.maximum(Z, 0.0)
        adata = adata.copy()
        adata.obsm[self.obsm_key] = Z.astype(np.float32)
        return adata

    # ------------------------------------------------------------------
    def _ensure_nonneg(self, X: np.ndarray, fit: bool) -> np.ndarray:
        if self.shift_negative:
            if fit:
                self._shift_ = float(max(0.0, -X.min()))
            if hasattr(self, "_shift_") and self._shift_ > 0:
                return X + self._shift_
        return X

    @property
    def metagenes(self) -> np.ndarray:
        """Consensus gene program matrix H, shape (n_components, n_genes)."""
        return self.components_
