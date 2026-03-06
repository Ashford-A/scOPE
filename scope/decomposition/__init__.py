"""Latent-factor decomposition methods."""

from scope.decomposition.base import BaseDecomposition
from scope.decomposition.cnmf import ConsensusNMFDecomposition
from scope.decomposition.factor_analysis import FactorAnalysisDecomposition
from scope.decomposition.ica import ICADecomposition
from scope.decomposition.nmf import NMFDecomposition
from scope.decomposition.pca import PCADecomposition
from scope.decomposition.svd import SVDDecomposition

__all__ = [
    "BaseDecomposition",
    "SVDDecomposition",
    "NMFDecomposition",
    "ICADecomposition",
    "PCADecomposition",
    "FactorAnalysisDecomposition",
    "ConsensusNMFDecomposition",
    "get_decomposition",
]

_REGISTRY = {
    "svd": SVDDecomposition,
    "nmf": NMFDecomposition,
    "ica": ICADecomposition,
    "pca": PCADecomposition,
    "fa": FactorAnalysisDecomposition,
    "cnmf": ConsensusNMFDecomposition,
}


def get_decomposition(method: str, **kwargs) -> BaseDecomposition:
    """Factory: instantiate a decomposition by name.

    Available methods: ``"svd"``, ``"nmf"``, ``"ica"``, ``"pca"``,
    ``"fa"`` (Factor Analysis), ``"cnmf"`` (Consensus NMF).
    """
    method = method.lower()
    if method not in _REGISTRY:
        raise ValueError(
            f"Unknown decomposition '{method}'. " f"Choose from: {list(_REGISTRY)}"
        )
    return _REGISTRY[method](**kwargs)
