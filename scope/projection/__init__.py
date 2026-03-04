"""Projection utilities — thin re-exports from the decomposition layer."""

from scope.decomposition.nmf import NMFDecomposition
from scope.decomposition.svd import SVDDecomposition

__all__ = ["SVDDecomposition", "NMFDecomposition"]
