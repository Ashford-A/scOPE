"""
scOPE: single-cell Oncological Prediction Explorer
===================================================
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("scope-bio")
except PackageNotFoundError:
    __version__ = "0.1.0-dev"

from scope.pipeline.bulk_pipeline import BulkPipeline
from scope.pipeline.sc_pipeline import SingleCellPipeline

__all__ = [
    "BulkPipeline",
    "SingleCellPipeline",
    "__version__",
]