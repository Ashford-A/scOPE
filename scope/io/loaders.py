"""Data loaders for common single-cell and bulk RNA-seq formats."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from scope.utils.logging import get_logger

log = get_logger(__name__)

PathLike = Union[str, Path]


# ---------------------------------------------------------------------------
# Generic entry-point
# ---------------------------------------------------------------------------


def load(
    path: PathLike,
    fmt: Optional[str] = None,
    **kwargs,
) -> AnnData:
    """Load expression data from *path* and return an :class:`anndata.AnnData`.

    Format is inferred from the file extension when *fmt* is not provided.

    Parameters
    ----------
    path:
        Path to the data file or directory (for 10x MTX).
    fmt:
        One of ``"h5ad"``, ``"csv"``, ``"tsv"``, ``"10x_mtx"``, ``"10x_h5"``,
        ``"loom"``. Inferred from extension if ``None``.
    **kwargs:
        Passed verbatim to the underlying loader.

    Returns
    -------
    AnnData
    """
    path = Path(path)
    if fmt is None:
        fmt = _infer_format(path)
    loaders = {
        "h5ad": load_h5ad,
        "csv": load_delimited,
        "tsv": load_delimited,
        "10x_mtx": load_10x_mtx,
        "10x_h5": load_10x_h5,
        "loom": load_loom,
    }
    if fmt not in loaders:
        raise ValueError(f"Unsupported format '{fmt}'. Choose from: {list(loaders)}")
    log.info("Loading data from '%s' (format=%s).", path, fmt)
    return loaders[fmt](path, **kwargs)


def _infer_format(path: Path) -> str:
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith(".h5ad"):
        return "h5ad"
    if suffixes.endswith(".loom"):
        return "loom"
    if suffixes.endswith(".h5") or suffixes.endswith(".hdf5"):
        return "10x_h5"
    if suffixes.endswith(".csv"):
        return "csv"
    if suffixes.endswith(".tsv") or suffixes.endswith(".txt"):
        return "tsv"
    if path.is_dir():
        return "10x_mtx"
    raise ValueError(
        f"Cannot infer format from '{path}'. Please specify fmt= explicitly."
    )


# ---------------------------------------------------------------------------
# Specific loaders
# ---------------------------------------------------------------------------


def load_h5ad(path: PathLike, backed: Optional[str] = None) -> AnnData:
    """Load an ``h5ad`` file.

    Parameters
    ----------
    path:
        Path to ``.h5ad`` file.
    backed:
        If ``"r"`` or ``"r+"``, open in backed mode (memory-mapped).
    """
    import anndata as ad

    adata = ad.read_h5ad(str(path), backed=backed)
    log.info("Loaded h5ad: %d cells × %d genes.", adata.n_obs, adata.n_vars)
    return adata


def load_delimited(
    path: PathLike,
    sep: Optional[str] = None,
    genes_as_columns: bool = True,
    obs_label_col: Optional[str] = None,
) -> AnnData:
    """Load a CSV/TSV expression matrix.

    Parameters
    ----------
    path:
        Path to the delimited file.
    sep:
        Delimiter — inferred from extension if ``None``.
    genes_as_columns:
        If ``True`` (default), the matrix is (samples × genes). If ``False``,
        it is (genes × samples) and will be transposed.
    obs_label_col:
        Column name to use as observation (sample/cell) names. If ``None``,
        the row index of the CSV is used.
    """
    path = Path(path)
    if sep is None:
        sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    df = pd.read_csv(path, sep=sep, index_col=0)
    if not genes_as_columns:
        df = df.T
    if obs_label_col is not None and obs_label_col in df.columns:
        df.index = df.pop(obs_label_col)
    adata = AnnData(X=df.values.astype(np.float32))
    adata.obs_names = list(df.index.astype(str))
    adata.var_names = list(df.columns.astype(str))
    log.info(
        "Loaded delimited file: %d samples × %d genes.", adata.n_obs, adata.n_vars
    )
    return adata


def load_10x_mtx(
    path: PathLike,
    var_names: str = "gene_symbols",
    cache: bool = True,
) -> AnnData:
    """Load a 10x Genomics MTX directory (barcodes, features/genes, matrix).

    Parameters
    ----------
    path:
        Directory containing ``matrix.mtx(.gz)``, ``barcodes.tsv(.gz)``,
        and ``features.tsv(.gz)`` (or ``genes.tsv(.gz)``).
    var_names:
        ``"gene_symbols"`` or ``"gene_ids"``.
    cache:
        Cache the result for faster re-loading.
    """
    import scanpy as sc

    adata = sc.read_10x_mtx(str(path), var_names=var_names, cache=cache)
    log.info(
        "Loaded 10x MTX: %d cells × %d genes.", adata.n_obs, adata.n_vars
    )
    return adata


def load_10x_h5(path: PathLike, genome: Optional[str] = None) -> AnnData:
    """Load a 10x Genomics HDF5 file (``.h5``).

    Parameters
    ----------
    path:
        Path to ``.h5`` file.
    genome:
        Genome name for multi-genome h5 files.
    """
    import scanpy as sc

    adata = sc.read_10x_h5(str(path), genome=genome)
    log.info(
        "Loaded 10x H5: %d cells × %d genes.", adata.n_obs, adata.n_vars
    )
    return adata


def load_loom(path: PathLike, **kwargs) -> AnnData:
    """Load a Loom file.

    Parameters
    ----------
    path:
        Path to ``.loom`` file.
    **kwargs:
        Passed to :func:`anndata.read_loom`.
    """
    import anndata as ad

    adata = ad.read_loom(str(path), **kwargs)
    log.info("Loaded Loom: %d cells × %d genes.", adata.n_obs, adata.n_vars)
    return adata


# ---------------------------------------------------------------------------
# Mutation label helpers
# ---------------------------------------------------------------------------


def load_mutation_labels(
    path: PathLike,
    sample_col: str = "sample_id",
    sep: str = ",",
) -> pd.DataFrame:
    """Load a CSV with sample-level mutation labels.

    Expected format: rows = samples, columns = sample_id + one column per
    gene/mutation (0 = wild-type, 1 = mutant).

    Parameters
    ----------
    path:
        Path to mutation label file.
    sample_col:
        Column used as row index.
    sep:
        CSV delimiter.

    Returns
    -------
    pd.DataFrame
        Index = sample IDs, columns = mutation names, values ∈ {0, 1}.
    """
    df = pd.read_csv(str(path), sep=sep)
    if sample_col in df.columns:
        df = df.set_index(sample_col)
    log.info(
        "Loaded mutation labels: %d samples, %d mutations.", len(df), df.shape[1]
    )
    return df.astype(int)
