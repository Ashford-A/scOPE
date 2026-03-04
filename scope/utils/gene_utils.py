"""Gene-universe alignment and mapping utilities."""

from __future__ import annotations

import warnings
from collections.abc import Sequence

import numpy as np
from anndata import AnnData

from scope.utils.logging import get_logger

log = get_logger(__name__)


# ---------------------------------------------------------------------------
# Gene universe helpers
# ---------------------------------------------------------------------------


def get_shared_genes(
    adata_bulk: AnnData,
    adata_sc: AnnData,
    gene_key: str = "gene_names",
) -> list[str]:
    """Return the intersection of gene names present in both datasets.

    Parameters
    ----------
    adata_bulk:
        Bulk AnnData (``var_names`` are gene names).
    adata_sc:
        Single-cell AnnData (``var_names`` are gene names).
    gene_key:
        Ignored if ``var_names`` are already gene symbols; provided for
        datasets that store gene IDs in ``var[gene_key]``.

    Returns
    -------
    List[str]
        Sorted list of shared gene symbols.
    """
    bulk_genes = set(adata_bulk.var_names)
    sc_genes = set(adata_sc.var_names)
    shared = sorted(bulk_genes & sc_genes)
    n_bulk_only = len(bulk_genes - sc_genes)
    n_sc_only = len(sc_genes - bulk_genes)
    log.info(
        "Gene universe: %d shared, %d bulk-only, %d sc-only",
        len(shared),
        n_bulk_only,
        n_sc_only,
    )
    if len(shared) == 0:
        raise ValueError(
            "No shared genes found between bulk and single-cell datasets. "
            "Ensure var_names are in the same gene symbol/ID space."
        )
    if len(shared) < 500:
        warnings.warn(
            f"Only {len(shared)} genes shared — this may be too few for "
            "reliable decomposition.",
            UserWarning,
            stacklevel=2,
        )
    return shared


def subset_to_shared_genes(
    adata_bulk: AnnData,
    adata_sc: AnnData,
    gene_key: str = "gene_names",
    inplace: bool = False,
) -> tuple[AnnData, AnnData]:
    """Subset both AnnDatas to their shared gene universe.

    Parameters
    ----------
    adata_bulk, adata_sc:
        Input AnnDatas.
    gene_key:
        See :func:`get_shared_genes`.
    inplace:
        If ``True``, modify the objects in place; otherwise return copies.

    Returns
    -------
    Tuple[AnnData, AnnData]
        ``(adata_bulk_sub, adata_sc_sub)`` sharing the same ``var_names``.
    """
    shared = get_shared_genes(adata_bulk, adata_sc, gene_key=gene_key)
    if not inplace:
        adata_bulk = adata_bulk.copy()
        adata_sc = adata_sc.copy()
    adata_bulk = adata_bulk[:, shared].copy()
    adata_sc = adata_sc[:, shared].copy()
    return adata_bulk, adata_sc


def align_gene_order(
    matrix: np.ndarray,
    source_genes: Sequence[str],
    target_genes: Sequence[str],
    fill_value: float = 0.0,
) -> np.ndarray:
    """Reorder / pad *matrix* columns so they match *target_genes*.

    Parameters
    ----------
    matrix:
        Shape ``(n_samples, len(source_genes))``.
    source_genes:
        Gene list corresponding to *matrix* columns.
    target_genes:
        Desired column ordering.
    fill_value:
        Value used for genes present in *target_genes* but absent from
        *source_genes*.

    Returns
    -------
    np.ndarray
        Shape ``(n_samples, len(target_genes))``.
    """
    src_idx = {g: i for i, g in enumerate(source_genes)}
    out = np.full((matrix.shape[0], len(target_genes)), fill_value, dtype=matrix.dtype)
    for j, gene in enumerate(target_genes):
        if gene in src_idx:
            out[:, j] = matrix[:, src_idx[gene]]
    return out


def filter_variable_genes(
    adata: AnnData,
    n_top_genes: int = 2000,
    layer: str | None = None,
) -> list[str]:
    """Select highly variable genes by coefficient of variation.

    Parameters
    ----------
    adata:
        Input AnnData.
    n_top_genes:
        Number of top-variable genes to select.
    layer:
        Key in ``adata.layers`` to use. Falls back to ``adata.X`` if None.

    Returns
    -------
    List[str]
        ``var_names`` of selected genes.
    """
    X = adata.layers[layer] if layer is not None else adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    mean = np.mean(X, axis=0)
    std = np.std(X, axis=0)
    with np.errstate(invalid="ignore"):
        cv = np.where(mean > 0, std / mean, 0.0)
    top_idx = np.argsort(cv)[::-1][:n_top_genes]
    selected = list(np.array(adata.var_names)[top_idx])
    log.info("Selected %d highly variable genes (by CV).", len(selected))
    return selected
