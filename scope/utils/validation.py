"""Input validation helpers used throughout scOPE.

Centralising checks here keeps the domain-logic modules clean and ensures
consistent, user-friendly error messages everywhere.
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData


def check_adata(
    adata: object,
    require_obs_names: bool = True,
    require_var_names: bool = True,
    min_obs: int = 2,
    min_vars: int = 10,
    name: str = "adata",
) -> AnnData:
    """Assert that *adata* is a valid :class:`anndata.AnnData`.

    Parameters
    ----------
    adata:
        Object to validate.
    require_obs_names:
        Raise if ``obs_names`` are not unique.
    require_var_names:
        Raise if ``var_names`` are not unique.
    min_obs:
        Minimum number of observations (samples / cells).
    min_vars:
        Minimum number of variables (genes).
    name:
        Variable name used in error messages.

    Returns
    -------
    AnnData
        The input object, unchanged.

    Raises
    ------
    TypeError
        If *adata* is not an :class:`~anndata.AnnData`.
    ValueError
        If any constraint is violated.
    """
    if not isinstance(adata, AnnData):
        raise TypeError(
            f"Expected '{name}' to be an AnnData object, got {type(adata).__name__}."
        )
    if adata.n_obs < min_obs:
        raise ValueError(
            f"'{name}' has only {adata.n_obs} observations; need at least {min_obs}."
        )
    if adata.n_vars < min_vars:
        raise ValueError(
            f"'{name}' has only {adata.n_vars} variables (genes); need at least {min_vars}."
        )
    if require_obs_names and adata.obs_names.duplicated().any():
        raise ValueError(f"'{name}.obs_names' contains duplicate identifiers.")
    if require_var_names and adata.var_names.duplicated().any():
        raise ValueError(f"'{name}.var_names' contains duplicate gene names.")
    return adata


def check_mutation_labels(
    labels: object,
    adata: Optional[AnnData] = None,
    min_mutations: int = 1,
    name: str = "mutation_labels",
) -> pd.DataFrame:
    """Assert that *labels* is a valid binary mutation label DataFrame.

    Parameters
    ----------
    labels:
        Object to validate.
    adata:
        If provided, check that ``labels.index`` and ``adata.obs_names``
        share at least one common identifier.
    min_mutations:
        Minimum number of mutation columns.
    name:
        Variable name used in error messages.

    Returns
    -------
    pd.DataFrame
        The input, unchanged.
    """
    if not isinstance(labels, pd.DataFrame):
        raise TypeError(
            f"Expected '{name}' to be a pd.DataFrame, got {type(labels).__name__}."
        )
    if labels.shape[1] < min_mutations:
        raise ValueError(
            f"'{name}' has {labels.shape[1]} mutation columns; need at least {min_mutations}."
        )
    # Check binary-ish (allow NaN, warn)
    unique_vals = set(labels.values.ravel())
    unique_vals.discard(np.nan)
    non_binary = unique_vals - {0, 1, 0.0, 1.0}
    if non_binary:
        import warnings
        warnings.warn(
            f"'{name}' contains non-binary values: {non_binary}. "
            "Expected 0/1 integer labels.",
            UserWarning,
            stacklevel=2,
        )
    if adata is not None:
        shared = set(labels.index) & set(adata.obs_names)
        if len(shared) == 0:
            raise ValueError(
                f"No overlapping identifiers between '{name}.index' and "
                "adata.obs_names. Check that sample IDs match."
            )
    return labels


def check_is_fitted(estimator: object, attributes: Sequence[str]) -> None:
    """Raise :class:`sklearn.exceptions.NotFittedError` if *estimator* is not fitted.

    Parameters
    ----------
    estimator:
        Estimator instance to check.
    attributes:
        Attribute names that should be present after fitting (trailing ``_``
        convention is checked automatically if the attribute list is given
        without the trailing underscore).

    Raises
    ------
    sklearn.exceptions.NotFittedError
    """
    from sklearn.exceptions import NotFittedError

    missing = [
        attr for attr in attributes if not hasattr(estimator, attr)
    ]
    if missing:
        cls_name = type(estimator).__name__
        raise NotFittedError(
            f"This {cls_name} instance is not fitted yet. "
            f"Missing attributes: {missing}. "
            f"Call 'fit' before using this estimator."
        )


def check_nonneg(
    X: np.ndarray,
    name: str = "X",
    raise_on_negative: bool = False,
) -> np.ndarray:
    """Warn (or raise) if *X* contains negative values.

    Parameters
    ----------
    X:
        Array to check.
    name:
        Variable name for the error/warning message.
    raise_on_negative:
        If ``True``, raise :class:`ValueError`; otherwise just warn.

    Returns
    -------
    np.ndarray
        The input array, unchanged.
    """
    import warnings

    if X.min() < 0:
        msg = (
            f"'{name}' contains {(X < 0).sum()} negative values. "
            "NMF requires non-negative input; consider setting "
            "shift_negative=True or using a different decomposition."
        )
        if raise_on_negative:
            raise ValueError(msg)
        warnings.warn(msg, UserWarning, stacklevel=2)
    return X


def check_gene_overlap(
    source_genes: Sequence[str],
    target_genes: Sequence[str],
    min_overlap: int = 100,
    warn_threshold: float = 0.5,
) -> None:
    """Validate the gene overlap between two gene lists.

    Parameters
    ----------
    source_genes:
        Gene list from the source dataset (e.g. sc).
    target_genes:
        Gene list from the target dataset (e.g. bulk).
    min_overlap:
        Raise if fewer than this many genes are shared.
    warn_threshold:
        Warn if the fraction of target genes covered falls below this value.
    """
    import warnings

    shared = set(source_genes) & set(target_genes)
    n_shared = len(shared)
    if n_shared < min_overlap:
        raise ValueError(
            f"Only {n_shared} genes shared between source and target gene sets "
            f"(minimum required: {min_overlap}). "
            "Ensure both datasets use the same gene symbol space."
        )
    coverage = n_shared / max(len(target_genes), 1)
    if coverage < warn_threshold:
        warnings.warn(
            f"Only {coverage:.1%} of target genes are covered by source data "
            f"({n_shared} / {len(target_genes)}). "
            "Low coverage may reduce projection quality.",
            UserWarning,
            stacklevel=2,
        )
