"""Hyperparameter search utilities for scOPE.

Provides convenience wrappers for sweeping decomposition and classifier
parameters while tracking performance metrics, without requiring the user
to write boilerplate loops.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable
from typing import Any

import pandas as pd
from anndata import AnnData

from scope.utils.logging import get_logger

log = get_logger(__name__)


# ---------------------------------------------------------------------------
# Parameter grid helpers
# ---------------------------------------------------------------------------


def _product_grid(param_grid: dict[str, list[Any]]) -> list[dict[str, Any]]:
    """Expand a parameter grid dict into a list of individual param dicts."""
    keys = list(param_grid.keys())
    values = list(param_grid.values())
    return [dict(zip(keys, combo)) for combo in itertools.product(*values)]


# ---------------------------------------------------------------------------
# Decomposition sweep
# ---------------------------------------------------------------------------


def sweep_n_components(
    adata_bulk: AnnData,
    mutation_labels: pd.DataFrame,
    n_components_list: Iterable[int],
    decomposition: str = "svd",
    classifier: str = "logistic",
    classifier_kwargs: dict | None = None,
    cv: int = 5,
    norm_method: str = "cpm",
    layer: str | None = None,
) -> pd.DataFrame:
    """Sweep over latent dimensionalities and report CV performance.

    Parameters
    ----------
    adata_bulk:
        Raw bulk AnnData.
    mutation_labels:
        Binary mutation label DataFrame.
    n_components_list:
        Sequence of component counts to try (e.g. ``[10, 25, 50, 100]``).
    decomposition:
        Decomposition method name.
    classifier:
        Classifier name.
    classifier_kwargs:
        Extra kwargs for the classifier.
    cv:
        Number of cross-validation folds.
    norm_method:
        Bulk normalisation method.
    layer:
        AnnData layer to use.

    Returns
    -------
    pd.DataFrame
        Columns: ``n_components``, ``mutation``, ``fold``, ``auroc``, ``auprc``.
    """
    from scope.pipeline.bulk_pipeline import BulkPipeline

    results = []
    for k in n_components_list:
        log.info("Sweeping n_components=%d …", k)
        pipe = BulkPipeline(
            norm_method=norm_method,
            decomposition=decomposition,
            n_components=k,
            classifier=classifier,
            classifier_kwargs=classifier_kwargs or {},
            layer=layer,
        )
        pipe.fit(adata_bulk, mutation_labels, cv=cv)
        if pipe.cv_results_ is not None:
            df = pipe.cv_results_.copy()
            df["n_components"] = k
            results.append(df)
    if not results:
        return pd.DataFrame()
    return pd.concat(results, ignore_index=True)


# ---------------------------------------------------------------------------
# Full grid search over pipeline parameters
# ---------------------------------------------------------------------------


def grid_search_pipeline(
    adata_bulk: AnnData,
    mutation_labels: pd.DataFrame,
    param_grid: dict[str, list[Any]],
    cv: int = 5,
    scoring_mutation: str | None = None,
    layer: str | None = None,
) -> pd.DataFrame:
    """Exhaustive grid search over scOPE ``BulkPipeline`` hyperparameters.

    Parameters
    ----------
    adata_bulk:
        Raw bulk AnnData.
    mutation_labels:
        Binary mutation label DataFrame.
    param_grid:
        Dict mapping ``BulkPipeline`` init parameter names to lists of values
        to try.  Example::

            {
                "decomposition": ["svd", "nmf"],
                "n_components": [25, 50],
                "classifier": ["logistic", "random_forest"],
            }

    cv:
        Cross-validation folds.
    scoring_mutation:
        If provided, restrict the summary to results for this mutation only.
        If ``None``, all mutations are included.
    layer:
        AnnData layer.

    Returns
    -------
    pd.DataFrame
        One row per (param combination × mutation × CV fold) with columns
        from the CV results plus all hyperparameter columns.
    """
    from scope.pipeline.bulk_pipeline import BulkPipeline

    combos = _product_grid(param_grid)
    log.info("Grid search: %d parameter combinations × %d folds.", len(combos), cv)
    all_results = []
    for i, params in enumerate(combos):
        log.info("[%d/%d] %s", i + 1, len(combos), params)
        try:
            pipe = BulkPipeline(layer=layer, **params)
            pipe.fit(adata_bulk, mutation_labels, cv=cv)
            if pipe.cv_results_ is not None:
                df = pipe.cv_results_.copy()
                for k, v in params.items():
                    df[k] = v
                all_results.append(df)
        except Exception as exc:
            log.warning("Combination %s failed: %s", params, exc)
    if not all_results:
        return pd.DataFrame()
    combined = pd.concat(all_results, ignore_index=True)
    if scoring_mutation is not None:
        combined = combined[combined["mutation"] == scoring_mutation]
    return combined


def summarise_grid_search(
    results: pd.DataFrame,
    param_cols: list[str] | None = None,
    metric: str = "auroc",
) -> pd.DataFrame:
    """Aggregate grid-search results to mean ± std per parameter combination.

    Parameters
    ----------
    results:
        Raw output of :func:`grid_search_pipeline`.
    param_cols:
        Columns that identify the hyperparameter combination.  If ``None``,
        auto-detected as any non-metric, non-index column.
    metric:
        Metric to aggregate (``"auroc"`` or ``"auprc"``).

    Returns
    -------
    pd.DataFrame
        Sorted descending by mean metric.
    """
    fixed_cols = {"mutation", "fold", "auroc", "auprc", "brier", "n_pos", "n_neg"}
    if param_cols is None:
        param_cols = [c for c in results.columns if c not in fixed_cols]
    group_cols = param_cols + ["mutation"]
    agg = (
        results.groupby(group_cols)[metric]
        .agg(mean="mean", std="std", n_folds="count")
        .reset_index()
        .sort_values("mean", ascending=False)
    )
    return agg
