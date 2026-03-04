"""Evaluation metrics and cross-validation utilities for scOPE classifiers."""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold

from scope.utils.logging import get_logger

log = get_logger(__name__)


# ---------------------------------------------------------------------------
# Per-mutation evaluation summary
# ---------------------------------------------------------------------------


def evaluate_classifier(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    mutation_name: str = "mutation",
) -> Dict[str, float]:
    """Compute a comprehensive set of binary classification metrics.

    Parameters
    ----------
    y_true:
        Ground-truth binary labels (0/1).
    y_prob:
        Predicted probability of the positive class.
    mutation_name:
        Name tag for logging.

    Returns
    -------
    dict with keys: ``auroc``, ``auprc``, ``brier``, ``n_pos``, ``n_neg``.
    """
    n_pos = int(y_true.sum())
    n_neg = int((1 - y_true).sum())
    if n_pos == 0 or n_neg == 0:
        log.warning("'%s': degenerate label vector — skipping metrics.", mutation_name)
        return {"auroc": float("nan"), "auprc": float("nan"), "brier": float("nan"),
                "n_pos": n_pos, "n_neg": n_neg}
    auroc = roc_auc_score(y_true, y_prob)
    auprc = average_precision_score(y_true, y_prob)
    brier = brier_score_loss(y_true, y_prob)
    log.info(
        "  %-20s  AUROC=%.3f  AUPRC=%.3f  Brier=%.3f  (pos=%d / %d)",
        mutation_name, auroc, auprc, brier, n_pos, n_pos + n_neg,
    )
    return {"auroc": auroc, "auprc": auprc, "brier": brier, "n_pos": n_pos, "n_neg": n_neg}


def evaluate_all(
    y_true_df: pd.DataFrame,
    y_prob_df: pd.DataFrame,
) -> pd.DataFrame:
    """Evaluate all mutations in a label matrix.

    Parameters
    ----------
    y_true_df:
        Binary label DataFrame (rows = samples, columns = mutations).
    y_prob_df:
        Predicted-probability DataFrame with matching column names
        (``mutation_prob_<name>`` → ``<name>`` mapping is handled).

    Returns
    -------
    pd.DataFrame
        Metrics indexed by mutation name.
    """
    rows = []
    for col in y_true_df.columns:
        prob_col = f"mutation_prob_{col}" if f"mutation_prob_{col}" in y_prob_df.columns else col
        if prob_col not in y_prob_df.columns:
            log.warning("No probability column for '%s' — skipping.", col)
            continue
        metrics = evaluate_classifier(y_true_df[col].values, y_prob_df[prob_col].values, col)
        metrics["mutation"] = col
        rows.append(metrics)
    return pd.DataFrame(rows).set_index("mutation")


def cross_validate_classifiers(
    Z: np.ndarray,
    labels: pd.DataFrame,
    classifier_factory,
    cv: int = 5,
    random_state: int = 42,
) -> pd.DataFrame:
    """Stratified k-fold cross-validation for all mutations.

    Parameters
    ----------
    Z:
        Latent embedding (n_samples × n_components).
    labels:
        Binary label DataFrame.
    classifier_factory:
        Callable returning a fresh unfitted classifier.
    cv:
        Number of folds.
    random_state:
        Seed for fold splitting.

    Returns
    -------
    pd.DataFrame
        Columns: ``mutation``, ``fold``, ``auroc``, ``auprc``.
    """
    skf = StratifiedKFold(n_splits=cv, shuffle=True, random_state=random_state)
    rows = []
    for mutation in labels.columns:
        y = labels[mutation].values
        if y.mean() < 0.05 or y.mean() > 0.95:
            continue
        for fold_idx, (train_idx, test_idx) in enumerate(skf.split(Z, y)):
            clf = classifier_factory()
            clf.fit(Z[train_idx], y[train_idx])
            y_prob = clf.predict_proba(Z[test_idx])[:, 1]
            metrics = evaluate_classifier(y[test_idx], y_prob, mutation_name=f"{mutation}/fold{fold_idx}")
            rows.append({"mutation": mutation, "fold": fold_idx, **metrics})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Curve data (for plotting)
# ---------------------------------------------------------------------------


def roc_curve_data(
    y_true: np.ndarray,
    y_prob: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Return (fpr, tpr, auroc) for a ROC curve."""
    fpr, tpr, _ = roc_curve(y_true, y_prob)
    auroc = roc_auc_score(y_true, y_prob)
    return fpr, tpr, auroc


def pr_curve_data(
    y_true: np.ndarray,
    y_prob: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Return (precision, recall, auprc) for a PR curve."""
    precision, recall, _ = precision_recall_curve(y_true, y_prob)
    auprc = average_precision_score(y_true, y_prob)
    return precision, recall, auprc
