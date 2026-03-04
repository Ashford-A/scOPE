"""Base interfaces for per-mutation classifiers."""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.preprocessing import StandardScaler

from scope.utils.logging import get_logger

log = get_logger(__name__)


class BaseMutationClassifier(ABC, BaseEstimator, ClassifierMixin):
    """Interface every scOPE classifier must implement.

    Classifiers operate on the latent embedding Z (n_samples × n_components),
    not on raw gene expression.
    """

    @abstractmethod
    def fit(self, Z: np.ndarray, y: np.ndarray) -> BaseMutationClassifier:
        """Train on latent embedding *Z* with binary labels *y*."""

    @abstractmethod
    def predict_proba(self, Z: np.ndarray) -> np.ndarray:
        """Return (n_samples, 2) probability array."""

    def predict(self, Z: np.ndarray) -> np.ndarray:
        return (self.predict_proba(Z)[:, 1] >= 0.5).astype(int)

    def score(self, Z: np.ndarray, y: np.ndarray) -> float:
        from sklearn.metrics import roc_auc_score
        proba = self.predict_proba(Z)[:, 1]
        return roc_auc_score(y, proba)


class PerMutationClassifierSet:
    """Trains and stores one classifier per mutation target.

    Parameters
    ----------
    classifier_factory:
        A callable that returns a fresh, unfitted :class:`BaseMutationClassifier`
        each time it is called.  This is typically a lambda or partial:
        ``lambda: LogisticMutationClassifier(C=0.1)``.
    min_positive_frac:
        Skip training for any mutation where fewer than this fraction of
        samples are mutant (avoids degenerate single-class problems).
    scale_features:
        If ``True``, prepend a ``StandardScaler`` to the classifier.
    """

    def __init__(
        self,
        classifier_factory,
        min_positive_frac: float = 0.05,
        scale_features: bool = True,
    ):
        self.classifier_factory = classifier_factory
        self.min_positive_frac = min_positive_frac
        self.scale_features = scale_features
        self.classifiers_: dict[str, BaseMutationClassifier] = {}
        self.skipped_: list[str] = []

    def fit(self, Z: np.ndarray, labels: pd.DataFrame) -> PerMutationClassifierSet:
        """Train one classifier per column of *labels*.

        Parameters
        ----------
        Z:
            Latent embedding, shape (n_samples, n_components).
        labels:
            Binary mutation label DataFrame; rows aligned with Z.
        """
        for mutation in labels.columns:
            y = labels[mutation].values
            pos_frac = y.mean()
            if pos_frac < self.min_positive_frac or pos_frac > (1 - self.min_positive_frac):
                log.warning(
                    "Skipping '%s': positive fraction=%.3f (threshold=%.3f).",
                    mutation, pos_frac, self.min_positive_frac,
                )
                self.skipped_.append(mutation)
                continue

            clf = self.classifier_factory()
            if self.scale_features:
                clf = _SklearnWrapper(clf)
            clf.fit(Z, y)
            self.classifiers_[mutation] = clf
            log.info(
                "Trained classifier for '%s' (n_pos=%d / n_tot=%d).",
                mutation, y.sum(), len(y),
            )
        return self

    def predict_proba(self, Z: np.ndarray) -> pd.DataFrame:
        """Return per-mutation probability DataFrame, shape (n_cells, n_mutations)."""
        result = {}
        for mutation, clf in self.classifiers_.items():
            result[f"mutation_prob_{mutation}"] = clf.predict_proba(Z)[:, 1]
        return pd.DataFrame(result)

    def cross_validate(
        self,
        Z: np.ndarray,
        labels: pd.DataFrame,
        cv: int = 5,
        scoring: str = "roc_auc",
    ) -> dict[str, dict]:
        """Run stratified k-fold CV for each mutation.

        Returns
        -------
        Dict[str, dict]
            ``{mutation: {"test_score": array, "mean": float, "std": float}}``
        """
        results = {}
        skf = StratifiedKFold(n_splits=cv, shuffle=True, random_state=42)
        for mutation in labels.columns:
            y = labels[mutation].values
            if y.mean() < self.min_positive_frac:
                continue
            clf = self.classifier_factory()
            cv_res = cross_validate(
                clf, Z, y, cv=skf, scoring=scoring, return_train_score=False
            )
            results[mutation] = {
                "test_score": cv_res["test_score"],
                "mean": cv_res["test_score"].mean(),
                "std": cv_res["test_score"].std(),
            }
            log.info(
                "CV '%s': %s mean=%.3f ± %.3f",
                mutation,
                scoring,
                results[mutation]["mean"],
                results[mutation]["std"],
            )
        return results


class _SklearnWrapper(BaseMutationClassifier):
    """Wraps any BaseMutationClassifier with a StandardScaler prefix."""

    def __init__(self, clf: BaseMutationClassifier):
        self._clf = clf
        self._scaler = StandardScaler()

    def fit(self, Z, y):
        Z_s = self._scaler.fit_transform(Z)
        self._clf.fit(Z_s, y)
        return self

    def predict_proba(self, Z):
        Z_s = self._scaler.transform(Z)
        return self._clf.predict_proba(Z_s)
