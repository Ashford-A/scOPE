"""Concrete mutation-prediction classifier implementations.

All classifiers wrap scikit-learn (or compatible) estimators and expose the
:class:`~scope.classification.base.BaseMutationClassifier` interface so they
can be used interchangeably in :class:`~scope.classification.base.PerMutationClassifierSet`.
"""

from __future__ import annotations

from typing import Literal, Optional

import numpy as np
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import (
    GradientBoostingClassifier,
    RandomForestClassifier,
)
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC

from scope.classification.base import BaseMutationClassifier
from scope.utils.logging import get_logger

log = get_logger(__name__)


# ---------------------------------------------------------------------------
# Logistic Regression
# ---------------------------------------------------------------------------


class LogisticMutationClassifier(BaseMutationClassifier):
    """Regularised logistic regression on the latent embedding.

    Parameters
    ----------
    C:
        Inverse regularisation strength.  Smaller = stronger regularisation.
    penalty:
        ``"l1"``, ``"l2"``, or ``"elasticnet"``.
    solver:
        Solver (must be compatible with *penalty*).
    max_iter:
        Maximum solver iterations.
    class_weight:
        ``"balanced"`` recommended for imbalanced mutation labels.
    random_state:
        Seed.
    """

    def __init__(
        self,
        C: float = 1.0,
        penalty: Literal["l1", "l2", "elasticnet"] = "l2",
        solver: str = "lbfgs",
        max_iter: int = 1000,
        class_weight: Optional[str] = "balanced",
        random_state: int = 42,
    ):
        self.C = C
        self.penalty = penalty
        self.solver = solver
        self.max_iter = max_iter
        self.class_weight = class_weight
        self.random_state = random_state

    def fit(self, Z: np.ndarray, y: np.ndarray) -> "LogisticMutationClassifier":
        self._model = LogisticRegression(
            C=self.C,
            penalty=self.penalty,
            solver=self.solver,
            max_iter=self.max_iter,
            class_weight=self.class_weight,
            random_state=self.random_state,
        )
        self._model.fit(Z, y)
        return self

    def predict_proba(self, Z: np.ndarray) -> np.ndarray:
        return self._model.predict_proba(Z)

    @property
    def coef_(self) -> np.ndarray:
        """Logistic regression coefficients (1 × n_components)."""
        return self._model.coef_


# ---------------------------------------------------------------------------
# Random Forest
# ---------------------------------------------------------------------------


class RandomForestMutationClassifier(BaseMutationClassifier):
    """Random forest on the latent embedding.

    Parameters
    ----------
    n_estimators:
        Number of trees.
    max_depth:
        Maximum tree depth. ``None`` → grow until pure.
    min_samples_leaf:
        Minimum samples per leaf.
    class_weight:
        ``"balanced"`` or ``"balanced_subsample"`` recommended.
    n_jobs:
        Parallelism (−1 = all cores).
    random_state:
        Seed.
    """

    def __init__(
        self,
        n_estimators: int = 300,
        max_depth: Optional[int] = None,
        min_samples_leaf: int = 1,
        class_weight: Optional[str] = "balanced",
        n_jobs: int = -1,
        random_state: int = 42,
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.min_samples_leaf = min_samples_leaf
        self.class_weight = class_weight
        self.n_jobs = n_jobs
        self.random_state = random_state

    def fit(self, Z: np.ndarray, y: np.ndarray) -> "RandomForestMutationClassifier":
        self._model = RandomForestClassifier(
            n_estimators=self.n_estimators,
            max_depth=self.max_depth,
            min_samples_leaf=self.min_samples_leaf,
            class_weight=self.class_weight,
            n_jobs=self.n_jobs,
            random_state=self.random_state,
        )
        self._model.fit(Z, y)
        return self

    def predict_proba(self, Z: np.ndarray) -> np.ndarray:
        return self._model.predict_proba(Z)

    @property
    def feature_importances_(self) -> np.ndarray:
        return self._model.feature_importances_


# ---------------------------------------------------------------------------
# Gradient Boosting (sklearn)
# ---------------------------------------------------------------------------


class GBMMutationClassifier(BaseMutationClassifier):
    """Gradient boosting machine on the latent embedding.

    Parameters
    ----------
    n_estimators:
        Number of boosting rounds.
    learning_rate:
        Step-size shrinkage.
    max_depth:
        Depth of individual trees.
    subsample:
        Fraction of samples used per tree (stochastic GBM).
    random_state:
        Seed.
    """

    def __init__(
        self,
        n_estimators: int = 200,
        learning_rate: float = 0.05,
        max_depth: int = 3,
        subsample: float = 0.8,
        random_state: int = 42,
    ):
        self.n_estimators = n_estimators
        self.learning_rate = learning_rate
        self.max_depth = max_depth
        self.subsample = subsample
        self.random_state = random_state

    def fit(self, Z: np.ndarray, y: np.ndarray) -> "GBMMutationClassifier":
        self._model = GradientBoostingClassifier(
            n_estimators=self.n_estimators,
            learning_rate=self.learning_rate,
            max_depth=self.max_depth,
            subsample=self.subsample,
            random_state=self.random_state,
        )
        self._model.fit(Z, y)
        return self

    def predict_proba(self, Z: np.ndarray) -> np.ndarray:
        return self._model.predict_proba(Z)


# ---------------------------------------------------------------------------
# XGBoost (optional dependency)
# ---------------------------------------------------------------------------


class XGBMutationClassifier(BaseMutationClassifier):
    """XGBoost classifier on the latent embedding.

    Requires the ``xgboost`` package (install with ``pip install scope-bio[full]``).

    Parameters
    ----------
    n_estimators:
        Boosting rounds.
    max_depth:
        Tree depth.
    learning_rate:
        Step size.
    subsample:
        Row subsampling ratio.
    colsample_bytree:
        Feature subsampling per tree.
    scale_pos_weight:
        Weight applied to the positive class. Useful for imbalanced labels.
        Set to ``None`` to compute automatically from y.
    n_jobs:
        Parallelism.
    random_state:
        Seed.
    """

    def __init__(
        self,
        n_estimators: int = 300,
        max_depth: int = 4,
        learning_rate: float = 0.05,
        subsample: float = 0.8,
        colsample_bytree: float = 0.8,
        scale_pos_weight: Optional[float] = None,
        n_jobs: int = -1,
        random_state: int = 42,
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.learning_rate = learning_rate
        self.subsample = subsample
        self.colsample_bytree = colsample_bytree
        self.scale_pos_weight = scale_pos_weight
        self.n_jobs = n_jobs
        self.random_state = random_state

    def fit(self, Z: np.ndarray, y: np.ndarray) -> "XGBMutationClassifier":
        try:
            from xgboost import XGBClassifier
        except ImportError as e:
            raise ImportError(
                "XGBoost is required for XGBMutationClassifier. "
                "Install it with: pip install scope-bio[full]"
            ) from e
        spw = self.scale_pos_weight
        if spw is None:
            n_neg = (y == 0).sum()
            n_pos = (y == 1).sum()
            spw = n_neg / n_pos if n_pos > 0 else 1.0
        self._model = XGBClassifier(
            n_estimators=self.n_estimators,
            max_depth=self.max_depth,
            learning_rate=self.learning_rate,
            subsample=self.subsample,
            colsample_bytree=self.colsample_bytree,
            scale_pos_weight=spw,
            n_jobs=self.n_jobs,
            random_state=self.random_state,
            eval_metric="logloss",
            use_label_encoder=False,
        )
        self._model.fit(Z, y)
        return self

    def predict_proba(self, Z: np.ndarray) -> np.ndarray:
        return self._model.predict_proba(Z)


# ---------------------------------------------------------------------------
# LightGBM (optional dependency)
# ---------------------------------------------------------------------------


class LGBMMutationClassifier(BaseMutationClassifier):
    """LightGBM classifier on the latent embedding.

    Requires ``lightgbm`` (``pip install scope-bio[full]``).
    """

    def __init__(
        self,
        n_estimators: int = 300,
        max_depth: int = -1,
        learning_rate: float = 0.05,
        num_leaves: int = 31,
        class_weight: Optional[str] = "balanced",
        n_jobs: int = -1,
        random_state: int = 42,
        verbose: int = -1,
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.learning_rate = learning_rate
        self.num_leaves = num_leaves
        self.class_weight = class_weight
        self.n_jobs = n_jobs
        self.random_state = random_state
        self.verbose = verbose

    def fit(self, Z: np.ndarray, y: np.ndarray) -> "LGBMMutationClassifier":
        try:
            from lightgbm import LGBMClassifier
        except ImportError as e:
            raise ImportError(
                "LightGBM is required for LGBMMutationClassifier. "
                "Install it with: pip install scope-bio[full]"
            ) from e
        self._model = LGBMClassifier(
            n_estimators=self.n_estimators,
            max_depth=self.max_depth,
            learning_rate=self.learning_rate,
            num_leaves=self.num_leaves,
            class_weight=self.class_weight,
            n_jobs=self.n_jobs,
            random_state=self.random_state,
            verbose=self.verbose,
        )
        self._model.fit(Z, y)
        return self

    def predict_proba(self, Z: np.ndarray) -> np.ndarray:
        return self._model.predict_proba(Z)


# ---------------------------------------------------------------------------
# SVM
# ---------------------------------------------------------------------------


class SVMMutationClassifier(BaseMutationClassifier):
    """Support vector machine (RBF kernel) with Platt calibration.

    Probabilities are obtained via Platt scaling (CalibratedClassifierCV).

    Parameters
    ----------
    C:
        Regularisation.
    kernel:
        Kernel type.
    gamma:
        Kernel coefficient.  ``"scale"`` = 1 / (n_features * X.var()).
    class_weight:
        ``"balanced"`` recommended.
    cv:
        Cross-validation folds for Platt calibration.
    """

    def __init__(
        self,
        C: float = 1.0,
        kernel: str = "rbf",
        gamma: str = "scale",
        class_weight: Optional[str] = "balanced",
        cv: int = 5,
        random_state: int = 42,
    ):
        self.C = C
        self.kernel = kernel
        self.gamma = gamma
        self.class_weight = class_weight
        self.cv = cv
        self.random_state = random_state

    def fit(self, Z: np.ndarray, y: np.ndarray) -> "SVMMutationClassifier":
        base = SVC(
            C=self.C,
            kernel=self.kernel,
            gamma=self.gamma,
            class_weight=self.class_weight,
            probability=False,  # CalibratedClassifierCV handles this
            random_state=self.random_state,
        )
        self._model = CalibratedClassifierCV(base, cv=self.cv, method="sigmoid")
        self._model.fit(Z, y)
        return self

    def predict_proba(self, Z: np.ndarray) -> np.ndarray:
        return self._model.predict_proba(Z)


# ---------------------------------------------------------------------------
# MLP
# ---------------------------------------------------------------------------


class MLPMutationClassifier(BaseMutationClassifier):
    """Multi-layer perceptron classifier.

    Parameters
    ----------
    hidden_layer_sizes:
        Tuple of hidden-layer widths.
    activation:
        Activation function: ``"relu"`` (default), ``"tanh"``, ``"logistic"``.
    dropout:
        Note: sklearn's MLP does not support dropout natively.  If you need
        regularised MLPs set ``alpha`` (L2) instead.
    alpha:
        L2 regularisation coefficient.
    learning_rate_init:
        Initial learning rate.
    max_iter:
        Maximum training epochs.
    early_stopping:
        Use a validation fraction for early stopping.
    random_state:
        Seed.
    """

    def __init__(
        self,
        hidden_layer_sizes: tuple = (128, 64),
        activation: str = "relu",
        alpha: float = 1e-4,
        learning_rate_init: float = 1e-3,
        max_iter: int = 500,
        early_stopping: bool = True,
        random_state: int = 42,
    ):
        self.hidden_layer_sizes = hidden_layer_sizes
        self.activation = activation
        self.alpha = alpha
        self.learning_rate_init = learning_rate_init
        self.max_iter = max_iter
        self.early_stopping = early_stopping
        self.random_state = random_state

    def fit(self, Z: np.ndarray, y: np.ndarray) -> "MLPMutationClassifier":
        self._model = MLPClassifier(
            hidden_layer_sizes=self.hidden_layer_sizes,
            activation=self.activation,
            alpha=self.alpha,
            learning_rate_init=self.learning_rate_init,
            max_iter=self.max_iter,
            early_stopping=self.early_stopping,
            random_state=self.random_state,
        )
        self._model.fit(Z, y)
        return self

    def predict_proba(self, Z: np.ndarray) -> np.ndarray:
        return self._model.predict_proba(Z)


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------

_CLASSIFIER_REGISTRY = {
    "logistic": LogisticMutationClassifier,
    "random_forest": RandomForestMutationClassifier,
    "gbm": GBMMutationClassifier,
    "xgboost": XGBMutationClassifier,
    "lightgbm": LGBMMutationClassifier,
    "svm": SVMMutationClassifier,
    "mlp": MLPMutationClassifier,
}


def get_classifier(name: str, **kwargs) -> BaseMutationClassifier:
    """Factory: instantiate a classifier by name.

    Parameters
    ----------
    name:
        One of ``"logistic"``, ``"random_forest"``, ``"gbm"``, ``"xgboost"``,
        ``"lightgbm"``, ``"svm"``, ``"mlp"``.
    **kwargs:
        Passed to the constructor.

    Returns
    -------
    BaseMutationClassifier
    """
    name = name.lower()
    if name not in _CLASSIFIER_REGISTRY:
        raise ValueError(
            f"Unknown classifier '{name}'. Choose from: {list(_CLASSIFIER_REGISTRY)}"
        )
    return _CLASSIFIER_REGISTRY[name](**kwargs)
