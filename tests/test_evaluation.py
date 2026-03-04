"""Tests for evaluation metrics."""

import numpy as np
import pandas as pd
import pytest

from scope.classification.models import LogisticMutationClassifier
from scope.evaluation.metrics import (
    cross_validate_classifiers,
    evaluate_all,
    evaluate_classifier,
    pr_curve_data,
    roc_curve_data,
)


@pytest.fixture()
def binary_preds(rng):
    y_true = rng.binomial(1, 0.4, 80)
    y_prob = np.clip(y_true + rng.normal(0, 0.3, 80), 0, 1)
    return y_true, y_prob


class TestEvaluateClassifier:
    def test_returns_metrics(self, binary_preds):
        y_true, y_prob = binary_preds
        m = evaluate_classifier(y_true, y_prob)
        assert 0 <= m["auroc"] <= 1
        assert 0 <= m["auprc"] <= 1
        assert 0 <= m["brier"] <= 1

    def test_degenerate_all_pos(self):
        m = evaluate_classifier(np.ones(20), np.ones(20))
        assert np.isnan(m["auroc"])

    def test_roc_curve(self, binary_preds):
        y_true, y_prob = binary_preds
        fpr, tpr, auroc = roc_curve_data(y_true, y_prob)
        assert fpr[0] == 0.0
        assert tpr[-1] == 1.0
        assert 0 <= auroc <= 1

    def test_pr_curve(self, binary_preds):
        y_true, y_prob = binary_preds
        precision, recall, auprc = pr_curve_data(y_true, y_prob)
        assert 0 <= auprc <= 1


class TestEvaluateAll:
    def test_evaluate_all(self, rng):
        n = 60
        y_true_df = pd.DataFrame(
            {"KRAS": rng.binomial(1, 0.4, n), "TP53": rng.binomial(1, 0.5, n)}
        )
        y_prob_df = pd.DataFrame(
            {
                "mutation_prob_KRAS": np.clip(rng.uniform(0, 1, n), 0, 1),
                "mutation_prob_TP53": np.clip(rng.uniform(0, 1, n), 0, 1),
            }
        )
        result = evaluate_all(y_true_df, y_prob_df)
        assert "KRAS" in result.index
        assert "auroc" in result.columns


class TestCrossValidate:
    def test_returns_df(self, rng):
        n = 60
        Z = rng.normal(size=(n, 10)).astype(np.float32)
        labels = pd.DataFrame(
            {"KRAS": rng.binomial(1, 0.4, n)},
            index=[f"S{i}" for i in range(n)],
        )
        result = cross_validate_classifiers(
            Z, labels, LogisticMutationClassifier, cv=3
        )
        assert "auroc" in result.columns
        assert len(result) == 3  # 3 folds × 1 mutation
