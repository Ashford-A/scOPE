"""Tests for utility modules."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scope.utils.validation import (
    check_adata,
    check_gene_overlap,
    check_is_fitted,
    check_mutation_labels,
    check_nonneg,
)
from scope.utils.gene_utils import (
    align_gene_order,
    filter_variable_genes,
    get_shared_genes,
    subset_to_shared_genes,
)


class TestCheckAdata:
    def test_valid(self, adata_bulk):
        check_adata(adata_bulk)  # should not raise

    def test_wrong_type(self):
        with pytest.raises(TypeError, match="AnnData"):
            check_adata(np.zeros((10, 5)))

    def test_too_few_obs(self):
        tiny = AnnData(np.zeros((1, 100)))
        with pytest.raises(ValueError, match="observations"):
            check_adata(tiny, min_obs=2)

    def test_duplicate_obs_names(self):
        adata = AnnData(np.zeros((4, 10)))
        adata.obs_names = ["A", "A", "B", "C"]
        with pytest.raises(ValueError, match="duplicate"):
            check_adata(adata)


class TestCheckMutationLabels:
    def test_valid(self, mutation_labels, adata_bulk):
        check_mutation_labels(mutation_labels, adata=adata_bulk)

    def test_wrong_type(self):
        with pytest.raises(TypeError, match="pd.DataFrame"):
            check_mutation_labels(np.zeros((10, 2)))

    def test_no_overlap(self, mutation_labels):
        adata = AnnData(np.zeros((5, 50)))
        adata.obs_names = ["NOMATCH_0", "NOMATCH_1", "NOMATCH_2", "NOMATCH_3", "NOMATCH_4"]
        with pytest.raises(ValueError, match="No overlapping"):
            check_mutation_labels(mutation_labels, adata=adata)


class TestCheckIsFixed:
    def test_fitted(self):
        class MyEstimator:
            coef_ = np.array([1.0])
        check_is_fitted(MyEstimator(), ["coef_"])

    def test_not_fitted(self):
        from sklearn.exceptions import NotFittedError
        class MyEstimator:
            pass
        with pytest.raises(NotFittedError):
            check_is_fitted(MyEstimator(), ["coef_"])


class TestCheckNonneg:
    def test_nonneg_ok(self):
        check_nonneg(np.array([0.0, 1.0, 2.0]))

    def test_negative_warns(self):
        with pytest.warns(UserWarning, match="negative"):
            check_nonneg(np.array([-1.0, 0.0, 1.0]))

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="negative"):
            check_nonneg(np.array([-1.0, 0.0]), raise_on_negative=True)


class TestGeneUtils:
    def test_get_shared_genes(self, adata_bulk, adata_sc):
        shared = get_shared_genes(adata_bulk, adata_sc)
        assert len(shared) > 0
        assert all(g in adata_bulk.var_names for g in shared)
        assert all(g in adata_sc.var_names for g in shared)

    def test_subset_to_shared(self, adata_bulk, adata_sc):
        b, s = subset_to_shared_genes(adata_bulk, adata_sc)
        assert list(b.var_names) == list(s.var_names)

    def test_align_gene_order(self):
        X = np.arange(6).reshape(2, 3).astype(float)
        src = ["A", "B", "C"]
        tgt = ["C", "A", "D"]
        out = align_gene_order(X, src, tgt)
        assert out.shape == (2, 3)
        np.testing.assert_array_equal(out[:, 0], X[:, 2])  # C
        np.testing.assert_array_equal(out[:, 1], X[:, 0])  # A
        np.testing.assert_array_equal(out[:, 2], 0.0)      # D (missing → 0)

    def test_filter_variable_genes(self, adata_bulk):
        selected = filter_variable_genes(adata_bulk, n_top_genes=100)
        assert len(selected) == 100
        assert all(g in adata_bulk.var_names for g in selected)

    def test_gene_overlap_raises_on_low(self):
        with pytest.raises(ValueError, match="shared"):
            check_gene_overlap(["A", "B"], ["C", "D", "E"], min_overlap=1)
