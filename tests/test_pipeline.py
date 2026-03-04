"""End-to-end integration tests for BulkPipeline and SingleCellPipeline."""

import numpy as np
import pytest

from scope.pipeline import BulkPipeline, SingleCellPipeline


class TestBulkPipeline:
    @pytest.mark.parametrize("decomp", ["svd", "pca", "nmf"])
    @pytest.mark.parametrize("clf", ["logistic", "random_forest"])
    def test_fit(self, adata_bulk, mutation_labels, decomp, clf):
        pipe = BulkPipeline(
            norm_method="cpm",
            decomposition=decomp,
            n_components=10,
            classifier=clf,
        )
        pipe.fit(adata_bulk, mutation_labels)
        assert hasattr(pipe, "classifier_set_")
        assert len(pipe.classifier_set_.classifiers_) > 0

    def test_save_load(self, adata_bulk, mutation_labels, tmp_path):
        pipe = BulkPipeline(n_components=5)
        pipe.fit(adata_bulk, mutation_labels)
        save_path = tmp_path / "pipeline.pkl"
        pipe.save(save_path)
        loaded = BulkPipeline.load(save_path)
        prob_orig = pipe.predict_bulk(adata_bulk)
        prob_load = loaded.predict_bulk(adata_bulk)
        import pandas as pd
        pd.testing.assert_frame_equal(prob_orig, prob_load)

    def test_cv(self, adata_bulk, mutation_labels):
        pipe = BulkPipeline(n_components=5)
        pipe.fit(adata_bulk, mutation_labels, cv=3)
        assert pipe.cv_results_ is not None
        assert len(pipe.cv_results_) > 0


class TestSingleCellPipeline:
    @pytest.fixture()
    def fitted_bulk_pipe(self, adata_bulk, mutation_labels):
        pipe = BulkPipeline(n_components=10, classifier="logistic")
        pipe.fit(adata_bulk, mutation_labels)
        return pipe

    def test_full_transform(self, fitted_bulk_pipe, adata_bulk, adata_sc):
        adata_bulk_pp = fitted_bulk_pipe.preprocessor_.transform(adata_bulk)

        sc_pipe = SingleCellPipeline(
            bulk_pipeline=fitted_bulk_pipe,
            alignment_method="z_score_bulk",
        )
        sc_pipe.fit(adata_bulk_pp, adata_sc)
        result = sc_pipe.transform(adata_sc)

        # Check latent embedding was added
        assert fitted_bulk_pipe.obsm_key_ in result.obsm

        # Check mutation probability columns were added to obs
        prob_cols = [c for c in result.obs.columns if c.startswith("mutation_prob_")]
        assert len(prob_cols) > 0

        # Probabilities should be in [0, 1]
        for col in prob_cols:
            vals = result.obs[col].values
            assert vals.min() >= 0.0 - 1e-6
            assert vals.max() <= 1.0 + 1e-6

    def test_partial_gene_overlap(self, fitted_bulk_pipe, adata_bulk, rng):
        """Ensure sc data with only partial gene overlap is handled gracefully."""
        import anndata as ad

        half_genes = list(adata_bulk.var_names[:250])
        X = rng.negative_binomial(5, 0.6, size=(50, 250)).astype(np.float32)
        adata_sc_partial = ad.AnnData(X=X)
        adata_sc_partial.obs_names = [f"PC{i}" for i in range(50)]
        adata_sc_partial.var_names = half_genes

        adata_bulk_pp = fitted_bulk_pipe.preprocessor_.transform(adata_bulk)
        sc_pipe = SingleCellPipeline(
            bulk_pipeline=fitted_bulk_pipe,
            alignment_method="none",
        )
        sc_pipe.fit(adata_bulk_pp, adata_sc_partial)
        result = sc_pipe.transform(adata_sc_partial)
        assert result.n_obs == 50
