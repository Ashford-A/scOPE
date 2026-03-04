"""Shared fixtures for scOPE tests."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


@pytest.fixture(scope="session")
def rng():
    return np.random.default_rng(42)


@pytest.fixture(scope="session")
def n_bulk_samples():
    return 80


@pytest.fixture(scope="session")
def n_sc_cells():
    return 200


@pytest.fixture(scope="session")
def n_genes():
    return 500


@pytest.fixture(scope="session")
def gene_names(n_genes):
    return [f"GENE{i:04d}" for i in range(n_genes)]


@pytest.fixture(scope="session")
def adata_bulk(rng, n_bulk_samples, n_genes, gene_names):
    """Synthetic bulk RNA-seq AnnData (samples × genes, raw counts)."""
    X = rng.negative_binomial(10, 0.5, size=(n_bulk_samples, n_genes)).astype(np.float32)
    adata = AnnData(X=X)
    adata.obs_names = [f"SAMPLE{i:03d}" for i in range(n_bulk_samples)]
    adata.var_names = gene_names
    return adata


@pytest.fixture(scope="session")
def adata_sc(rng, n_sc_cells, n_genes, gene_names):
    """Synthetic scRNA-seq AnnData (cells × genes, raw counts)."""
    X = rng.negative_binomial(5, 0.6, size=(n_sc_cells, n_genes)).astype(np.float32)
    adata = AnnData(X=X)
    adata.obs_names = [f"CELL{i:04d}" for i in range(n_sc_cells)]
    adata.var_names = gene_names
    return adata


@pytest.fixture(scope="session")
def mutation_labels(rng, n_bulk_samples, adata_bulk):
    """Binary mutation label DataFrame aligned to adata_bulk."""
    df = pd.DataFrame(
        {
            "KRAS": rng.binomial(1, 0.35, n_bulk_samples),
            "TP53": rng.binomial(1, 0.50, n_bulk_samples),
            "EGFR": rng.binomial(1, 0.15, n_bulk_samples),
        },
        index=adata_bulk.obs_names,
    )
    return df
