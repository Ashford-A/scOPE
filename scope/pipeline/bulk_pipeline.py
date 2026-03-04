"""BulkPipeline: orchestrates Phase 1 of scOPE.

Phase 1 learns a latent space from bulk RNA-seq and trains per-mutation
classifiers in that space.

Typical usage
-------------
>>> bulk_pipe = BulkPipeline(
...     norm_method="cpm",
...     decomposition="svd",
...     n_components=50,
...     classifier="logistic",
... )
>>> bulk_pipe.fit(adata_bulk, mutation_labels)
"""

from __future__ import annotations

import pickle
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn.base import BaseEstimator

from scope.classification.base import PerMutationClassifierSet
from scope.classification.models import get_classifier
from scope.decomposition import BaseDecomposition, get_decomposition
from scope.evaluation.metrics import cross_validate_classifiers, evaluate_all
from scope.preprocessing.bulk import BulkPreprocessor
from scope.utils.logging import get_logger

log = get_logger(__name__)

PathLike = Union[str, Path]


class BulkPipeline(BaseEstimator):
    """End-to-end bulk RNA-seq phase of scOPE.

    Parameters
    ----------
    norm_method:
        Bulk normalisation strategy: ``"cpm"``, ``"tpm"``, ``"median_ratio"``,
        ``"tmm"``, or ``"none"``.
    log1p:
        Log(x+1) transform after normalisation.
    center, scale:
        Gene-wise centering / scaling.
    decomposition:
        Decomposition method: ``"svd"`` (default), ``"nmf"``, ``"ica"``,
        ``"pca"``.
    n_components:
        Latent dimension.
    decomposition_kwargs:
        Extra keyword arguments forwarded to the decomposition constructor.
    classifier:
        Classifier name: ``"logistic"``, ``"random_forest"``, ``"gbm"``,
        ``"xgboost"``, ``"lightgbm"``, ``"svm"``, ``"mlp"``.
    classifier_kwargs:
        Extra keyword arguments forwarded to the classifier constructor.
    min_positive_frac:
        Skip classifier training for mutations with fewer positive samples
        than this fraction.
    scale_features:
        Prepend a standard scaler to each classifier.
    gene_lengths:
        Array of gene lengths (bp) required only for ``norm_method="tpm"``.
    layer:
        AnnData layer to use as the expression matrix.  ``None`` → ``adata.X``.
    """

    def __init__(
        self,
        norm_method: str = "cpm",
        log1p: bool = True,
        center: bool = True,
        scale: bool = True,
        decomposition: str = "svd",
        n_components: int = 50,
        decomposition_kwargs: dict | None = None,
        classifier: str = "logistic",
        classifier_kwargs: dict | None = None,
        min_positive_frac: float = 0.05,
        scale_features: bool = True,
        gene_lengths: np.ndarray | None = None,
        layer: str | None = None,
    ):
        self.norm_method = norm_method
        self.log1p = log1p
        self.center = center
        self.scale = scale
        self.decomposition = decomposition
        self.n_components = n_components
        self.decomposition_kwargs = decomposition_kwargs or {}
        self.classifier = classifier
        self.classifier_kwargs = classifier_kwargs or {}
        self.min_positive_frac = min_positive_frac
        self.scale_features = scale_features
        self.gene_lengths = gene_lengths
        self.layer = layer

    # ------------------------------------------------------------------
    def fit(
        self,
        adata_bulk: AnnData,
        mutation_labels: pd.DataFrame,
        cv: int | None = None,
    ) -> BulkPipeline:
        """Learn latent space and train mutation classifiers.

        Parameters
        ----------
        adata_bulk:
            Raw (or pre-filtered) bulk RNA-seq AnnData.
            ``obs_names`` must align with ``mutation_labels.index``.
        mutation_labels:
            Binary DataFrame (samples × mutations).  Index must match
            ``adata_bulk.obs_names``.
        cv:
            If provided, run stratified k-fold cross-validation and store
            results in ``self.cv_results_``.

        Returns
        -------
        self
        """
        log.info("=== BulkPipeline.fit ===")

        # ── 1. Align sample order ──────────────────────────────────────────
        shared_samples = adata_bulk.obs_names.intersection(mutation_labels.index)
        if len(shared_samples) == 0:
            raise ValueError(
                "No overlapping sample IDs between adata_bulk.obs_names and "
                "mutation_labels.index."
            )
        if len(shared_samples) < adata_bulk.n_obs:
            log.warning(
                "Subsetting to %d / %d samples with mutation labels.",
                len(shared_samples), adata_bulk.n_obs,
            )
        adata_bulk = adata_bulk[shared_samples].copy()
        mutation_labels = mutation_labels.loc[shared_samples]

        # ── 2. Preprocessing ──────────────────────────────────────────────
        log.info("Preprocessing bulk data (norm=%s, log1p=%s).", self.norm_method, self.log1p)
        self.preprocessor_ = BulkPreprocessor(
            norm_method=self.norm_method,
            log1p=self.log1p,
            center=self.center,
            scale=self.scale,
            gene_lengths=self.gene_lengths,
            layer_in=self.layer,
            layer_out=self.layer,
        )
        adata_pp = self.preprocessor_.fit_transform(adata_bulk)

        # ── 3. Decomposition ──────────────────────────────────────────────
        log.info("Decomposition: %s (k=%d).", self.decomposition, self.n_components)
        self.decomposer_: BaseDecomposition = get_decomposition(
            self.decomposition,
            n_components=self.n_components,
            layer=self.layer,
            **self.decomposition_kwargs,
        )
        adata_pp = self.decomposer_.fit_transform(adata_pp)
        Z_bulk = adata_pp.obsm[self.decomposer_.obsm_key]
        self.obsm_key_ = self.decomposer_.obsm_key

        # Store gene names for downstream checks
        self.gene_names_ = list(adata_bulk.var_names)

        # ── 4. Classifiers ────────────────────────────────────────────────
        log.info("Training classifiers (%s).", self.classifier)

        def _factory():
            return get_classifier(self.classifier, **self.classifier_kwargs)

        self.classifier_set_ = PerMutationClassifierSet(
            classifier_factory=_factory,
            min_positive_frac=self.min_positive_frac,
            scale_features=self.scale_features,
        )
        self.classifier_set_.fit(Z_bulk, mutation_labels)

        # ── 5. (Optional) cross-validation ────────────────────────────────
        if cv is not None:
            log.info("Running %d-fold cross-validation.", cv)
            self.cv_results_ = cross_validate_classifiers(
                Z_bulk, mutation_labels, _factory, cv=cv
            )
        else:
            self.cv_results_ = None

        self.mutation_labels_ = mutation_labels
        self.is_fitted_ = True
        log.info("BulkPipeline.fit complete.")
        return self

    # ------------------------------------------------------------------
    def transform_bulk(self, adata_bulk: AnnData) -> AnnData:
        """Preprocess and embed a new bulk dataset into the learned space."""
        adata_pp = self.preprocessor_.transform(adata_bulk)
        adata_pp = self.decomposer_.transform(adata_pp)
        return adata_pp

    def predict_bulk(self, adata_bulk: AnnData) -> pd.DataFrame:
        """Return per-mutation probabilities for new bulk samples."""
        adata_emb = self.transform_bulk(adata_bulk)
        Z = adata_emb.obsm[self.obsm_key_]
        return self.classifier_set_.predict_proba(Z)

    # ------------------------------------------------------------------
    def evaluate(self, adata_bulk: AnnData, mutation_labels: pd.DataFrame) -> pd.DataFrame:
        """Evaluate on a held-out bulk dataset."""
        proba_df = self.predict_bulk(adata_bulk)
        shared = adata_bulk.obs_names.intersection(mutation_labels.index)
        labels_sub = mutation_labels.loc[shared]
        proba_sub = proba_df.loc[proba_df.index.isin(shared)] if hasattr(proba_df.index, 'isin') else proba_df
        return evaluate_all(labels_sub, proba_sub)

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self, path: PathLike) -> None:
        """Serialise the fitted pipeline to *path* with pickle.

        Parameters
        ----------
        path:
            File path (e.g. ``"bulk_pipeline.pkl"``).
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "wb") as fh:
            pickle.dump(self, fh)
        log.info("BulkPipeline saved to '%s'.", path)

    @classmethod
    def load(cls, path: PathLike) -> BulkPipeline:
        """Load a previously saved pipeline.

        Parameters
        ----------
        path:
            Path to the ``.pkl`` file.
        """
        with open(str(path), "rb") as fh:
            obj = pickle.load(fh)
        log.info("BulkPipeline loaded from '%s'.", path)
        return obj
