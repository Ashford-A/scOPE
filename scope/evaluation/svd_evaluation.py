"""SVD component evaluation and visualisation for scOPE.

Provides :class:`SVDEvaluator`, which takes a fitted
:class:`~scope.pipeline.bulk_pipeline.BulkPipeline` and a latent embedding
``Z_bulk`` and produces a comprehensive suite of plots and tables that
explain *which* SVD components drive per-mutation classification and *which
genes* define those components.

Typical usage
-------------
>>> from scope.evaluation.svd_evaluation import SVDEvaluator
>>>
>>> # After bulk_pipe.fit(adata_bulk, mutation_labels) ...
>>> adata_pp = bulk_pipe.transform_bulk(adata_bulk)
>>> Z_bulk = adata_pp.obsm[bulk_pipe.obsm_key_]
>>>
>>> ev = SVDEvaluator(bulk_pipe, Z_bulk, mutation="KRAS")
>>> ev.run_all(output_dir="figures/svd_eval_KRAS")
"""

from __future__ import annotations

import contextlib
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, permutation_test_score
from sklearn.pipeline import Pipeline

matplotlib.use("Agg")

from scope.decomposition.svd import SVDDecomposition
from scope.evaluation.metrics import roc_curve_data
from scope.pipeline.bulk_pipeline import BulkPipeline
from scope.utils.logging import get_logger

log = get_logger(__name__)

# Optional deps -----------------------------------------------------------
try:
    import shap as _shap_lib

    _HAS_SHAP = True
except ImportError:
    _HAS_SHAP = False

try:
    import umap as _umap_lib

    _HAS_UMAP = True
except ImportError:
    _HAS_UMAP = False

try:
    from statsmodels.stats.multitest import multipletests

    _HAS_STATSMODELS = True
except ImportError:
    _HAS_STATSMODELS = False

_PAL = sns.color_palette("Set2")
sns.set_theme(style="whitegrid", font_scale=1.1)


# =========================================================================
class SVDEvaluator:
    """Interpret SVD components in the context of a trained mutation classifier.

    Computes and plots:

    * Component importance (|LR coef| × σ) — scree overlay
    * Gene loading heatmap for top-k important components
    * Signed gene bar charts per component
    * Pairwise latent scatter coloured by mutation label
    * Violin plots of component scores by mutation status (Mann-Whitney p)
    * Spearman component-label correlation heatmap (FDR-corrected)
    * ROC curve + component ablation AUC curve
    * Permutation test (real AUC vs. null)
    * Gene biplot (sample scores + gene arrows)
    * SHAP summary (requires ``shap``)
    * UMAP of Z_bulk (requires ``umap-learn``)
    * Calibration curve
    * Component cross-correlation (detects batch bleed)
    * Gene program CSV table

    Parameters
    ----------
    bulk_pipeline:
        A fitted :class:`~scope.pipeline.bulk_pipeline.BulkPipeline` whose
        decomposition is ``"svd"``.
    Z_bulk:
        Bulk latent embedding, shape ``(n_samples, k)``.  Obtain via::

            adata_pp = bulk_pipe.transform_bulk(adata_bulk)
            Z_bulk = adata_pp.obsm[bulk_pipe.obsm_key_]
    mutation:
        Name of the mutation to evaluate (must be a key in
        ``bulk_pipeline.classifier_set_.classifiers_``).
    """

    def __init__(
        self,
        bulk_pipeline: BulkPipeline,
        Z_bulk: np.ndarray,
        mutation: str,
    ):
        if bulk_pipeline.decomposition != "svd":
            raise ValueError(
                f"SVDEvaluator requires decomposition='svd', "
                f"got '{bulk_pipeline.decomposition}'."
            )
        self._bp = bulk_pipeline
        self._svd: SVDDecomposition = bulk_pipeline.decomposer_  # type: ignore[assignment]
        self.Z_bulk = np.asarray(Z_bulk, dtype=np.float64)
        self.mutation = mutation

        # Ground-truth labels
        self.y: np.ndarray = bulk_pipeline.mutation_labels_[mutation].values.astype(int)
        # Gene names and loadings
        self.gene_names: list[str] = bulk_pipeline.gene_names_
        # V: (n_genes, k)
        self.V: np.ndarray = self._svd.components_.T
        self.Sigma: np.ndarray = self._svd.singular_values_
        self.evr: np.ndarray = self._svd._explained_variance_ratio_
        self.k: int = self.Sigma.shape[0]

        # Extract LR coefficients from the per-mutation classifier.
        # The stored object may be a raw estimator or a sklearn Pipeline.
        clf_obj = bulk_pipeline.classifier_set_.classifiers_[mutation]
        self._clf = clf_obj
        coef_raw = self._extract_coef(clf_obj)  # (1, k) or (k,)
        self.coef_abs: np.ndarray = np.abs(coef_raw).mean(axis=0)  # (k,)
        self.importance: np.ndarray = self.coef_abs * self.Sigma  # weighted

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_coef(clf) -> np.ndarray:
        """Return coef_ from a plain estimator or sklearn Pipeline."""
        if isinstance(clf, Pipeline):
            # Walk through steps to find the final estimator with coef_
            for _, step in reversed(clf.steps):
                if hasattr(step, "coef_"):
                    return step.coef_
            raise AttributeError("No step with coef_ found in Pipeline.")
        if hasattr(clf, "coef_"):
            return clf.coef_
        raise AttributeError(
            f"Classifier {type(clf).__name__} has no coef_. "
            "Use logistic/SVM classifiers for SVD interpretation, "
            "or use SHAP for tree-based models."
        )

    def _top_genes(self, component: int, n: int = 20) -> pd.DataFrame:
        loadings = self.V[:, component]
        idx = np.argsort(np.abs(loadings))[::-1][:n]
        return pd.DataFrame(
            {
                "gene": [self.gene_names[i] for i in idx],
                "loading": loadings[idx],
            }
        )

    def _ranked_components(self, top_k: int) -> np.ndarray:
        return np.argsort(self.importance)[::-1][:top_k]

    def _savefig(self, fig, name: str, out: Path, dpi: int = 150) -> None:
        out.mkdir(parents=True, exist_ok=True)
        p = out / name
        fig.savefig(p, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        log.info("  saved → %s", p)

    def _cv_proba(self, Z: np.ndarray, cv: int = 5) -> np.ndarray:
        """Cross-validated predicted probabilities for Z."""
        skf = StratifiedKFold(n_splits=cv, shuffle=True, random_state=42)
        y_prob = np.zeros(len(self.y))
        for tr, te in skf.split(Z, self.y):
            clf = LogisticRegression(max_iter=1000, random_state=42)
            clf.fit(Z[tr], self.y[tr])
            p = clf.predict_proba(Z[te])
            y_prob[te] = p[:, 1] if p.shape[1] == 2 else p.max(1)
        return y_prob

    # ------------------------------------------------------------------
    # 1. Classifier-weighted scree plot
    # ------------------------------------------------------------------

    def plot_weighted_scree(self, output_dir: Path, top_k: int = 40) -> None:
        """Scree plot with bars colour-coded by LR component importance.

        Extends the existing :func:`~scope.visualization.embeddings.plot_scree`
        by overlaying |coef|×σ importance as bar colour.
        """
        k = min(self.k, top_k)
        evr = self.evr[:k]
        cumevr = np.cumsum(evr)
        imp = self.importance[:k]
        imp_norm = imp / (imp.max() + 1e-12)

        fig, ax1 = plt.subplots(figsize=(11, 4))
        cmap = plt.cm.RdYlBu_r
        ax1.bar(
            np.arange(1, k + 1),
            evr * 100,
            color=cmap(imp_norm),
            edgecolor="white",
            linewidth=0.3,
        )
        ax1.set_xlabel("SVD component")
        ax1.set_ylabel("Explained variance (%)", color="k")
        ax1.set_xticks(
            np.arange(1, k + 1, max(1, k // 20)),
        )

        ax2 = ax1.twinx()
        ax2.plot(
            np.arange(1, k + 1), cumevr * 100, "k--o", ms=3, lw=1.5, label="Cumul. EVR"
        )
        ax2.set_ylabel("Cumulative EVR (%)")
        for thresh in [80, 90]:
            ax2.axhline(thresh, color="lightgrey", ls=":", lw=1)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, imp.max()))
        sm.set_array([])
        plt.colorbar(sm, ax=ax1, label=f"|coef|×σ  ({self.mutation})")
        ax1.set_title(
            f"Classifier-weighted scree — {self.mutation}\n"
            f"bar colour = importance for LR prediction"
        )
        self._savefig(fig, "weighted_scree.png", output_dir)

    # ------------------------------------------------------------------
    # 2. Component importance bar chart
    # ------------------------------------------------------------------

    def plot_component_importance(self, output_dir: Path, top_k: int = 20) -> None:
        """Ranked bar chart of |coef|×σ and decomposed |coef| vs σ."""
        order = self._ranked_components(top_k)
        labels = [f"SVD{i + 1}" for i in order]

        fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(12, 4), sharey=True)

        ax_l.barh(labels[::-1], self.importance[order][::-1], color=_PAL[0])
        ax_l.set_xlabel("|coef| × σ")
        ax_l.set_title(f"Weighted importance  ({self.mutation})")

        x = np.arange(top_k)
        ax_r.bar(x, self.coef_abs[order][::-1], label="|LR coef|", color=_PAL[1])
        sigma_norm = self.Sigma[order] / (self.Sigma.max() + 1e-12)
        ax_r.bar(x, sigma_norm[::-1], alpha=0.45, label="σ (norm.)", color=_PAL[2])
        ax_r.set_xticks(x)
        ax_r.set_xticklabels(labels[::-1], rotation=45, ha="right")
        ax_r.set_title("Raw |coef| vs σ")
        ax_r.legend(frameon=False)

        fig.suptitle(f"SVD Component Importance — {self.mutation}", y=1.02)
        self._savefig(fig, "component_importance.png", output_dir)

    # ------------------------------------------------------------------
    # 3. Gene loading heatmap
    # ------------------------------------------------------------------

    def plot_gene_loading_heatmap(
        self,
        output_dir: Path,
        top_components: int = 10,
        top_genes_per_comp: int = 10,
    ) -> None:
        """Hierarchically clustered heatmap of gene loadings × top components."""
        from scipy.cluster import hierarchy as sch

        order = self._ranked_components(top_components)
        gene_set: dict[str, bool] = {}
        for ci in order:
            for g in self._top_genes(ci, n=top_genes_per_comp)["gene"]:
                gene_set[g] = True
        genes = list(gene_set)
        g_idx = [self.gene_names.index(g) for g in genes]

        mat = self.V[np.ix_(g_idx, order)]  # (genes_sub, top_k)
        linkage = sch.linkage(mat, method="ward")
        row_order = sch.leaves_list(linkage)
        mat = mat[row_order]
        genes_ord = [genes[r] for r in row_order]
        col_labels = [f"SVD{i + 1}" for i in order]

        fig_h = max(7, len(genes) * 0.32)
        fig_w = max(7, top_components * 0.8)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        sns.heatmap(
            mat,
            xticklabels=col_labels,
            yticklabels=genes_ord,
            center=0,
            cmap="RdBu_r",
            linewidths=0.25,
            linecolor="white",
            cbar_kws={"shrink": 0.6, "label": "Gene loading  V"},
            ax=ax,
        )
        ax.set_title(
            f"Gene loadings — top {top_components} LR-important components\n"
            f"({self.mutation})"
        )
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
        plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)
        self._savefig(fig, "gene_loading_heatmap.png", output_dir)

    # ------------------------------------------------------------------
    # 4. Per-component signed gene bar charts
    # ------------------------------------------------------------------

    def plot_top_genes_per_component(
        self,
        output_dir: Path,
        top_components: int = 6,
        top_genes: int = 20,
    ) -> None:
        """One subplot per top-important component with signed gene loadings."""
        order = self._ranked_components(top_components)
        ncols = 3
        nrows = int(np.ceil(top_components / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5, nrows * 5))
        axes = np.array(axes).flatten()

        for ax_i, ci in enumerate(order):
            df = self._top_genes(ci, n=top_genes)
            colours = ["#e74c3c" if v > 0 else "#3498db" for v in df["loading"]]
            axes[ax_i].barh(df["gene"][::-1], df["loading"][::-1], color=colours[::-1])
            axes[ax_i].axvline(0, color="black", lw=0.7)
            axes[ax_i].set_title(
                f"SVD{ci + 1}  imp={self.importance[ci]:.3g}", fontsize=9
            )
            axes[ax_i].set_xlabel("Gene loading")

        for ax_i in range(top_components, len(axes)):
            axes[ax_i].set_visible(False)

        fig.suptitle(f"Top gene loadings per component  ({self.mutation})", y=1.01)
        plt.tight_layout()
        self._savefig(fig, "top_genes_per_component.png", output_dir)

    # ------------------------------------------------------------------
    # 5. Latent scatter
    # ------------------------------------------------------------------

    def plot_latent_scatter(self, output_dir: Path, n_pairs: int = 3) -> None:
        """Pairwise scatter of top components, coloured by mutation label."""
        order = self._ranked_components(n_pairs + 1)
        pairs = [
            (order[i], order[j])
            for i in range(len(order))
            for j in range(i + 1, len(order))
        ][:n_pairs]

        classes = np.unique(self.y)
        palette = dict(zip(classes.tolist(), _PAL))

        fig, axes = plt.subplots(1, len(pairs), figsize=(6 * len(pairs), 5))
        if len(pairs) == 1:
            axes = [axes]

        for ax, (ci, cj) in zip(axes, pairs):
            for cls in classes:
                m = self.y == cls
                ax.scatter(
                    self.Z_bulk[m, ci],
                    self.Z_bulk[m, cj],
                    c=[palette[cls]],
                    alpha=0.7,
                    s=28,
                    edgecolors="none",
                    label=f"{self.mutation}={cls}",
                )
            ax.set_xlabel(f"Z  SVD{ci + 1}")
            ax.set_ylabel(f"Z  SVD{cj + 1}")
            ax.legend(markerscale=1.5, frameon=False)

        fig.suptitle(f"Bulk latent space — top components  ({self.mutation})")
        plt.tight_layout()
        self._savefig(fig, "latent_scatter.png", output_dir)

    # ------------------------------------------------------------------
    # 6. Separation violins
    # ------------------------------------------------------------------

    def plot_separation_violins(
        self, output_dir: Path, top_components: int = 12
    ) -> None:
        """Violin + Mann-Whitney p per component, split by mutation label."""
        order = self._ranked_components(top_components)
        ncols = 4
        nrows = int(np.ceil(top_components / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.5))
        axes = np.array(axes).flatten()
        classes = np.unique(self.y)

        for ax_i, ci in enumerate(order):
            groups = [self.Z_bulk[self.y == cls, ci] for cls in classes]
            parts = axes[ax_i].violinplot(groups, showmedians=True, showextrema=False)
            for pc, col in zip(parts["bodies"], _PAL):
                pc.set_facecolor(col)
                pc.set_alpha(0.72)
            axes[ax_i].set_xticks(range(1, len(classes) + 1))
            axes[ax_i].set_xticklabels([str(c) for c in classes])
            axes[ax_i].set_title(f"SVD{ci + 1}", fontsize=9)
            axes[ax_i].set_ylabel("Z score")

            if len(groups) == 2:
                _, pval = stats.mannwhitneyu(
                    groups[0], groups[1], alternative="two-sided"
                )
                ptext = f"p={pval:.1e}" if pval < 0.05 else f"p={pval:.2f}"
                axes[ax_i].set_xlabel(ptext, fontsize=8)

        for ax_i in range(top_components, len(axes)):
            axes[ax_i].set_visible(False)

        fig.suptitle(f"Component score distributions by {self.mutation} status", y=1.01)
        plt.tight_layout()
        self._savefig(fig, "separation_violins.png", output_dir)

    # ------------------------------------------------------------------
    # 7. Component–label Spearman correlation heatmap
    # ------------------------------------------------------------------

    def plot_component_label_correlation(
        self, output_dir: Path, top_k: int = 30
    ) -> None:
        """Spearman ρ between each SVD component and mutation label (FDR-corrected)."""
        k = min(self.k, top_k)
        classes = np.unique(self.y)
        rho_mat = np.zeros((k, len(classes)))
        p_mat = np.zeros_like(rho_mat)

        for ci in range(k):
            for j, cls in enumerate(classes):
                y_bin = (self.y == cls).astype(float)
                r, p = stats.spearmanr(self.Z_bulk[:, ci], y_bin)
                rho_mat[ci, j] = r
                p_mat[ci, j] = p

        if _HAS_STATSMODELS:
            _, q_flat, _, _ = multipletests(p_mat.flatten(), method="fdr_bh")
            q_mat = q_flat.reshape(p_mat.shape)
        else:
            q_mat = p_mat

        annot = np.array(
            [
                [
                    f"{rho_mat[r, c]:.2f}\n(q={q_mat[r, c]:.2f})"
                    for c in range(len(classes))
                ]
                for r in range(k)
            ]
        )
        fig, ax = plt.subplots(figsize=(max(4, len(classes) * 1.6), max(6, k * 0.35)))
        sns.heatmap(
            rho_mat,
            xticklabels=[f"{self.mutation}={cls}" for cls in classes],
            yticklabels=[f"SVD{i + 1}" for i in range(k)],
            center=0,
            cmap="RdBu_r",
            annot=annot,
            fmt="",
            annot_kws={"fontsize": 7},
            cbar_kws={"label": "Spearman ρ"},
            ax=ax,
        )
        ax.set_title(f"Component–label Spearman ρ  ({self.mutation})")
        plt.tight_layout()
        self._savefig(fig, "component_label_correlation.png", output_dir)

    # ------------------------------------------------------------------
    # 8. ROC + ablation
    # ------------------------------------------------------------------

    def plot_roc_ablation(
        self, output_dir: Path, cv: int = 5, n_ablation_steps: int = 10
    ) -> None:
        """Full-model ROC (CV) + AUC vs. n_components retained (ablation)."""
        order = self._ranked_components(self.k)

        # ---- full model ROC
        y_prob = self._cv_proba(self.Z_bulk, cv=cv)
        fpr, tpr, full_auc = roc_curve_data(self.y, y_prob)

        # ---- ablation
        step = max(1, self.k // n_ablation_steps)
        keep_sizes = sorted(set(list(range(1, self.k + 1, step)) + [self.k]))
        abl_aucs, abl_n = [], []
        for n_keep in keep_sizes:
            Z_sub = self.Z_bulk[:, order[:n_keep]]
            fold_aucs = []
            skf = StratifiedKFold(n_splits=cv, shuffle=True, random_state=42)
            for tr, te in skf.split(Z_sub, self.y):
                clf = LogisticRegression(max_iter=1000, random_state=42)
                clf.fit(Z_sub[tr], self.y[tr])
                p = clf.predict_proba(Z_sub[te])
                with contextlib.suppress(Exception):
                    fold_aucs.append(
                        roc_auc_score(
                            self.y[te], p[:, 1] if p.shape[1] == 2 else p.max(1)
                        )
                    )
            abl_aucs.append(float(np.nanmean(fold_aucs)) if fold_aucs else np.nan)
            abl_n.append(n_keep)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        ax1.plot(fpr, tpr, lw=2, color=_PAL[0], label=f"AUC = {full_auc:.3f}")
        ax1.plot([0, 1], [0, 1], "k--", lw=1)
        ax1.set_xlabel("FPR")
        ax1.set_ylabel("TPR")
        ax1.set_title(f"ROC (CV={cv})  {self.mutation}")
        ax1.legend(frameon=False)

        ax2.plot(abl_n, abl_aucs, "o-", color=_PAL[1], lw=2, ms=5)
        ax2.axhline(full_auc, color=_PAL[0], ls="--", lw=1.5, label="Full model")
        ax2.axhline(0.5, color="lightgrey", ls=":", lw=1, label="Random")
        ax2.set_xlabel("# SVD components retained (top-ranked by importance)")
        ax2.set_ylabel("Mean CV AUROC")
        ax2.set_title(f"AUROC vs. n_components  ({self.mutation})")
        ax2.legend(frameon=False)

        plt.tight_layout()
        self._savefig(fig, "roc_ablation.png", output_dir)

    # ------------------------------------------------------------------
    # 9. Permutation test
    # ------------------------------------------------------------------

    def plot_permutation_test(
        self, output_dir: Path, n_permutations: int = 200, cv: int = 5
    ) -> None:
        """Permutation test: observed AUROC vs. null distribution."""
        clf = LogisticRegression(max_iter=1000, random_state=42)
        score, perm_scores, pvalue = permutation_test_score(
            clf,
            self.Z_bulk,
            self.y,
            scoring="roc_auc",
            cv=StratifiedKFold(cv, shuffle=True, random_state=42),
            n_permutations=n_permutations,
            random_state=42,
            n_jobs=-1,
        )
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.hist(
            perm_scores,
            bins=30,
            color=_PAL[2],
            alpha=0.75,
            label=f"Permuted (n={n_permutations})",
        )
        ax.axvline(
            score,
            color="red",
            lw=2,
            label=f"Observed AUROC = {score:.3f}  (p = {pvalue:.4f})",
        )
        ax.set_xlabel("AUROC (permuted labels)")
        ax.set_ylabel("Count")
        ax.set_title(f"Permutation test  ({self.mutation})")
        ax.legend(frameon=False)
        self._savefig(fig, "permutation_test.png", output_dir)

    # ------------------------------------------------------------------
    # 10. Gene biplot
    # ------------------------------------------------------------------

    def plot_gene_biplot(self, output_dir: Path, top_genes: int = 25) -> None:
        """Biplot: sample scores on top-2 components + gene loading arrows."""
        order = self._ranked_components(2)
        c1, c2 = int(order[0]), int(order[1])

        fig, ax = plt.subplots(figsize=(8, 7))
        classes = np.unique(self.y)
        for cls, col in zip(classes, _PAL):
            m = self.y == cls
            ax.scatter(
                self.Z_bulk[m, c1],
                self.Z_bulk[m, c2],
                c=[col],
                alpha=0.65,
                s=32,
                edgecolors="none",
                label=f"{self.mutation}={cls}",
            )

        # Top genes by combined loading magnitude
        mag = np.sqrt(self.V[:, c1] ** 2 + self.V[:, c2] ** 2)
        top_idx = np.argsort(mag)[::-1][:top_genes]
        scale = np.percentile(np.abs(self.Z_bulk[:, [c1, c2]]), 90) * 0.9
        vmax = mag[top_idx].max() + 1e-12

        for i in top_idx:
            vx = self.V[i, c1] / vmax * scale
            vy = self.V[i, c2] / vmax * scale
            ax.annotate(
                "",
                xy=(vx, vy),
                xytext=(0, 0),
                arrowprops={"arrowstyle": "->", "color": "grey", "lw": 0.7},
            )
            ax.text(
                vx * 1.06,
                vy * 1.06,
                self.gene_names[i],
                fontsize=7,
                color="dimgray",
                ha="center",
            )

        ax.set_xlabel(f"Z  SVD{c1 + 1}")
        ax.set_ylabel(f"Z  SVD{c2 + 1}")
        ax.legend(frameon=False)
        ax.set_title(f"Gene biplot — top-2 components  ({self.mutation})")
        self._savefig(fig, "gene_biplot.png", output_dir)

    # ------------------------------------------------------------------
    # 11. SHAP summary
    # ------------------------------------------------------------------

    def plot_shap(self, output_dir: Path) -> None:
        """SHAP dot + bar summary (requires ``shap``)."""
        if not _HAS_SHAP:
            log.warning("shap not installed — skipping SHAP plots.")
            return
        # Unwrap Pipeline if needed
        clf_final = self._clf
        if isinstance(clf_final, Pipeline):
            Z = clf_final[:-1].transform(self.Z_bulk)
            clf_final = clf_final[-1]
        else:
            Z = self.Z_bulk

        feature_names = [f"SVD{i + 1}" for i in range(self.k)]
        explainer = _shap_lib.LinearExplainer(clf_final, Z)
        shap_vals = explainer.shap_values(Z)
        if isinstance(shap_vals, list):
            shap_vals = np.array(shap_vals).mean(axis=0)

        _shap_lib.summary_plot(shap_vals, Z, feature_names=feature_names, show=False)
        plt.title(f"SHAP dot summary  ({self.mutation})")
        self._savefig(plt.gcf(), "shap_summary_dot.png", output_dir)

        _shap_lib.summary_plot(
            shap_vals, Z, feature_names=feature_names, plot_type="bar", show=False
        )
        plt.title(f"SHAP bar summary  ({self.mutation})")
        self._savefig(plt.gcf(), "shap_summary_bar.png", output_dir)

    # ------------------------------------------------------------------
    # 12. UMAP of Z_bulk
    # ------------------------------------------------------------------

    def plot_umap(self, output_dir: Path, n_neighbors: int = 15) -> None:
        """UMAP of Z_bulk coloured by mutation label (requires ``umap-learn``)."""
        if not _HAS_UMAP:
            log.warning("umap-learn not installed — skipping UMAP.")
            return
        reducer = _umap_lib.UMAP(n_neighbors=n_neighbors, random_state=42)
        emb = reducer.fit_transform(self.Z_bulk)

        fig, ax = plt.subplots(figsize=(7, 6))
        for cls, col in zip(np.unique(self.y), _PAL):
            m = self.y == cls
            ax.scatter(
                emb[m, 0],
                emb[m, 1],
                c=[col],
                s=30,
                alpha=0.7,
                edgecolors="none",
                label=f"{self.mutation}={cls}",
            )
        ax.set_xlabel("UMAP-1")
        ax.set_ylabel("UMAP-2")
        ax.set_title(f"UMAP of Z_bulk  ({self.mutation})")
        ax.legend(frameon=False)
        self._savefig(fig, "umap_zbulk.png", output_dir)

    # ------------------------------------------------------------------
    # 13. Calibration curve
    # ------------------------------------------------------------------

    def plot_calibration(self, output_dir: Path, cv: int = 5, n_bins: int = 10) -> None:
        """Reliability curve — checks if probabilities are well-calibrated."""
        from sklearn.calibration import calibration_curve

        y_prob = self._cv_proba(self.Z_bulk, cv=cv)
        pos_label = int(np.unique(self.y)[-1])
        frac_pos, mean_pred = calibration_curve(
            self.y, y_prob, n_bins=n_bins, pos_label=pos_label
        )
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.plot(mean_pred, frac_pos, "o-", color=_PAL[0], lw=2, label="LR (CV)")
        ax.plot([0, 1], [0, 1], "k--", lw=1, label="Perfect")
        ax.set_xlabel("Mean predicted probability")
        ax.set_ylabel("Fraction positives")
        ax.set_title(f"Calibration curve  ({self.mutation})")
        ax.legend(frameon=False)
        self._savefig(fig, "calibration_curve.png", output_dir)

    # ------------------------------------------------------------------
    # 14. Component cross-correlation
    # ------------------------------------------------------------------

    def plot_component_crosscorr(self, output_dir: Path, top_k: int = 20) -> None:
        """Pearson correlation among SVD component scores.

        SVD components are theoretically orthogonal; deviations after
        ``scale_by_singular_values=True`` and batch effects can introduce
        spurious correlations that may inflate classifier confidence.
        """
        k = min(self.k, top_k)
        order = self._ranked_components(k)
        Z_sub = self.Z_bulk[:, order]
        labels = [f"SVD{i + 1}" for i in order]
        corr = np.corrcoef(Z_sub.T)

        fig, ax = plt.subplots(figsize=(max(6, k * 0.55), max(5, k * 0.55)))
        mask = np.tril(np.ones_like(corr, dtype=bool))  # lower triangle
        sns.heatmap(
            corr,
            xticklabels=labels,
            yticklabels=labels,
            center=0,
            cmap="RdBu_r",
            mask=~mask,
            annot=(k <= 15),
            fmt=".2f",
            annot_kws={"fontsize": 7},
            cbar_kws={"label": "Pearson r"},
            ax=ax,
        )
        ax.set_title(
            "SVD component cross-correlation\n"
            "(non-zero lower triangle may indicate batch bleed)"
        )
        plt.tight_layout()
        self._savefig(fig, "component_crosscorr.png", output_dir)

    # ------------------------------------------------------------------
    # 15. Gene program table
    # ------------------------------------------------------------------

    def export_gene_program_table(
        self, output_dir: Path, top_components: int = 10, top_genes: int = 30
    ) -> pd.DataFrame:
        """Export a tidy CSV: component × gene loadings with importance metadata."""
        order = self._ranked_components(top_components)
        rows = []
        for rank, ci in enumerate(order):
            df = self._top_genes(ci, n=top_genes)
            for _, row in df.iterrows():
                rows.append(
                    {
                        "importance_rank": rank + 1,
                        "component": f"SVD{ci + 1}",
                        "sigma": float(self.Sigma[ci]),
                        "coef_abs": float(self.coef_abs[ci]),
                        "weighted_importance": float(self.importance[ci]),
                        "evr": float(self.evr[ci]),
                        "gene": row["gene"],
                        "loading": float(row["loading"]),
                        "loading_abs": float(abs(row["loading"])),
                        "loading_direction": (
                            "positive" if row["loading"] > 0 else "negative"
                        ),
                    }
                )
        out_df = pd.DataFrame(rows)
        output_dir.mkdir(parents=True, exist_ok=True)
        p = output_dir / "gene_program_table.csv"
        out_df.to_csv(p, index=False)
        log.info("  saved → %s", p)
        return out_df

    # ------------------------------------------------------------------
    # run_all
    # ------------------------------------------------------------------

    def run_all(
        self,
        output_dir: str | Path = "figures/svd_eval",
        cv: int = 5,
        n_permutations: int = 200,
    ) -> None:
        """Run every evaluation and save all outputs to *output_dir*."""
        out = Path(output_dir)
        log.info(
            "=== SVDEvaluator running for mutation '%s' → %s ===", self.mutation, out
        )

        self.plot_weighted_scree(out)
        self.plot_component_importance(out)
        self.plot_gene_loading_heatmap(out)
        self.plot_top_genes_per_component(out)
        self.plot_latent_scatter(out)
        self.plot_separation_violins(out)
        self.plot_component_label_correlation(out)
        self.plot_roc_ablation(out, cv=cv)
        self.plot_permutation_test(out, n_permutations=n_permutations, cv=cv)
        self.plot_gene_biplot(out)
        self.plot_shap(out)
        self.plot_umap(out)
        self.plot_calibration(out, cv=cv)
        self.plot_component_crosscorr(out)
        self.export_gene_program_table(out)

        log.info("=== SVDEvaluator complete ===")
