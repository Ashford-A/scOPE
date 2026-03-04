"""Visualisation utilities for scOPE embeddings and mutation probabilities."""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData
from matplotlib.axes import Axes
from matplotlib.colors import Normalize
from matplotlib.figure import Figure

from scope.utils.logging import get_logger

log = get_logger(__name__)


# ---------------------------------------------------------------------------
# 2-D embedding helpers
# ---------------------------------------------------------------------------


def compute_umap(
    adata: AnnData,
    obsm_key: str = "X_svd",
    n_neighbors: int = 15,
    min_dist: float = 0.3,
    n_components: int = 2,
    random_state: int = 42,
    umap_key: str = "X_umap",
) -> AnnData:
    """Compute a UMAP embedding from the latent representation.

    Requires ``umap-learn`` (``pip install scope-bio[full]``).

    Parameters
    ----------
    adata:
        AnnData with latent embedding in ``adata.obsm[obsm_key]``.
    obsm_key:
        Key to read the pre-computed latent coords from.
    n_neighbors, min_dist:
        Standard UMAP hyperparameters.
    n_components:
        Number of UMAP dimensions (2 for visualisation).
    random_state:
        Seed.
    umap_key:
        Key to write the UMAP coords into ``adata.obsm``.

    Returns
    -------
    AnnData with ``adata.obsm[umap_key]`` filled.
    """
    try:
        import umap
    except ImportError as e:
        raise ImportError(
            "umap-learn is required. Install it with: pip install scope-bio[full]"
        ) from e

    Z = adata.obsm[obsm_key]
    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=n_components,
        random_state=random_state,
    )
    coords = reducer.fit_transform(Z)
    adata = adata.copy()
    adata.obsm[umap_key] = coords.astype(np.float32)
    log.info("UMAP computed: %d cells → %dD.", adata.n_obs, n_components)
    return adata


def compute_tsne(
    adata: AnnData,
    obsm_key: str = "X_svd",
    perplexity: float = 30.0,
    n_components: int = 2,
    random_state: int = 42,
    tsne_key: str = "X_tsne",
) -> AnnData:
    """Compute a t-SNE embedding from the latent representation.

    Parameters
    ----------
    adata:
        AnnData with latent embedding.
    obsm_key:
        Key of the latent representation.
    perplexity:
        t-SNE perplexity.
    n_components:
        Number of output dimensions (typically 2).
    random_state:
        Seed.
    tsne_key:
        Key to write t-SNE coords into ``adata.obsm``.

    Returns
    -------
    AnnData with ``adata.obsm[tsne_key]`` filled.
    """
    from sklearn.manifold import TSNE

    Z = adata.obsm[obsm_key]
    tsne = TSNE(
        n_components=n_components,
        perplexity=perplexity,
        random_state=random_state,
    )
    coords = tsne.fit_transform(Z)
    adata = adata.copy()
    adata.obsm[tsne_key] = coords.astype(np.float32)
    log.info("t-SNE computed: %d cells → %dD.", adata.n_obs, n_components)
    return adata


# ---------------------------------------------------------------------------
# Scatter plots
# ---------------------------------------------------------------------------


def plot_embedding(
    adata: AnnData,
    color_key: Optional[str] = None,
    obsm_key: str = "X_umap",
    ax: Optional[Axes] = None,
    figsize: Tuple[float, float] = (6.0, 5.0),
    cmap: str = "viridis",
    point_size: float = 3.0,
    alpha: float = 0.7,
    title: Optional[str] = None,
    colorbar_label: Optional[str] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
) -> Tuple[Figure, Axes]:
    """Scatter plot of a 2-D embedding, coloured by a scalar or categorical key.

    Parameters
    ----------
    adata:
        AnnData with embedding in ``adata.obsm[obsm_key]``.
    color_key:
        Column in ``adata.obs`` used for colouring.  If ``None``, all points
        are coloured uniformly.
    obsm_key:
        Key of the 2-D embedding.
    ax:
        Matplotlib axes; created if ``None``.
    figsize:
        Figure size in inches.
    cmap:
        Colormap for continuous variables.
    point_size:
        Scatter point size.
    alpha:
        Point transparency.
    title:
        Axes title.
    colorbar_label:
        Label for the colourbar (continuous) or legend (categorical).
    vmin, vmax:
        Colourbar limits for continuous variables.

    Returns
    -------
    (fig, ax)
    """
    coords = adata.obsm[obsm_key]
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    if color_key is None:
        ax.scatter(coords[:, 0], coords[:, 1], s=point_size, alpha=alpha, rasterized=True)
    else:
        values = adata.obs[color_key].values
        if pd.api.types.is_numeric_dtype(values):
            sc = ax.scatter(
                coords[:, 0], coords[:, 1],
                c=values.astype(float),
                s=point_size,
                alpha=alpha,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                rasterized=True,
            )
            cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
            if colorbar_label:
                cb.set_label(colorbar_label)
        else:
            # Categorical
            categories = pd.Categorical(values)
            palette = plt.get_cmap("tab20").colors
            for i, cat in enumerate(categories.categories):
                mask = categories == cat
                ax.scatter(
                    coords[mask, 0], coords[mask, 1],
                    s=point_size,
                    alpha=alpha,
                    color=palette[i % len(palette)],
                    label=str(cat),
                    rasterized=True,
                )
            ax.legend(
                title=colorbar_label or color_key,
                markerscale=3,
                bbox_to_anchor=(1.02, 1),
                loc="upper left",
                frameon=False,
            )
    ax.set_xlabel(f"{obsm_key}-1")
    ax.set_ylabel(f"{obsm_key}-2")
    ax.set_title(title or (color_key or ""))
    ax.spines[["top", "right"]].set_visible(False)
    return fig, ax


def plot_mutation_probabilities(
    adata: AnnData,
    mutations: Optional[List[str]] = None,
    obsm_key: str = "X_umap",
    figsize_per_panel: Tuple[float, float] = (4.5, 4.0),
    cmap: str = "RdBu_r",
    vmin: float = 0.0,
    vmax: float = 1.0,
    point_size: float = 3.0,
    ncols: int = 3,
) -> Figure:
    """Grid of scatter panels showing per-cell mutation probabilities.

    Parameters
    ----------
    adata:
        Single-cell AnnData with mutation probabilities in ``adata.obs``
        (columns named ``mutation_prob_<gene>``).
    mutations:
        List of gene names (without the ``mutation_prob_`` prefix) to plot.
        If ``None``, all ``mutation_prob_*`` columns are used.
    obsm_key:
        Embedding key.
    figsize_per_panel:
        Size of each individual panel.
    cmap:
        Colormap.  ``"RdBu_r"`` works well for probabilities centred at 0.5.
    vmin, vmax:
        Probability axis limits.
    point_size:
        Point size.
    ncols:
        Number of columns in the grid.

    Returns
    -------
    matplotlib.figure.Figure
    """
    prob_cols = [c for c in adata.obs.columns if c.startswith("mutation_prob_")]
    if mutations is not None:
        prob_cols = [f"mutation_prob_{m}" for m in mutations if f"mutation_prob_{m}" in adata.obs.columns]
    if not prob_cols:
        raise ValueError("No mutation_prob_* columns found in adata.obs.")

    ncols = min(ncols, len(prob_cols))
    nrows = int(np.ceil(len(prob_cols) / ncols))
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows),
        squeeze=False,
    )
    coords = adata.obsm[obsm_key]
    for idx, col in enumerate(prob_cols):
        row, c = divmod(idx, ncols)
        ax = axes[row][c]
        values = adata.obs[col].values.astype(float)
        sc = ax.scatter(
            coords[:, 0], coords[:, 1],
            c=values, cmap=cmap, vmin=vmin, vmax=vmax,
            s=point_size, alpha=0.8, rasterized=True,
        )
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        gene = col.replace("mutation_prob_", "")
        ax.set_title(f"P({gene} mut)")
        ax.set_xlabel(f"{obsm_key}-1")
        ax.set_ylabel(f"{obsm_key}-2")
        ax.spines[["top", "right"]].set_visible(False)

    # Hide unused axes
    for idx in range(len(prob_cols), nrows * ncols):
        row, c = divmod(idx, ncols)
        axes[row][c].set_visible(False)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Scree / singular value plot
# ---------------------------------------------------------------------------


def plot_scree(
    scree_data: dict,
    ax: Optional[Axes] = None,
    figsize: Tuple[float, float] = (7.0, 4.0),
    max_components: Optional[int] = None,
) -> Tuple[Figure, Axes]:
    """Plot singular values and cumulative explained variance.

    Parameters
    ----------
    scree_data:
        Dict returned by :meth:`~scope.decomposition.svd.SVDDecomposition.scree_data`.
    ax:
        Axes to draw on.  If ``None``, a new figure is created.
    figsize:
        Figure size.
    max_components:
        Truncate plot at this component number.

    Returns
    -------
    (fig, ax)
    """
    comp = scree_data["component"]
    evr = scree_data["explained_variance_ratio"]
    cumevr = scree_data["cumulative_evr"]
    if max_components is not None:
        comp = comp[:max_components]
        evr = evr[:max_components]
        cumevr = cumevr[:max_components]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    ax2 = ax.twinx()
    ax.bar(comp, evr, color="#4393c3", alpha=0.8, label="Explained variance ratio")
    ax2.plot(comp, cumevr, color="#d6604d", marker="o", markersize=3, label="Cumulative EVR")
    ax2.axhline(0.9, ls="--", color="grey", lw=0.8, label="90% threshold")
    ax.set_xlabel("Component")
    ax.set_ylabel("Explained variance ratio", color="#4393c3")
    ax2.set_ylabel("Cumulative EVR", color="#d6604d")
    ax.set_title("Scree plot")
    ax.spines[["top"]].set_visible(False)
    fig.tight_layout()
    return fig, ax


# ---------------------------------------------------------------------------
# Heatmap of mutation probabilities across cell clusters
# ---------------------------------------------------------------------------


def plot_mutation_heatmap(
    adata: AnnData,
    cluster_key: str = "leiden",
    mutations: Optional[List[str]] = None,
    figsize: Tuple[float, float] = (8.0, 5.0),
    cmap: str = "YlOrRd",
    vmin: float = 0.0,
    vmax: float = 1.0,
) -> Tuple[Figure, Axes]:
    """Heatmap: mean mutation probability per cell cluster.

    Parameters
    ----------
    adata:
        Single-cell AnnData with ``adata.obs`` containing cluster labels and
        ``mutation_prob_*`` columns.
    cluster_key:
        Key in ``adata.obs`` for cluster assignments.
    mutations:
        List of mutation names; ``None`` → all available.
    figsize:
        Figure size.
    cmap:
        Colormap.
    vmin, vmax:
        Colour axis limits.

    Returns
    -------
    (fig, ax)
    """
    import seaborn as sns

    prob_cols = [c for c in adata.obs.columns if c.startswith("mutation_prob_")]
    if mutations is not None:
        prob_cols = [f"mutation_prob_{m}" for m in mutations if f"mutation_prob_{m}" in adata.obs.columns]

    df = adata.obs[[cluster_key] + prob_cols].copy()
    mean_prob = df.groupby(cluster_key)[prob_cols].mean()
    mean_prob.columns = [c.replace("mutation_prob_", "") for c in mean_prob.columns]

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        mean_prob.T,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        annot=True,
        fmt=".2f",
        ax=ax,
        linewidths=0.4,
        linecolor="white",
    )
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Mutation")
    ax.set_title("Mean mutation probability per cluster")
    fig.tight_layout()
    return fig, ax
