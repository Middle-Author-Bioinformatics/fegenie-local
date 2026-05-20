#!/usr/bin/env python3
"""
fegenie_report.py

Generate a self-contained, interactive HTML report from a FeGenie heatmap-format
CSV file. The report contains:

  1. A dot plot mirroring the heatmap-format CSV (genomes on x-axis, categories
     on y-axis, dot size and color encode the value).
  2. A clustered heatmap with dendrograms on both axes (rows = categories,
     columns = genomes), reflecting the same numeric values (integers or
     fractions) from the CSV.

The output HTML is fully self-contained (Plotly bundled inline) so it can be
downloaded and opened in any modern browser without a network connection.

Usage:
    python fegenie_report.py FeGenie-heatmap-data.csv [-o report.html]

Dependencies:
    pandas, numpy, scipy, plotly
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_heatmap_csv(csv_path: Path) -> pd.DataFrame:
    """Load a FeGenie heatmap CSV.

    Expected layout:
        first column = category labels (header is typically 'X')
        remaining columns = genomes (one per column)
        cell values    = counts / fractions (numeric)
    """
    df = pd.read_csv(csv_path, index_col=0)
    df.index.name = "Category"
    df.columns.name = "Genome"

    # Coerce to numeric; non-numeric becomes NaN -> 0 for plotting safety,
    # but we keep NaN-aware logic where it matters.
    df = df.apply(pd.to_numeric, errors="coerce")
    return df


# ---------------------------------------------------------------------------
# Dot plot
# ---------------------------------------------------------------------------
def build_dotplot(df: pd.DataFrame) -> go.Figure:
    """Dot plot: x = genomes, y = categories. Size + color encode value."""
    # Long-format for scatter
    long_df = (
        df.reset_index()
        .melt(id_vars="Category", var_name="Genome", value_name="Value")
    )

    # Keep zeros visible as faint markers so the grid reads cleanly,
    # but scale size from the non-zero distribution.
    vmax = float(np.nanmax(df.values)) if df.size else 1.0
    vmax = vmax if vmax > 0 else 1.0

    sizes = long_df["Value"].fillna(0).clip(lower=0)
    # Min marker size 4 (for zero), scaled up to 36 at vmax.
    marker_sizes = 4 + (sizes / vmax) * 32

    hover = [
        f"<b>{cat}</b><br>Genome: {gen}<br>Value: {val:g}"
        for cat, gen, val in zip(
            long_df["Category"], long_df["Genome"], long_df["Value"].fillna(0)
        )
    ]

    fig = go.Figure(
        go.Scatter(
            x=long_df["Genome"],
            y=long_df["Category"],
            mode="markers",
            marker=dict(
                size=marker_sizes,
                color=long_df["Value"].fillna(0),
                colorscale="Viridis",
                cmin=0,
                cmax=vmax,
                colorbar=dict(title="Value", thickness=14, len=0.75),
                line=dict(width=0.5, color="rgba(40,40,40,0.6)"),
            ),
            hovertemplate="%{text}<extra></extra>",
            text=hover,
        )
    )

    fig.update_layout(
        title="Dot plot (mirrors heatmap CSV values)",
        xaxis=dict(title="Genome", tickangle=-30, automargin=True),
        yaxis=dict(
            title="Category",
            automargin=True,
            categoryorder="array",
            categoryarray=list(df.index)[::-1],  # first row at top
        ),
        height=max(420, 38 * len(df.index) + 160),
        margin=dict(l=120, r=40, t=70, b=120),
        plot_bgcolor="#fafafa",
    )
    fig.update_xaxes(showgrid=True, gridcolor="#ececec")
    fig.update_yaxes(showgrid=True, gridcolor="#ececec")
    return fig


# ---------------------------------------------------------------------------
# Dendrogram + heatmap
# ---------------------------------------------------------------------------
def _cluster_transform(matrix: np.ndarray, mode: str) -> np.ndarray:
    """Transform values used only for clustering, not heatmap display.

    The raw FeGenie matrix is often sparse and count-dominated; a few abundant
    categories can otherwise control the genome dendrogram. These transforms
    make the clustering reflect profile similarity more than absolute magnitude.
    """
    data = np.nan_to_num(matrix.astype(float), nan=0.0)

    if mode == "none":
        return data

    if mode == "log1p":
        return np.log1p(np.clip(data, a_min=0, a_max=None))

    if mode == "presence_absence":
        return (data > 0).astype(float)

    if mode == "column_fraction":
        col_sums = data.sum(axis=0, keepdims=True)
        return np.divide(data, col_sums, out=np.zeros_like(data), where=col_sums != 0)

    if mode == "column_fraction_row_zscore":
        col_sums = data.sum(axis=0, keepdims=True)
        frac = np.divide(data, col_sums, out=np.zeros_like(data), where=col_sums != 0)
        row_mean = frac.mean(axis=1, keepdims=True)
        row_sd = frac.std(axis=1, keepdims=True)
        return np.divide(frac - row_mean, row_sd, out=np.zeros_like(frac), where=row_sd != 0)

    raise ValueError(f"Unknown cluster transform: {mode}")


def _safe_linkage(
    matrix: np.ndarray,
    axis: int,
    metric: str = "correlation",
    method: str = "average",
) -> np.ndarray | None:
    """Compute linkage along an axis. Returns None if not enough data."""
    data = matrix if axis == 0 else matrix.T
    n = data.shape[0]
    if n < 2:
        return None

    # Replace NaNs with 0 for distance computation.
    data = np.nan_to_num(data, nan=0.0)

    try:
        dists = pdist(data, metric=metric)
        # Correlation/cosine distances can become NaN for constant all-zero rows.
        # Treat those cases as maximally uninformative but finite so linkage works.
        dists = np.nan_to_num(dists, nan=1.0, posinf=1.0, neginf=0.0)
        if not np.any(dists):
            # Add tiny jitter so linkage produces a stable order.
            dists = dists + 1e-9
        return linkage(dists, method=method)
    except Exception:
        return None

def _dendrogram_traces(Z: np.ndarray, orientation: str, leaf_count: int):
    """Build line traces for a dendrogram from a SciPy linkage matrix.

    orientation: 'top' (columns) or 'left' (rows).
    Leaves are placed at integer positions 0..leaf_count-1 in the order given
    by scipy.cluster.hierarchy.leaves_list(Z).
    Returns (traces, max_height).
    """
    leaves = leaves_list(Z)
    # position[node_id] = x-or-y coordinate of that node along the leaf axis
    n = leaf_count
    position = np.zeros(2 * n - 1)
    position[:n] = np.argsort(leaves)  # leaf i is at slot argsort(leaves)[i]
    # Reassign so leaf at index k of `leaves` lives at coordinate k:
    # position[leaf_id] = k where leaves[k] == leaf_id
    inv = np.empty_like(leaves)
    inv[leaves] = np.arange(n)
    position[:n] = inv

    height = np.zeros(2 * n - 1)
    traces = []

    for i, (a, b, dist, _count) in enumerate(Z):
        a, b = int(a), int(b)
        node_id = n + i
        xa, xb = position[a], position[b]
        ya, yb = height[a], height[b]
        position[node_id] = (xa + xb) / 2.0
        height[node_id] = dist

        if orientation == "top":
            # Horizontal axis = leaf positions, vertical axis = distance
            xs = [xa, xa, xb, xb]
            ys = [ya, dist, dist, yb]
        else:  # 'left'
            # Vertical axis = leaf positions, horizontal axis = distance
            xs = [ya, dist, dist, yb]
            ys = [xa, xa, xb, xb]

        traces.append(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines",
                line=dict(color="#444", width=1),
                hoverinfo="skip",
                showlegend=False,
            )
        )

    max_height = float(height.max()) if len(height) else 1.0
    if max_height == 0:
        max_height = 1.0
    return traces, max_height, leaves


def build_clustered_heatmap(
    df: pd.DataFrame,
    cluster_transform: str = "log1p",
    cluster_metric: str = "correlation",
    cluster_method: str = "average",
) -> go.Figure:
    """Heatmap with hierarchical-clustering dendrograms on top and left."""
    matrix = df.values.astype(float)
    cluster_matrix = _cluster_transform(matrix, cluster_transform)
    row_labels = list(df.index)
    col_labels = list(df.columns)

    Z_rows = _safe_linkage(cluster_matrix, axis=0, metric=cluster_metric, method=cluster_method)
    Z_cols = _safe_linkage(cluster_matrix, axis=1, metric=cluster_metric, method=cluster_method)

    # Reorder rows / columns according to dendrogram leaves
    if Z_rows is not None:
        row_order = leaves_list(Z_rows)
    else:
        row_order = np.arange(len(row_labels))
    if Z_cols is not None:
        col_order = leaves_list(Z_cols)
    else:
        col_order = np.arange(len(col_labels))

    ordered_rows = [row_labels[i] for i in row_order]
    ordered_cols = [col_labels[i] for i in col_order]
    ordered_matrix = matrix[np.ix_(row_order, col_order)]

    # Layout: 2x2 grid. Top-left empty, top-right = column dendrogram,
    # bottom-left = row dendrogram, bottom-right = heatmap.
    fig = make_subplots(
        rows=2,
        cols=2,
        column_widths=[0.12, 0.88],
        row_heights=[0.18, 0.82],
        horizontal_spacing=0.01,
        vertical_spacing=0.01,
        shared_xaxes=False,
        shared_yaxes=False,
    )

    # --- Column dendrogram (top-right) ---
    if Z_cols is not None:
        col_traces, col_max_h, _ = _dendrogram_traces(
            Z_cols, orientation="top", leaf_count=len(col_labels)
        )
        for tr in col_traces:
            fig.add_trace(tr, row=1, col=2)
        fig.update_xaxes(
            range=[-0.5, len(col_labels) - 0.5],
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            row=1,
            col=2,
        )
        fig.update_yaxes(
            range=[0, col_max_h * 1.05],
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            row=1,
            col=2,
        )

    # --- Row dendrogram (bottom-left) ---
    if Z_rows is not None:
        row_traces, row_max_h, _ = _dendrogram_traces(
            Z_rows, orientation="left", leaf_count=len(row_labels)
        )
        for tr in row_traces:
            fig.add_trace(tr, row=2, col=1)
        # Distance grows from right (closest to heatmap) to left, so reverse x.
        fig.update_xaxes(
            range=[row_max_h * 1.05, 0],
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            row=2,
            col=1,
        )
        fig.update_yaxes(
            range=[-0.5, len(row_labels) - 0.5],
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            row=2,
            col=1,
        )

    # --- Heatmap (bottom-right) ---
    # Show numeric annotations inside cells so fractional values are readable.
    text_vals = np.where(
        np.isnan(ordered_matrix),
        "",
        np.vectorize(lambda v: f"{v:g}")(ordered_matrix),
    )

    fig.add_trace(
        go.Heatmap(
            z=ordered_matrix,
            x=ordered_cols,
            y=ordered_rows,
            colorscale="Viridis",
            colorbar=dict(
                title=dict(text="Value", side="top"),
                orientation="h",
                thickness=14,
                len=0.45,
                x=0.58,
                xanchor="center",
                y=-0.32,
                yanchor="top",
            ),
            hovertemplate=(
                "<b>%{y}</b><br>Genome: %{x}<br>Value: %{z:g}<extra></extra>"
            ),
            text=text_vals,
            texttemplate="%{text}",
            textfont=dict(size=11, color="white"),
            zmin=float(np.nanmin(ordered_matrix)) if ordered_matrix.size else 0,
            zmax=float(np.nanmax(ordered_matrix)) if ordered_matrix.size else 1,
        ),
        row=2,
        col=2,
    )
    fig.update_xaxes(
        tickangle=-30,
        automargin=True,
        showgrid=False,
        zeroline=False,
        row=2,
        col=2,
    )
    fig.update_yaxes(
        automargin=True,
        showgrid=False,
        zeroline=False,
        side="right",
        tickfont=dict(size=12),
        row=2,
        col=2,
    )

    fig.update_layout(
        title="Clustered heatmap with dendrograms",
        height=max(520, 42 * len(row_labels) + 220),
        margin=dict(l=110, r=260, t=70, b=230),
        plot_bgcolor="white",
        showlegend=False,
    )
    return fig


# ---------------------------------------------------------------------------
# HTML assembly
# ---------------------------------------------------------------------------
HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>FeGenie report &mdash; {source}</title>
<style>
  :root {{
    --fg: #1f2937;
    --muted: #6b7280;
    --bg: #ffffff;
    --panel: #f8fafc;
    --border: #e5e7eb;
    --accent: #2563eb;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    color: var(--fg);
    background: var(--bg);
    margin: 0;
    padding: 32px 24px 64px;
    line-height: 1.5;
  }}
  .container {{ max-width: 1200px; margin: 0 auto; }}
  header h1 {{
    font-size: 1.6rem;
    margin: 0 0 4px;
  }}
  header p {{
    color: var(--muted);
    margin: 0 0 24px;
    font-size: 0.95rem;
  }}
  section {{
    background: var(--panel);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 18px 20px;
    margin-bottom: 28px;
  }}
  section h2 {{
    margin: 0 0 4px;
    font-size: 1.15rem;
  }}
  section p.desc {{
    color: var(--muted);
    margin: 0 0 14px;
    font-size: 0.9rem;
  }}
  footer {{
    color: var(--muted);
    font-size: 0.8rem;
    text-align: center;
    margin-top: 32px;
  }}
  details summary {{
    cursor: pointer;
    color: var(--accent);
    font-size: 0.9rem;
  }}
  table {{
    border-collapse: collapse;
    margin-top: 10px;
    font-size: 0.85rem;
  }}
  th, td {{
    border: 1px solid var(--border);
    padding: 4px 8px;
    text-align: right;
  }}
  th:first-child, td:first-child {{ text-align: left; }}
  thead th {{ background: #eef2ff; }}
</style>
</head>
<body>
<div class="container">
  <header>
    <h1>FeGenie heatmap report</h1>
    <p>Source: <code>{source}</code> &middot; {n_rows} categories &times; {n_cols} genomes</p>
  </header>

  <section>
    <h2>Dot plot</h2>
    <p class="desc">Mirrors the heatmap-format CSV. Marker size and color both encode the cell value; zeros are kept visible as small dots.</p>
    {dotplot_div}
  </section>

  <section>
    <h2>Clustered heatmap with dendrograms</h2>
    <p class="desc">Rows (categories) and columns (genomes) are reordered by hierarchical clustering. By default, clustering uses log1p-transformed values with average linkage and correlation distance, while the heatmap itself still displays the raw input values. Category labels are placed on the right to avoid overlap with the row dendrogram.</p>
    {heatmap_div}
  </section>

  <section>
    <details>
      <summary>Show raw data table</summary>
      {table_html}
    </details>
  </section>

  <footer>
    Generated by <code>Middle Author Bioinformatics</code>. Open this file in any modern browser &mdash; no internet required.
  </footer>
</div>
</body>
</html>
"""


def build_html_report(
    df: pd.DataFrame,
    source_name: str,
    cluster_transform: str = "log1p",
    cluster_metric: str = "correlation",
    cluster_method: str = "average",
) -> str:
    dot_fig = build_dotplot(df)
    heat_fig = build_clustered_heatmap(
        df,
        cluster_transform=cluster_transform,
        cluster_metric=cluster_metric,
        cluster_method=cluster_method,
    )

    # First figure includes the Plotly JS bundle inline; subsequent figures
    # reuse it so the file stays self-contained but isn't duplicated.
    dot_div = dot_fig.to_html(
        include_plotlyjs="inline", full_html=False, div_id="dotplot"
    )
    heat_div = heat_fig.to_html(
        include_plotlyjs=False, full_html=False, div_id="heatmap"
    )

    table_html = df.to_html(
        classes="rawdata", border=0, float_format=lambda v: f"{v:g}"
    )

    return HTML_TEMPLATE.format(
        source=source_name,
        n_rows=df.shape[0],
        n_cols=df.shape[1],
        dotplot_div=dot_div,
        heatmap_div=heat_div,
        table_html=table_html,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate an HTML report (dot plot + clustered heatmap) "
        "from a FeGenie heatmap-format CSV file."
    )
    parser.add_argument("csv", type=Path, help="Input heatmap-format CSV")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output HTML path (default: <input>.report.html)",
    )
    parser.add_argument(
        "--cluster-transform",
        choices=[
            "none",
            "log1p",
            "presence_absence",
            "column_fraction",
            "column_fraction_row_zscore",
        ],
        default="log1p",
        help=(
            "Transform used only for dendrogram clustering; the heatmap still "
            "shows raw values. Default: log1p."
        ),
    )
    parser.add_argument(
        "--cluster-metric",
        default="correlation",
        help="Distance metric passed to scipy.spatial.distance.pdist. Default: correlation.",
    )
    parser.add_argument(
        "--cluster-method",
        default="average",
        help="Linkage method passed to scipy.cluster.hierarchy.linkage. Default: average.",
    )
    args = parser.parse_args(argv)

    if not args.csv.is_file():
        print(f"Error: {args.csv} not found", file=sys.stderr)
        return 2

    out_path = args.output or args.csv.with_suffix(".report.html")

    df = load_heatmap_csv(args.csv)
    if df.empty:
        print("Error: input CSV produced an empty data frame", file=sys.stderr)
        return 2

    html = build_html_report(
        df,
        source_name=args.csv.name,
        cluster_transform=args.cluster_transform,
        cluster_metric=args.cluster_metric,
        cluster_method=args.cluster_method,
    )
    out_path.write_text(html, encoding="utf-8")
    print(f"Wrote {out_path}  ({df.shape[0]} categories x {df.shape[1]} genomes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
