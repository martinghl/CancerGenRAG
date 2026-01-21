#!/usr/bin/env python3
"""
plot_unsup_comparison.py
------------------------
Generate large, clear comparison plots between unsupervised method rankings and a
GWAS-based priority score.

USAGE (example):
python plot_unsup_comparison.py \
  --unsup-dir outputs/unsup \
  --methods GEO,MC-3,THURSTONE,RRA,BARD,CEMC_SPEARMAN,CEMC_KENDALL \
  --priority-csv /path/to/LC_GWAS_Priority_Scores.csv \
  --priority-id-col variant_id \
  --priority-col priority_score \
  --topn 12 \
  --ext png

Notes
-----
- The script tries hard to auto-detect the correct ID and rank columns in each
  method's aggregated CSV. You can keep the CLI exactly as above.
- For matching, it first tries "rsID" extraction (e.g., finds "rs12345" inside
  a string like "PGS000111|rs12345"). If that fails, it falls back to raw ID
  equality.
- Figures are large and high-DPI by default. Layout is adaptive and uses
  tight bounding boxes to avoid squashed subplots.
"""

import argparse
import os
import re
import sys
import glob
from math import ceil
from typing import Optional, Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


# ----------------------------- CLI & utilities ----------------------------- #

def parse_args():
    p = argparse.ArgumentParser(description="Plot unsupervised method comparisons vs a priority score.")
    p.add_argument("--unsup-dir", required=True, help="Directory containing unsupervised method outputs.")
    p.add_argument("--methods", required=True,
                   help="Comma-separated method names, e.g., GEO,MC-3,THURSTONE,RRA,BARD,CEMC_SPEARMAN,CEMC_KENDALL")
    p.add_argument("--priority-csv", required=True, help="CSV with priority scores.")
    p.add_argument("--priority-id-col", required=True, help="Column in priority CSV with variant ID (e.g., variant_id).")
    p.add_argument("--priority-col", required=True, help="Column in priority CSV with numeric priority scores.")
    p.add_argument("--topn", type=int, default=12, help="Top-N items (lowest ranks) per method to plot.")
    p.add_argument("--ext", default="png", choices=["png", "pdf", "svg"], help="Output image format.")
    # Optional (sane defaults so your current call doesn't need to change)
    p.add_argument("--dpi", type=int, default=600, help="Figure DPI for saved images (default: 600).")
    p.add_argument("--outfile", default=None, help="Optional output filename (without extension). Default: auto.")
    return p.parse_args()


def norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", s.lower())


def extract_first_rsid(val: str) -> Optional[str]:
    """Return the first 'rs1234' found in val (case-insensitive), else None."""
    if val is None:
        return None
    m = re.search(r"(rs\d+)", str(val), flags=re.IGNORECASE)
    return m.group(1).lower() if m else None


def find_aggregated_csv_for_method(unsup_dir: str, method: str) -> Optional[str]:
    """
    Try a few common locations/filenames, then fall back to a recursive scan.
    We prioritize files that appear under a folder matching the method name and
    that contain 'aggreg'/'consensus' in the filename.
    """
    m_upper = method.upper()
    m_lower = method.lower()
    m_norm = norm(method)

    candidates = [
        os.path.join(unsup_dir, m_upper, "aggregated.csv"),
        os.path.join(unsup_dir, m_upper, "aggregated_ranks.csv"),
        os.path.join(unsup_dir, m_upper, "consensus.csv"),
        os.path.join(unsup_dir, m_upper, "aggregated_results.csv"),
        os.path.join(unsup_dir, m_upper, "output.csv"),
        os.path.join(unsup_dir, m_lower, "aggregated.csv"),
        os.path.join(unsup_dir, m_lower, "consensus.csv"),
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c

    # Fallback: recursive scan in unsup_dir
    all_csvs = glob.glob(os.path.join(unsup_dir, "**", "*.csv"), recursive=True)
    scored: List[Tuple[int, str]] = []
    for path in all_csvs:
        base = os.path.basename(path).lower()
        path_norm = norm(path)
        score = 0
        if m_norm in path_norm:
            score += 2
        if "aggreg" in base or "consensus" in base:
            score += 2
        # small boost if parent folder looks like the method
        parent = os.path.basename(os.path.dirname(path))
        if norm(parent) == m_norm:
            score += 3
        if score > 0:
            scored.append((score, path))

    if not scored:
        return None
    # Pick the highest-scored candidate
    scored.sort(key=lambda x: (-x[0], x[1]))
    return scored[0][1]


def pick_id_and_rank_cols(df: pd.DataFrame) -> Tuple[Optional[str], Optional[str]]:
    """
    Try to infer the ID and RANK columns from a method aggregated CSV.
    """
    cols = [c.lower() for c in df.columns]

    # Rank column candidates (ordered by priority)
    rank_cands = ["aggregated_rank", "consensus_rank", "rank", "rra_rank", "borda_rank", "rank_score", "score", "value"]
    rank_col = None
    for rc in rank_cands:
        if rc in cols:
            rank_col = df.columns[cols.index(rc)]
            break

    # ID column candidates
    id_cands = ["gene_variant_id", "variant_id", "hm_rsid", "rsid", "id", "variant"]
    id_col = None
    for ic in id_cands:
        if ic in cols:
            id_col = df.columns[cols.index(ic)]
            break

    return id_col, rank_col


def safe_float_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def merge_on_rsid_then_fallback(method_df: pd.DataFrame,
                                priority_df: pd.DataFrame,
                                method_id_col: str,
                                priority_id_col: str) -> pd.DataFrame:
    m = method_df.copy()
    p = priority_df.copy()

    m["__rsid__"] = m[method_id_col].astype(str).apply(extract_first_rsid)
    p["__rsid__"] = p[priority_id_col].astype(str).apply(extract_first_rsid)

    # Prefer rsid join if it gives reasonable matches
    merged = m.merge(p, on="__rsid__", how="inner", suffixes=("_m", "_p"))
    if len(merged) >= 2:
        return merged

    # Fallback: direct ID equality
    merged = method_df.merge(priority_df, left_on=method_id_col, right_on=priority_id_col, how="inner",
                             suffixes=("_m", "_p"))
    return merged


# --------------------------------- Plotting -------------------------------- #

def make_plots(per_method_frames: List[Tuple[str, pd.DataFrame, str, str]],
               topn: int,
               priority_col: str,
               ext: str,
               dpi: int,
               outfile: Optional[str] = None):
    # Filter out empty frames
    usable = [(name, df, rank_col, id_col) for (name, df, rank_col, id_col) in per_method_frames if df is not None and len(df) > 0]
    if not usable:
        print("[ERROR] No usable data to plot. Aborting.")
        sys.exit(2)

    n = len(usable)
    ncols = 3 if n >= 3 else n
    nrows = int(ceil(n / ncols))
    fig_w = max(5 * ncols, 10)  # ensure wide enough
    fig_h = max(4.2 * nrows, 5)
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), constrained_layout=True)
    if nrows * ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = np.array([axes])
    elif ncols == 1:
        axes = axes.reshape(-1, 1)

    # Common styling
    plt.rcParams.update({
        "font.size": 12,
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10
    })

    for idx, (name, df, rank_col, id_col) in enumerate(usable):
        r = idx // ncols
        c = idx % ncols
        ax = axes[r, c]

        # Keep only needed cols, parse types
        df_plot = df[[rank_col, priority_col]].copy()
        df_plot[rank_col] = safe_float_series(df_plot[rank_col])
        df_plot[priority_col] = safe_float_series(df_plot[priority_col])
        df_plot = df_plot.dropna(subset=[rank_col, priority_col])
        df_plot = df_plot.sort_values(by=rank_col, ascending=True).head(topn)

        if len(df_plot) == 0:
            ax.text(0.5, 0.5, "No data", ha="center", va="center")
            ax.set_axis_off()
            continue

        x = df_plot[rank_col].to_numpy()
        y = df_plot[priority_col].to_numpy()

        # Spearman correlation for rank-vs-priority score
        try:
            rho, pval = spearmanr(x, y, nan_policy="omit")
        except Exception:
            rho, pval = np.nan, np.nan

        # Scatter
        ax.scatter(x, y, s=26, alpha=0.85, linewidths=0.5, edgecolors="none")

        # Simple (unweighted) linear trendline for visual guidance
        if len(df_plot) >= 2:
            try:
                z = np.polyfit(x, y, 1)
                xp = np.linspace(x.min(), x.max(), 200)
                yp = z[0] * xp + z[1]
                ax.plot(xp, yp, linewidth=1.4, linestyle="-", alpha=0.9)
            except Exception:
                pass

        ax.set_title(f"{name} (ρ={rho:.2f}, n={len(df_plot)})")
        ax.set_xlabel("Aggregated Rank (lower is better)")
        ax.set_ylabel(priority_col)

        ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.4)
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

    # Turn off any unused subplots
    total_axes = nrows * ncols
    for k in range(n, total_axes):
        r = k // ncols
        c = k % ncols
        axes[r, c].set_axis_off()

    # Save
    base = outfile if outfile else f"unsup_comparison_top{topn}"
    outpath = f"{base}.{ext}"
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    print(f"[OK] Saved: {outpath}")


# --------------------------------- Main ----------------------------------- #

def main():
    args = parse_args()

    # Read priority CSV
    if not os.path.isfile(args.priority_csv):
        print(f"[ERROR] priority-csv not found: {args.priority_csv}")
        sys.exit(1)

    priority_df = pd.read_csv(args.priority_csv)
    if args.priority_id_col not in priority_df.columns:
        print(f"[ERROR] priority-id-col '{args.priority_id_col}' not found in {args.priority_csv}. "
              f"Columns: {list(priority_df.columns)}")
        sys.exit(1)
    if args.priority_col not in priority_df.columns:
        print(f"[ERROR] priority-col '{args.priority_col}' not found in {args.priority_csv}. "
              f"Columns: {list(priority_df.columns)}")
        sys.exit(1)

    # Normalize methods list
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    if not methods:
        print("[ERROR] No methods provided.")
        sys.exit(1)

    per_method_frames = []  # list of (name, merged_df, rank_col, id_col)
    for method in methods:
        csv_path = find_aggregated_csv_for_method(args.unsup_dir, method)
        if not csv_path:
            print(f"[WARN] No aggregated file found for method '{method}' under {args.unsup_dir}")
            per_method_frames.append((method, None, None, None))
            continue

        try:
            mdf = pd.read_csv(csv_path)
        except Exception as e:
            print(f"[WARN] Could not read CSV for method '{method}': {csv_path}\n  -> {e}")
            per_method_frames.append((method, None, None, None))
            continue

        id_col, rank_col = pick_id_and_rank_cols(mdf)
        if id_col is None or rank_col is None:
            print(f"[WARN] Could not infer id/rank columns for '{method}' file: {csv_path}")
            per_method_frames.append((method, None, None, None))
            continue

        merged = merge_on_rsid_then_fallback(
            mdf[[id_col, rank_col]].copy(),  # only the needed cols
            priority_df[[args.priority_id_col, args.priority_col]].copy(),
            id_col, args.priority_id_col
        )

        if len(merged) == 0:
            print(f"[WARN] Merged 0 rows for method '{method}'. Check ID alignment. ({csv_path})")
            per_method_frames.append((method, None, None, None))
            continue

        per_method_frames.append((method, merged, rank_col, id_col))

    make_plots(per_method_frames, args.topn, args.priority_col, args.ext, args.dpi, outfile=args.outfile)


if __name__ == "__main__":
    main()
