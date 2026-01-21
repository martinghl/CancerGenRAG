#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compute_metrics_from_ranks.py

Compute R2 / MSE / MAE per model by aligning predicted scores from ranks_{method}.csv
with labels from an aligned-priority CSV (e.g., LC_GWAS_Priority_Scores_aligned.csv).
This avoids length-mismatch issues by joining on IDs.

Usage
-----
python compute_metrics_from_ranks.py \
  --methods "RF,DT,LR,NN,SVR,Ridge,Lasso" \
  --ranks-dir outputs/sup \
  --aligned-priority outputs/LC_GWAS_Priority_Scores_aligned.csv \
  --ranks-id-col gene_variant_id \
  --aligned-id-col hm_rsID \
  --out-dir outputs/metrics

Notes
-----
- ranks_{method}.csv must include columns: gene_variant_id, score (rank optional).
- aligned priority CSV must include: <aligned-id-col>, priority_score
- Metrics are computed on the intersection after inner-join by ID.
"""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import json

from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

ALIAS_MAP = {
    "RF": ["RF","random_forest","random-forest","RandomForest","rf","config_rf","randomforest"],
    "DT": ["DT","decision_tree","decision-tree","DecisionTree","dt","config_dt","tree"],
    "LR": ["LR","linear_regression","LinearRegression","lr","config_lr","logistic_regression","LogisticRegression","logreg"],
    "NN": ["NN","neural_network","neural-network","NeuralNetwork","mlp","MLP","nn","config_nn"],
    "SVR": ["SVR","svr","support_vector_regression","svm_regression","svm","config_svr"],
    "Ridge": ["Ridge","ridge","config_ridge"],
    "Lasso": ["Lasso","lasso","config_lasso"],
}

def _find_rank_file(ranks_dir: Path, alias: str) -> Path:
    aliases = ALIAS_MAP.get(alias, []) + [alias]
    # direct
    for token in aliases:
        p = ranks_dir / f"ranks_{token}.csv"
        if p.exists():
            return p
        p2 = ranks_dir / token / "ranks.csv"
        if p2.exists():
            return p2
    # fuzzy
    for p in sorted(ranks_dir.glob("ranks_*.csv")):
        if p.stem.replace("ranks_","").lower() in [a.lower() for a in aliases]:
            return p
    raise FileNotFoundError(f"ranks for '{alias}' not found in {ranks_dir}")

def parse_args():
    ap = argparse.ArgumentParser(description="Compute R2/MSE/MAE per model by ID-aligned join.")
    ap.add_argument("--methods", required=True, help="Comma-separated aliases e.g. 'RF,DT,LR,...'")
    ap.add_argument("--ranks-dir", required=True)
    ap.add_argument("--aligned-priority", required=True)
    ap.add_argument("--ranks-id-col", default="gene_variant_id")
    ap.add_argument("--aligned-id-col", default="hm_rsID")
    ap.add_argument("--out-dir", default="outputs/metrics")
    return ap.parse_args()

def main():
    args = parse_args()
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    prio = pd.read_csv(args.aligned_priority)
    if args.aligned_id_col not in prio.columns:
        raise ValueError(f"aligned priority missing column: {args.aligned_id_col}")
    if "priority_score" not in prio.columns:
        raise ValueError("aligned priority missing column: priority_score")
    prio["_id_norm"] = prio[args.aligned_id_col].astype(str).str.strip().str.lower()
    prio = prio[["_id_norm", "priority_score"]].dropna(subset=["_id_norm"])

    summary_rows = []

    for alias in [m.strip() for m in args.methods.split(",") if m.strip()]:
        rf = _find_rank_file(Path(args.ranks_dir), alias)
        df = pd.read_csv(rf)
        if args.ranks_id_col not in df.columns or "score" not in df.columns:
            raise ValueError(f"{rf} must include '{args.ranks_id_col}' and 'score'")

        df["_id_norm"] = df[args.ranks_id_col].astype(str).str.strip().str.lower()
        df2 = df.merge(prio, on="_id_norm", how="inner", suffixes=("_pred","_true"))
        n = len(df2)

        metrics_path = out_dir / f"metrics_{alias}.json"
        if n >= 2:
            y_true = pd.to_numeric(df2["priority_score"], errors="coerce").to_numpy()
            y_pred = pd.to_numeric(df2["score"], errors="coerce").to_numpy()
            # drop NaNs
            msk = ~np.isnan(y_true) & ~np.isnan(y_pred)
            y_true = y_true[msk]; y_pred = y_pred[msk]
            if len(y_true) >= 2:
                r2 = float(r2_score(y_true, y_pred))
                mse = float(mean_squared_error(y_true, y_pred))
                mae = float(mean_absolute_error(y_true, y_pred))
                metrics = {"method": alias, "n_overlap": int(len(y_true)), "r2": r2, "mse": mse, "mae": mae}
            else:
                metrics = {"method": alias, "n_overlap": int(len(y_true)), "note": "not enough overlap after NaN drop"}
        else:
            metrics = {"method": alias, "n_overlap": n, "note": "not enough overlap to compute metrics"}

        with open(metrics_path, "w") as f:
            json.dump(metrics, f, indent=2)
        summary_rows.append(metrics)

        print(f"[{alias}] overlap={metrics.get('n_overlap')}  "
              f"r2={metrics.get('r2')}  mse={metrics.get('mse')}  mae={metrics.get('mae')}")

    # write summary CSV
    pd.DataFrame(summary_rows).to_csv(out_dir / "metrics_summary_from_ranks.csv", index=False)
    print(f"[ok] wrote summary: {out_dir / 'metrics_summary_from_ranks.csv'}")

if __name__ == "__main__":
    main()
