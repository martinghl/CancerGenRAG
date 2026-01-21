#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
train_eval_topk_end2end.py

End-to-end utility to:
  1) Retrain each method by calling your MAIN script (fresh run).
  2) Read the NEWLY written ranks_{method}.csv (must include: gene_variant_id, score [, rank]).
  3) If --aligned-priority is provided, ID-join with priority_score to compute R2/MSE/MAE per method (robust to length mismatch).
  4) Emit Top-K per method; if aligned priority provided, also emit Top-K (labeled-only / unlabeled-only) using used_for_fit flag.

This file requires NO patching inside your original main. It orchestrates the run and evaluation around it.

Usage
-----
python train_eval_topk_end2end.py \
  --main-script ./main_sup_old_patched_score_labels_v2.py \
  --config-glob "configs/config_*.json" \
  --methods "RF,DT,LR,NN,SVR,Ridge,Lasso" \
  --ranks-dir outputs/sup \
  --topk 50 \
  --out-dir outputs/end2end \
  --seeds 123 \
  --aligned-priority outputs/LC_GWAS_Priority_Scores_aligned.csv \
  --aligned-id-col hm_rsID

Outputs
-------
- <out-dir>/topk/top50_{method}.csv
- <out-dir>/topk/top50_{method}_labeled.csv (if --aligned-priority given)
- <out-dir>/topk/top50_{method}_unlabeled.csv (if --aligned-priority given)
- <out-dir>/metrics/metrics_{method}.json
- <out-dir>/metrics/metrics_summary.csv

Notes
-----
- Your main script must write ranks_{method}.csv for each method (with score column). If rank is absent, we will rank by descending score.
- Metrics are computed ONLY on IDs present in both ranks and aligned priority (inner join).
"""

import argparse
import os
import sys
import subprocess
from pathlib import Path
import json

import pandas as pd
import numpy as np
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

def _find_rank_file(ranks_dir: Path, method_alias: str) -> Path:
    aliases = ALIAS_MAP.get(method_alias, []) + [method_alias]
    for token in aliases:
        p1 = ranks_dir / f"ranks_{token}.csv"
        if p1.exists():
            return p1
        p2 = ranks_dir / token / "ranks.csv"
        if p2.exists():
            return p2
    # fuzzy match
    for p in sorted(ranks_dir.glob("ranks_*.csv")):
        if p.stem.replace("ranks_","").lower() in [a.lower() for a in aliases]:
            return p
    raise FileNotFoundError(f"Fresh rank file for '{method_alias}' not found under {ranks_dir}")

def _run_all(main_script: Path, config_glob: str, seeds: str):
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"
    cfgs = sorted(Path(".").glob(config_glob))
    if not cfgs:
        raise ValueError(f"No configs matched glob: {config_glob}")
    for cfg in cfgs:
        for seed in [int(x) for x in seeds.split(",")]:
            print(f"[run] {cfg.stem} seed={seed}")
            cmd = [sys.executable, str(main_script), "--config", str(cfg),
                   "--action", "train", "--model_filepath", f"outputs/models/{cfg.stem}_seed{seed}.pkl"]
            subprocess.run(cmd, check=True, env=env)

def main():
    ap = argparse.ArgumentParser(description="Retrain, evaluate (ID-join), and emit Top-K per method.")
    ap.add_argument("--main-script", required=True)
    ap.add_argument("--config-glob", required=True)
    ap.add_argument("--methods", required=True)
    ap.add_argument("--ranks-dir", required=True)
    ap.add_argument("--topk", type=int, default=50)
    ap.add_argument("--out-dir", default="outputs/end2end")
    ap.add_argument("--seeds", default="123")
    ap.add_argument("--aligned-priority", default=None, help="Aligned priority CSV (id + priority_score)")
    ap.add_argument("--aligned-id-col", default="hm_rsID")
    ap.add_argument("--ranks-id-col", default="gene_variant_id")
    args = ap.parse_args()

    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    out_topk = out_dir / "topk"; out_topk.mkdir(parents=True, exist_ok=True)
    out_metrics = out_dir / "metrics"; out_metrics.mkdir(parents=True, exist_ok=True)

    # Step 1: retrain & predict (fresh ranks)
    _run_all(Path(args.main_script), args.config_glob, args.seeds)

    # Step 2: prepare aligned labels if provided
    prio = None
    used_set = None
    if args.aligned_priority:
        prio = pd.read_csv(args.aligned_priority)
        if args.aligned_id_col not in prio.columns or "priority_score" not in prio.columns:
            raise ValueError("Aligned priority must have columns: "
                             f"{args.aligned_id_col}, priority_score")
        prio["_id_norm"] = prio[args.aligned_id_col].astype(str).str.strip().str.lower()
        prio = prio[["_id_norm","priority_score"]]
        used_set = set(prio["_id_norm"].tolist())
        print(f"[info] aligned priority loaded; ids={len(used_set)}")

    # Step 3: per-method processing
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    summary = []

    for alias in methods:
        rf = _find_rank_file(Path(args.ranks_dir), alias)
        df = pd.read_csv(rf)
        # ensure id + rank
        if args.ranks_id_col not in df.columns:
            raise ValueError(f"{rf} is missing '{args.ranks_id_col}'")
        if "rank" not in df.columns:
            if "score" not in df.columns:
                raise ValueError(f"{rf} must contain 'rank' or 'score'")
            # derive rank by descending score
            df["rank"] = pd.to_numeric(df["score"], errors="coerce")
            df["rank"] = df["rank"].rank(ascending=False, method="average").astype(int)

        # annotate used_for_fit if labels provided
        if used_set is not None:
            ids_norm = df[args.ranks_id_col].astype(str).str.strip().str.lower()
            df["used_for_fit"] = ids_norm.isin(used_set).values

        # write Top-K (overall)
        top_all = df.sort_values("rank", ascending=True).head(args.topk)
        top_all.to_csv(out_topk / f"top{args.topk}_{alias}.csv", index=False)
        print(f"[ok] wrote: {out_topk / f'top{args.topk}_{alias}.csv'}")

        # labeled/unlabeled views
        if used_set is not None:
            top_lab = df[df["used_for_fit"]].sort_values("rank", ascending=True).head(args.topk)
            top_unlab = df[~df["used_for_fit"]].sort_values("rank", ascending=True).head(args.topk)
            top_lab.to_csv(out_topk / f"top{args.topk}_{alias}_labeled.csv", index=False)
            top_unlab.to_csv(out_topk / f"top{args.topk}_{alias}_unlabeled.csv", index=False)
            print(f"[ok] wrote: {out_topk / f'top{args.topk}_{alias}_labeled.csv'}")
            print(f"[ok] wrote: {out_topk / f'top{args.topk}_{alias}_unlabeled.csv'}")

        # metrics by ID-join
        metrics = {"method": alias}
        if prio is not None and "score" in df.columns:
            tmp = df[[args.ranks_id_col, "score"]].copy()
            tmp["_id_norm"] = tmp[args.ranks_id_col].astype(str).str.strip().str.lower()
            joined = tmp.merge(prio, on="_id_norm", how="inner")
            n = len(joined)
            metrics["n_overlap"] = int(n)
            if n >= 2:
                y_true = pd.to_numeric(joined["priority_score"], errors="coerce").to_numpy()
                y_pred = pd.to_numeric(joined["score"], errors="coerce").to_numpy()
                msk = ~np.isnan(y_true) & ~np.isnan(y_pred)
                y_true = y_true[msk]; y_pred = y_pred[msk]
                if len(y_true) >= 2:
                    metrics["r2"] = float(r2_score(y_true, y_pred))
                    metrics["mse"] = float(mean_squared_error(y_true, y_pred))
                    metrics["mae"] = float(mean_absolute_error(y_true, y_pred))
                else:
                    metrics["note"] = "not enough overlap after NaN drop"
            else:
                metrics["note"] = "not enough overlap to compute metrics"

        # write per-method metrics json
        with open(out_metrics / f"metrics_{alias}.json", "w") as f:
            json.dump(metrics, f, indent=2)

        summary.append(metrics)

    # write summary csv
    pd.DataFrame(summary).to_csv(out_metrics / "metrics_summary.csv", index=False)
    print(f"[ok] wrote metrics summary: {out_metrics / 'metrics_summary.csv'}")

if __name__ == "__main__":
    main()
