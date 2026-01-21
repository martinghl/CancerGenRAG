#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_labels_coverage.py

Purpose
-------
Given (A) a priority score CSV and (B) a data .pkl (e.g., data_dict.pkl),
report how many variants are effectively "labeled" under typical rules and
how much ID overlap exists between the two sources (priority vs model universe).

What it does
------------
1) Load priority CSV:
   - Use --priority-id-col (default: variant_id) and --priority-score-col (default: priority_score)
   - Clean IDs and deduplicate.
   - Count distinct IDs, non-null scores, zeros vs non-zeros (if numeric).

2) Load data .pkl:
   - Expect a dict-like {PGSxxxx: DataFrame} or a single DataFrame.
   - Search for ID columns in precedence order (--pkl-id-cols, default: hm_rsID,rsID,variant_id).
   - Union IDs across all tables; clean and deduplicate to form the "universe" used by the model.

3) Overlap analysis:
   - Intersect priority IDs with the universe IDs.
   - If priority score is numeric, estimate "labeled (y!=0)" count on the intersection.

4) Write a compact report to stdout and also save CSV/JSON artifacts under --out-dir.

Usage
-----
python check_labels_coverage.py \
  --priority-csv data/LC_GWAS_Priority_Scores.csv \
  --priority-id-col variant_id \
  --priority-score-col priority_score \
  --data-pkl /path/to/data_dict.pkl \
  --out-dir outputs/check_labels

Optional flags:
  --pkl-id-cols "hm_rsID,rsID,variant_id"    (precedence order when scanning DataFrames)
  --zero-as-unlabeled                        (treat score==0 as unlabeled; default True)
  --sample-k 20                              (how many sample IDs to write out for each set)
"""
import argparse
import os
import json
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional, Union

import pandas as pd
import numpy as np

# -------------------- helpers --------------------
def _to_str_id(x) -> Optional[str]:
    """Normalize an ID to a clean lowercase string, or None if invalid."""
    if pd.isna(x):
        return None
    try:
        s = str(x).strip()
    except Exception:
        return None
    if s == "" or s.lower() in {"nan", "none", "null"}:
        return None
    return s  # keep case; caller may lower() if desired

def _clean_id_series(series: pd.Series, lower: bool = True) -> pd.Series:
    out = series.map(_to_str_id)
    out = out.dropna()
    if lower:
        out = out.str.lower()
    # remove pure "chr:pos" with spaces trimmed; keep as-is, but caller should be aware
    return out

def _load_priority(path: Path, id_col: str, score_col: str, lower_ids: bool = True) -> Tuple[pd.DataFrame, Set[str]]:
    df = pd.read_csv(path)
    if id_col not in df.columns:
        raise ValueError(f"priority id_col '{id_col}' not found in {path}")
    if score_col not in df.columns:
        raise ValueError(f"priority score_col '{score_col}' not found in {path}")
    ids = _clean_id_series(df[id_col], lower=lower_ids)
    df = df.assign(_norm_id=ids)
    df = df.dropna(subset=["_norm_id"]).drop_duplicates("_norm_id")
    id_set = set(df["_norm_id"].tolist())
    return df, id_set

def _scan_df_for_idcol(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None

def _load_universe_from_pkl(path: Path, candidates: List[str], lower_ids: bool = True) -> Tuple[Set[str], Dict[str, int]]:
    obj = pd.read_pickle(path)
    id_counts: Dict[str, int] = {}
    uni: Set[str] = set()

    def _add_ids(series: pd.Series):
        nonlocal id_counts, uni
        s = _clean_id_series(series, lower=lower_ids)
        for v in s:
            uni.add(v)
            id_counts[v] = id_counts.get(v, 0) + 1

    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(v, pd.DataFrame):
                id_col = _scan_df_for_idcol(v, candidates)
                if id_col is not None:
                    _add_ids(v[id_col])
    elif isinstance(obj, pd.DataFrame):
        id_col = _scan_df_for_idcol(obj, candidates)
        if id_col is not None:
            _add_ids(obj[id_col])
    else:
        raise ValueError(f"Unsupported pkl payload type: {type(obj)}")

    return uni, id_counts

def _numeric_series(s: pd.Series) -> Optional[pd.Series]:
    try:
        out = pd.to_numeric(s, errors="coerce")
        return out
    except Exception:
        return None

# -------------------- main --------------------
def parse_args():
    p = argparse.ArgumentParser(description="Check label coverage vs priority CSV & data.pkl (ID overlap, counts).")
    p.add_argument("--priority-csv", type=str, required=True)
    p.add_argument("--priority-id-col", type=str, default="variant_id")
    p.add_argument("--priority-score-col", type=str, default="priority_score")
    p.add_argument("--data-pkl", type=str, required=True)
    p.add_argument("--pkl-id-cols", type=str, default="hm_rsID,rsID,variant_id",
                   help="Comma-separated precedence of ID columns to look for inside DataFrames.")
    p.add_argument("--zero-as-unlabeled", action="store_true", default=True)
    p.add_argument("--no-zero-as-unlabeled", dest="zero_as_unlabeled", action="store_false")
    p.add_argument("--out-dir", type=str, default="outputs/check_labels")
    p.add_argument("--sample-k", type=int, default=20)
    return p.parse_args()

def main():
    args = parse_args()
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    # 1) Priority CSV
    prio_df, prio_ids = _load_priority(Path(args.priority_csv), args.priority_id_col, args.priority_score_col, lower_ids=True)
    prio_scores_num = _numeric_series(prio_df[args.priority_score_col])

    prio_stats = {
        "priority_rows": int(len(prio_df)),
        "priority_unique_ids": int(len(prio_ids)),
        "priority_score_nonnull": int(prio_scores_num.notna().sum()) if prio_scores_num is not None else None,
        "priority_score_zeros": int((prio_scores_num == 0).sum()) if prio_scores_num is not None else None,
        "priority_score_nonzeros": int((prio_scores_num != 0).sum()) if prio_scores_num is not None else None,
    }

    # 2) Universe from data.pkl
    candidates = [c.strip() for c in args.pkl_id_cols.split(",") if c.strip()]
    uni_ids, id_counts = _load_universe_from_pkl(Path(args.data_pkl), candidates, lower_ids=True)
    uni_stats = {
        "universe_ids": int(len(uni_ids)),
        "id_occurrence_top5": sorted([(k, v) for k, v in id_counts.items()], key=lambda x: -x[1])[:5],
    }

    # 3) Overlap
    inter = prio_ids.intersection(uni_ids)
    inter_stats = {
        "intersection": int(len(inter)),
    }

    # 4) Estimate labeled counts on the intersection (y != 0 if numeric)
    labeled_guess = None
    zeros_on_inter = None
    if prio_scores_num is not None:
        prio_on_inter = prio_df.set_index("_norm_id").loc[list(inter)][args.priority_score_col]
        prio_on_inter_num = pd.to_numeric(prio_on_inter, errors="coerce")
        zeros_on_inter = int((prio_on_inter_num == 0).sum())
        labeled_guess = int((prio_on_inter_num != 0).sum())

    # 5) Save samples
    inter_sample = list(sorted(inter))[:args.sample_k]
    prio_only = list(sorted(prio_ids - uni_ids))[:args.sample_k]
    uni_only = list(sorted(uni_ids - prio_ids))[:args.sample_k]

    pd.Series(inter_sample, name="id").to_csv(out_dir / "sample_intersection.csv", index=False)
    pd.Series(prio_only, name="id_priority_only").to_csv(out_dir / "sample_priority_only.csv", index=False)
    pd.Series(uni_only, name="id_universe_only").to_csv(out_dir / "sample_universe_only.csv", index=False)

    # 6) Write summary JSON + echo to stdout
    summary = {
        "priority": prio_stats,
        "universe": uni_stats,
        "overlap": inter_stats,
        "labeled_estimate_on_intersection": {
            "labeled_guess_y_ne_0": labeled_guess,
            "zeros_on_intersection": zeros_on_inter,
        },
        "params": {
            "priority_csv": str(Path(args.priority_csv).resolve()),
            "data_pkl": str(Path(args.data_pkl).resolve()),
            "pkl_id_cols": candidates,
            "zero_as_unlabeled": bool(args.zero_as_unlabeled),
        }
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("\n==== LABEL COVERAGE REPORT ====")
    print(json.dumps(summary, indent=2))
    print(f"\nWrote samples & summary to: {out_dir}\n")

if __name__ == "__main__":
    main()
