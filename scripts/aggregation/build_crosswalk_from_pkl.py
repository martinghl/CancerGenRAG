#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_crosswalk_from_pkl.py

Create a variant ID crosswalk from a data .pkl that may contain one or more
DataFrames (e.g., {PGSxxxx: df}). The goal is to map between columns like:
  - hm_rsID, rsID, variant_id, chr, pos
and output a unified table where each row corresponds to one record from the pkl,
carrying as many identifier columns as available.

Usage
-----
python build_crosswalk_from_pkl.py \
  --data-pkl /path/to/data_dict.pkl \
  --out-csv outputs/crosswalk.csv \
  --id-cols "hm_rsID,rsID,variant_id" \
  --extra-cols "chr,pos"

Notes
-----
- We keep rows as-is (no grouping); this allows many-to-one relationships.
- You can later group/deduplicate on your preferred key (e.g., rsID) yourself.
"""
import argparse
from pathlib import Path
from typing import List, Dict, Any, Optional

import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description="Build a variant ID crosswalk CSV from a .pkl payload.")
    p.add_argument("--data-pkl", type=str, required=True)
    p.add_argument("--out-csv", type=str, required=True)
    p.add_argument("--id-cols", type=str, default="hm_rsID,rsID,variant_id")
    p.add_argument("--extra-cols", type=str, default="chr,pos")
    return p.parse_args()

def main():
    args = parse_args()
    obj = pd.read_pickle(Path(args.data_pkl))

    id_cols = [c.strip() for c in args.id_cols.split(",") if c.strip()]
    extra_cols = [c.strip() for c in args.extra_cols.split(",") if c.strip()]

    rows = []
    def _normalize(s: pd.Series) -> pd.Series:
        return s.astype(str).str.strip()

    if isinstance(obj, dict):
        for key, df in obj.items():
            if not isinstance(df, pd.DataFrame):
                continue
            cols = [c for c in id_cols + extra_cols if c in df.columns]
            if not cols:
                continue
            tmp = df[cols].copy()
            tmp = tmp.apply(_normalize)
            tmp.insert(0, "_source", str(key))
            rows.append(tmp)
    elif isinstance(obj, pd.DataFrame):
        cols = [c for c in id_cols + extra_cols if c in obj.columns]
        if cols:
            tmp = obj[cols].copy()
            tmp = tmp.apply(_normalize)
            tmp.insert(0, "_source", "root")
            rows.append(tmp)
    else:
        raise ValueError(f"Unsupported pkl type: {type(obj)}")

    if not rows:
        raise RuntimeError("No matching columns found in the .pkl for the provided id/extra columns.")

    out = pd.concat(rows, ignore_index=True)
    out.to_csv(Path(args.out_csv), index=False)
    print(f"[ok] wrote crosswalk to {args.out_csv}  (rows={len(out)})")

if __name__ == "__main__":
    main()
