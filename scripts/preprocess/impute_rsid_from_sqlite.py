#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
impute_rsid_from_sqlite.py

Fill missing rsIDs in CSV files using a local SQLite lookup created by
`build_dbsnp_sqlite.py`.

Expected columns (default):
- hm_chr (chrom)
- hm_pos (position, 1-based)
- hm_rsID (rsID)

The script will:
- read all CSVs in --input-dir
- for rows with missing hm_rsID, query (hm_chr, hm_pos) in the SQLite DB
- write updated CSVs to --output-dir (same file names by default)

Example
-------
python scripts/preprocess/impute_rsid_from_sqlite.py \
  --db data/variants.db \
  --input-dir path/to/Annotated \
  --output-dir path/to/Annotated_with_rsid
"""
from __future__ import annotations

import argparse
import glob
import os
import sqlite3
from pathlib import Path
from typing import Dict, Tuple, Optional

import pandas as pd
from tqdm import tqdm


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _norm_chr(x: str) -> str:
    s = str(x).strip()
    # keep as-is (db may store '1' or 'chr1'); user can pre-normalize if needed
    return s


def _lookup(conn: sqlite3.Connection, chrom: str, pos: int) -> Optional[str]:
    c = conn.cursor()
    c.execute("SELECT rsid FROM variants WHERE chrom=? AND pos=? LIMIT 1", (chrom, int(pos)))
    row = c.fetchone()
    return row[0] if row else None


def impute_file(
    csv_path: Path,
    conn: sqlite3.Connection,
    chr_col: str,
    pos_col: str,
    rsid_col: str,
) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    if chr_col not in df.columns or pos_col not in df.columns:
        raise ValueError(f"{csv_path.name}: missing required columns: {chr_col}, {pos_col}")
    if rsid_col not in df.columns:
        df[rsid_col] = pd.NA

    # Identify rows needing imputation
    need = df[rsid_col].isna() | (df[rsid_col].astype(str).str.strip() == "") | (df[rsid_col].astype(str).str.lower() == "nan")
    if need.sum() == 0:
        return df

    # Cache lookups per (chrom,pos)
    cache: Dict[Tuple[str, int], Optional[str]] = {}
    idxs = df.index[need].tolist()
    for i in tqdm(idxs, desc=f"impute:{csv_path.name}", unit="row"):
        chrom = _norm_chr(df.at[i, chr_col])
        pos = int(df.at[i, pos_col])
        key = (chrom, pos)
        if key not in cache:
            cache[key] = _lookup(conn, chrom, pos)
        if cache[key]:
            df.at[i, rsid_col] = cache[key]
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True, help="SQLite DB path created by build_dbsnp_sqlite.py")
    ap.add_argument("--input-dir", required=True, help="Directory containing input CSVs")
    ap.add_argument("--output-dir", required=True, help="Directory to write updated CSVs")
    ap.add_argument("--pattern", default="*.csv", help="Glob pattern under input-dir (default: *.csv)")
    ap.add_argument("--chr-col", default="hm_chr", help="Chromosome column name")
    ap.add_argument("--pos-col", default="hm_pos", help="Position column name (1-based)")
    ap.add_argument("--rsid-col", default="hm_rsID", help="rsID column name")
    args = ap.parse_args()

    db_path = Path(args.db)
    if not db_path.exists():
        raise FileNotFoundError(f"DB not found: {db_path}")

    in_dir = Path(args.input_dir)
    out_dir = Path(args.output_dir)
    _ensure_dir(out_dir)

    files = sorted([Path(p) for p in glob.glob(str(in_dir / args.pattern))])
    if not files:
        raise FileNotFoundError(f"No files matched: {in_dir / args.pattern}")

    conn = sqlite3.connect(str(db_path))
    try:
        for fp in files:
            df2 = impute_file(fp, conn, args.chr_col, args.pos_col, args.rsid_col)
            out_fp = out_dir / fp.name
            df2.to_csv(out_fp, index=False)
            print(f"[ok] wrote: {out_fp}")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
