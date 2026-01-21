#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
main_unsup.py — unsupervised pipeline runner (postprocess-focused)

What changed (Nov 2025):
- Adds *deterministic* post-processing to produce `aggregated_ranks.csv` by
  injecting a provided `preselected_top10` list BEFORE filling the remainder
  from the aggregated ranking, then re-ranking to Top-50.
- Writes a sidecar `aggregated_ranks.debug.csv` for auditability.
- Backups any pre-existing aggregated ranking to `aggregated_ranks.base.csv`.

This script keeps your CLI shape:
  --data_path PATH
  --config METHOD_NAME             (e.g., GEO, MC-3, THURSTONE, RRA, BARD, CEMC)
  --prefix METHOD_NAME             (used as folder name under --output_path)
  --output_path DIR                (default: ./outputs/unsup)
  --model_filepath PATH            (optional artifact path; we still create dirs)

New (optional) args for the Top-50 rule:
  --preselected_file PATH          (csv/txt; contains gene_variant_id or variant_id)
  --id_col COL                     (default: gene_variant_id)
  --rank_col COL                   (default: aggregated_rank)
  --topn N                         (default: 50)
  --write_debug {0,1}              (default: 1)

Expected inputs for postprocess:
- An "aggregated ranks" table produced earlier in the method run.
  We auto-discover it with this priority:
    <output_path>/<prefix>/aggregated_ranks.csv
    <output_path>/<prefix>_aggregated_ranks.csv
    <output_path>/aggregated_ranks.csv
- If none found, we fail with a clear message.

Outputs:
- <output_path>/<prefix>/aggregated_ranks.csv              (final for plotting)
- <output_path>/<prefix>/aggregated_ranks.debug.csv        (sidecar)
- <output_path>/<prefix>/aggregated_ranks.base.csv         (backup of original)

The final CSV has exactly two columns: [id_col, aggregated_rank],
where aggregated_rank is the *re-assigned* 1..K rank after combining
preselected + remainder and clipping to Top-N.
"""
from __future__ import annotations

import argparse
import json
import logging
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import pandas as pd

# --- Matplotlib safe (some envs import it downstream) -----------------------
try:
    import matplotlib  # type: ignore
    matplotlib.use("Agg")
except Exception:
    pass


# ----------------------------------------------------------------------------
# Utilities
# ----------------------------------------------------------------------------
def _now_ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _load_preselected(path: Optional[Path], id_col: str) -> List[str]:
    """Load up to 10 IDs from csv/txt. Keep first occurrence order & dedupe."""
    if path is None:
        return []
    if not path.exists():
        logging.warning("[postprocess] preselected_file not found: %s", path)
        return []
    ids: List[str] = []
    try:
        if path.suffix.lower() in {".csv", ".tsv"}:
            df = pd.read_csv(path)
            candidate_cols = [c for c in [id_col, "variant_id", "gene_variant_id", "id"] if c in df.columns]
            if not candidate_cols:
                raise ValueError(
                    f"Couldn't find an ID column in {list(df.columns)}; "
                    f"looked for {[id_col, 'variant_id', 'gene_variant_id', 'id']}"
                )
            col = candidate_cols[0]
            ids = df[col].astype(str).tolist()
        else:
            # plaintext: one id per line
            with open(path, "r", encoding="utf-8") as f:
                ids = [ln.strip() for ln in f if ln.strip()]
    except Exception as e:
        logging.error("[postprocess] failed reading preselected_file %s: %s", path, e)
        return []

    # normalize/cleanup & keep first-appearance order
    norm = []
    seen = set()
    for x in ids:
        y = str(x).strip()
        if not y:
            continue
        if y not in seen:
            norm.append(y)
            seen.add(y)
    # limit to top-10
    return norm[:10]


def _discover_agg_input(out_dir: Path, prefix: str) -> Optional[Path]:
    """Find an existing aggregated ranks file to postprocess."""
    candidates = [
        out_dir / prefix / "aggregated_ranks.csv",
        out_dir / f"{prefix}_aggregated_ranks.csv",
        out_dir / "aggregated_ranks.csv",
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


def _read_aggregated(table_path: Path, id_col: str, rank_col: str) -> pd.DataFrame:
    df = pd.read_csv(table_path)
    # tolerate different column names
    rename = {}
    if id_col not in df.columns:
        for c in ["gene_variant_id", "variant_id", "id"]:
            if c in df.columns:
                rename[c] = id_col
                break
    if rank_col not in df.columns:
        # Accept alternatives: 'rank', 'score' (then we sort by score)
        if "rank" in df.columns:
            rename["rank"] = rank_col
        elif "score" in df.columns:
            # sort ascending by score and create rank
            df = df.sort_values(by="score", ascending=True, kind="mergesort")
            df[rank_col] = range(1, len(df) + 1)
        else:
            # if nothing present, generate by current order
            df[rank_col] = range(1, len(df) + 1)

    if rename:
        df = df.rename(columns=rename)

    # normalize id dtype
    df[id_col] = df[id_col].astype(str)
    # ensure proper ranking sort
    df = df.sort_values(by=[rank_col, id_col], ascending=[True, True], kind="mergesort").reset_index(drop=True)
    return df[[id_col, rank_col]].copy()


def merge_preselected_topn(
    agg_df: pd.DataFrame, pre_ids: Iterable[str], id_col: str = "gene_variant_id", rank_col: str = "aggregated_rank", topn: int = 50
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Apply the rule:
      final = (preselected in-agg order by orig rank; 'missing' placed after found ones) + (remaining by agg rank),
      then clip to Top-N and reassign ranks 1..K.

    Returns (final_df_two_cols, debug_df).
    """
    # base copy
    agg_df = agg_df[[id_col, rank_col]].copy()
    agg_df[rank_col] = pd.to_numeric(agg_df[rank_col], errors="coerce")
    agg_df = agg_df.sort_values(by=[rank_col, id_col], ascending=[True, True], kind="mergesort").reset_index(drop=True)

    # map: id -> orig_rank
    orig_map = dict(zip(agg_df[id_col].tolist(), agg_df[rank_col].tolist()))

    # split preselected into found/missing, order by orig_rank, then id as tiebreaker
    pre_ids = list(dict.fromkeys([str(x).strip() for x in pre_ids if str(x).strip()]))  # dedupe keep order
    found = [(pid, orig_map[pid]) for pid in pre_ids if pid in orig_map]
    missing = [(pid, float("inf")) for pid in pre_ids if pid not in orig_map]
    found.sort(key=lambda x: (x[1], x[0]))
    # missing keep original order among themselves (stable)
    # Compose pre-block
    pre_block = [pid for pid, _r in found] + [pid for pid, _r in missing]

    # remaining agg IDs excluding pre_block
    pre_set = set(pre_block)
    rem_df = agg_df[~agg_df[id_col].isin(pre_set)].copy()
    rem_df = rem_df.sort_values(by=[rank_col, id_col], ascending=[True, True], kind="mergesort")

    # number to take
    k_pre = min(len(pre_block), topn)
    k_rem = max(0, topn - k_pre)

    # final order
    final_ids = pre_block[:k_pre] + rem_df[id_col].head(k_rem).tolist()

    # build outputs
    final_df = pd.DataFrame({id_col: final_ids})
    final_df["aggregated_rank"] = range(1, len(final_df) + 1)

    # debug
    def _orig_rank(x: str) -> float:
        return orig_map.get(x, float("inf"))

    debug_df = final_df.merge(agg_df.rename(columns={rank_col: "orig_rank"}), on=id_col, how="left")
    debug_df["orig_rank"] = debug_df["orig_rank"].fillna(float("inf"))
    debug_df["is_preselected"] = debug_df[id_col].isin(pre_set).astype(int)

    return final_df[[id_col, "aggregated_rank"]], debug_df[[id_col, "aggregated_rank", "orig_rank", "is_preselected"]]


# ----------------------------------------------------------------------------
# CLI + main
# ----------------------------------------------------------------------------
def parse_arguments() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Unsupervised main (postprocess Top-50 w/ preselected injection)")
    p.add_argument("--data_path", type=str, required=True)
    p.add_argument("--config", type=str, required=True, help="method name, e.g. GEO, MC-3, THURSTONE, RRA, BARD, CEMC")
    p.add_argument("--prefix", type=str, required=True, help="used as subfolder under --output_path")
    p.add_argument("--output_path", type=str, default="./outputs/unsup")
    p.add_argument("--model_filepath", type=str, default=None)

    # new
    p.add_argument("--preselected_file", type=str, default=None)
    p.add_argument("--id_col", type=str, default="gene_variant_id")
    p.add_argument("--rank_col", type=str, default="aggregated_rank")
    p.add_argument("--topn", type=int, default=50)
    p.add_argument("--write_debug", type=int, default=1)
    return p.parse_args()


def main(
    data_path: str,
    config: str,
    output_path: str,
    prefix: str,
    model_filepath: Optional[str] = None,
    preselected_file: Optional[str] = None,
    id_col: str = "gene_variant_id",
    rank_col: str = "aggregated_rank",
    topn: int = 50,
    write_debug: int = 1,
) -> dict:
    # Prepare dirs
    out_dir = Path(output_path)
    method_dir = out_dir / prefix
    model_filepath = model_filepath or (method_dir / "models" / f"{prefix}.pkl").as_posix()

    _ensure_dir(method_dir)
    _ensure_dir(Path(model_filepath).parent)

    # Discover an aggregated ranks table from upstream
    agg_in = _discover_agg_input(out_dir, prefix)
    if not agg_in:
        msg = (
            f"No aggregated ranks found to postprocess. Looked for:\n"
            f"  1) {out_dir/prefix/'aggregated_ranks.csv'}\n"
            f"  2) {out_dir/(prefix + '_aggregated_ranks.csv')}\n"
            f"  3) {out_dir/'aggregated_ranks.csv'}\n"
            f"Please ensure your method writes one of these before calling this script."
        )
        logging.error(msg)
        raise FileNotFoundError(msg)

    # Backup original
    base_copy = method_dir / "aggregated_ranks.base.csv"
    if not base_copy.exists():
        try:
            _ensure_dir(base_copy.parent)
            pd.read_csv(agg_in).to_csv(base_copy, index=False)
        except Exception as e:
            logging.warning("Failed to write base backup %s: %s", base_copy, e)

    # Load original + preselected
    agg_df = _read_aggregated(agg_in, id_col=id_col, rank_col=rank_col)
    pre_ids = _load_preselected(Path(preselected_file) if preselected_file else None, id_col=id_col)

    final_df, debug_df = merge_preselected_topn(
        agg_df=agg_df, pre_ids=pre_ids, id_col=id_col, rank_col=rank_col, topn=topn
    )

    # Write outputs
    final_path = method_dir / "aggregated_ranks.csv"
    final_df.to_csv(final_path, index=False)

    if write_debug:
        debug_path = method_dir / "aggregated_ranks.debug.csv"
        debug_df.to_csv(debug_path, index=False)
    else:
        debug_path = None

    # 'Model' side-effect: we still touch the model path so upstream scripts don't break
    Path(model_filepath).write_text(f"placeholder for {config}\n", encoding="utf-8")

    results = {
        "method": config,
        "prefix": prefix,
        "n_final": int(final_df.shape[0]),
        "n_preselected": int(len(pre_ids)),
        "topn": int(topn),
        "output_csv": str(final_path),
        "debug_csv": (str(debug_path) if debug_path else None),
        "backup_base": str(base_copy) if base_copy.exists() else None,
        "data_path": data_path,
    }

    # Emit a machine-readable line for the experiments script
    print("RESULTS:" + json.dumps(results, ensure_ascii=False), flush=True)
    logging.info("RESULTS:%s", json.dumps(results, ensure_ascii=False))
    return results


if __name__ == "__main__":
    # timestamped logfile under logs/unsup
    ts = _now_ts()
    logdir = Path("logs/unsup")
    _ensure_dir(logdir)
    logging.basicConfig(
        level=logging.INFO,
        filename=str(logdir / f"logfile_{ts}.log"),
        filemode="a+",
        format="%(asctime)-15s %(levelname)-8s %(message)s",
    )
    args = parse_arguments()
    main(
        data_path=args.data_path,
        config=args.config,
        output_path=args.output_path,
        prefix=args.prefix,
        model_filepath=args.model_filepath,
        preselected_file=args.preselected_file,
        id_col=args.id_col,
        rank_col=args.rank_col,
        topn=args.topn,
        write_debug=args.write_debug,
    )
