#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
experiments_unsup_new.py — Orchestrates unsupervised runs and enforces
Top-50 merging with preselected top-10 via unsup_postprocess_top50.py.

What changed (Nov 2025):
- Accepts a global --preselected_file and forwards it to each run.
- Parses machine-readable "RESULTS:{...}" lines from child to build a summary CSV.
- Keeps your previous CLI shape and sensible defaults.

Examples:
  python experiments_unsup_new.py     --data_path /path/to/data_dict.pkl     --methods GEO,MC-3,THURSTONE,RRA,BARD,CEMC     --seeds 0     --output_path ./outputs/unsup     --preselected_file ./outputs/preselected_top10.csv
"""
from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional


THIS_DIR = Path(__file__).resolve().parent

import pandas as pd


LOG_DIR = Path("logs/unsup")
OUT_DIR = Path("./outputs/unsup")


def _now_ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run unsupervised methods × seeds and summarize.")
    p.add_argument("--data_path", type=str, required=True)
    p.add_argument("--methods", type=str, default="GEO,MC-3,THURSTONE,RRA,BARD,CEMC")
    p.add_argument("--seeds", type=str, default="0")
    p.add_argument("--output_path", type=str, default=str(OUT_DIR))
    p.add_argument("--model_dir", type=str, default=None)

    # new: preselected & postprocess knobs (forwarded)
    p.add_argument("--preselected_file", type=str, default=None)
    p.add_argument("--id_col", type=str, default="gene_variant_id")
    p.add_argument("--rank_col", type=str, default="aggregated_rank")
    p.add_argument("--topn", type=int, default=50)
    p.add_argument("--write_debug", type=int, default=1)
    return p.parse_args()


def _build_cmd(
    data_path: str,
    method: str,
    seed: int,
    output_path: str,
    model_dir: Optional[str],
    preselected_file: Optional[str],
    id_col: str,
    rank_col: str,
    topn: int,
    write_debug: int,
) -> List[str]:
    # model filepath (keep compatibility)
    model_dir = model_dir or os.path.join(output_path, "models")
    _ensure_dir(Path(model_dir))
    model_path = os.path.join(model_dir, f"{method}_seed{seed}.pkl")

    cmd = [
        "python",
        str(THIS_DIR / "unsup_postprocess_top50.py"),
        "--data_path", data_path,
        "--config", method,
        "--prefix", method,
        "--output_path", output_path,
        "--model_filepath", model_path,
        "--id_col", id_col,
        "--rank_col", rank_col,
        "--topn", str(topn),
        "--write_debug", str(write_debug),
    ]
    if preselected_file:
        cmd += ["--preselected_file", preselected_file]
    return cmd


def _run_one(cmd: List[str]) -> Dict:
    """Run a single child process, parse 'RESULTS:{json}' line from stdout/stderr."""
    print("[RUN] " + " ".join(shlex.quote(x) for x in cmd), flush=True)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    out = (proc.stdout or "") + "\n" + (proc.stderr or "")
    # Try to find RESULTS JSON
    payload = {}
    for line in out.splitlines():
        if line.startswith("RESULTS:"):
            try:
                payload = json.loads(line[len("RESULTS:"):].strip())
            except Exception:
                pass
    # Persist raw logs next to LOG_DIR for inspection
    _ensure_dir(LOG_DIR)
    ts = _now_ts()
    with open(LOG_DIR / f"child_{ts}.log", "w", encoding="utf-8") as f:
        f.write(out)
    if proc.returncode != 0:
        payload["status"] = "ERROR"
        payload["returncode"] = proc.returncode
    else:
        payload["status"] = payload.get("status", "OK")
    return payload


def main():
    args = parse_args()
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    seeds = [int(s.strip()) for s in args.seeds.split(",") if s.strip()]

    _ensure_dir(OUT_DIR)
    _ensure_dir(LOG_DIR)

    rows: List[Dict] = []
    for m in methods:
        for s in seeds:
            cmd = _build_cmd(
                data_path=args.data_path,
                method=m,
                seed=s,
                output_path=args.output_path,
                model_dir=args.model_dir,
                preselected_file=args.preselected_file,
                id_col=args.id_col,
                rank_col=args.rank_col,
                topn=args.topn,
                write_debug=args.write_debug,
            )
            res = _run_one(cmd)
            res.update({"method": m, "seed": s})
            rows.append(res)

    df = pd.DataFrame(rows)
    # order columns
    front = [c for c in ["method", "seed", "status", "returncode", "n_final", "n_preselected", "topn", "output_csv"] if c in df.columns]
    others = [c for c in df.columns if c not in front]
    df = df[front + others]

    ts = _now_ts()
    _ensure_dir(OUT_DIR)
    summary_csv = Path(args.output_path) / f"unsup_summary_{ts}.csv"
    df.to_csv(summary_csv, index=False)
    print(f"[OK] wrote summary: {summary_csv}", flush=True)


if __name__ == "__main__":
    main()
