#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PGSChat (patched) — variant narrative from aggregated rankings
-------------------------------------------------------------
- Inputs: one or more aggregated rank CSVs with columns:
    * gene_variant_id
    * aggregated_rank
- Output:
    * top_variants_<timestamp>.csv  (global Top-N with inferred rsID/gene)
    * analysis_out_PGS_<timestamp>.txt (LLM narrative)  OR *_PROMPT_ONLY.txt when --dry-run/LLM unavailable
- Notes:
    * Disease is set via --disease (default: "lung cancer").
    * Uses OpenAI Python SDK (Responses API preferred, Chat Completions fallback).
    * Reads API key from OPENAI_API_KEY; safe to run with --dry-run to skip API calls.

Example:
    python PGSChat_fixed.py \      --agg-csv /path/to/agg1.csv /path/to/agg2.csv \      --disease "lung cancer" \      --topn 50 \      --output-dir /path/to/outputs/pgs \      --model gpt-4o-mini \      --dry-run
"""
from pathlib import Path
from textwrap import dedent
import argparse
import datetime as _dt
import os
import re
import pandas as pd
from typing import List, Tuple, Optional


def _now_tag() -> str:
    return _dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def _infer_rsid_and_gene(gene_variant_id: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Pull rsID and a plausible gene symbol from a free-form gene_variant_id string.
    Heuristics for patterns like GENE|rs123, GENE_rs123, rs123|GENE, etc.
    """
    s = str(gene_variant_id or "")
    m = re.search(r"(rs\d+)", s, flags=re.IGNORECASE)
    rsid = m.group(1) if m else None

    # Tokenize on common separators
    tokens = re.split(r"[|:_\-./\s]+", s)
    gene = None
    for t in tokens:
        if not t or re.match(r"^rs\d+$", t, flags=re.IGNORECASE):
            continue
        # Prefer tokens that look like HGNC-like symbols
        if re.match(r"^[A-Za-z0-9]{2,20}$", t):
            gene = t
            break

    return (rsid, gene)


def _read_and_merge(csv_paths: List[str]) -> pd.DataFrame:
    """
    Read multiple CSVs containing columns ['gene_variant_id','aggregated_rank'].
    Outer-merge on 'gene_variant_id'. Produce per-file rank columns and a 'global_rank' = min across ranks.
    """
    merged = None
    rank_cols = []
    for p in csv_paths:
        pth = Path(p)
        if not pth.exists():
            raise FileNotFoundError(f"CSV not found: {pth}")
        df = pd.read_csv(pth)
        if "gene_variant_id" not in df.columns or "aggregated_rank" not in df.columns:
            raise ValueError(f"{pth.name} must contain columns: gene_variant_id, aggregated_rank")
        col_name = f"aggregated_rank__{pth.stem}"
        df = df[["gene_variant_id", "aggregated_rank"]].rename(columns={"aggregated_rank": col_name})
        rank_cols.append(col_name)
        merged = df if merged is None else pd.merge(merged, df, on="gene_variant_id", how="outer")

    # Compute global rank = min across available rank cols
    merged[rank_cols] = merged[rank_cols].apply(pd.to_numeric, errors="coerce")
    merged["global_rank"] = merged[rank_cols].min(axis=1, skipna=True)

    # Parse rsID & gene
    parsed = merged["gene_variant_id"].apply(_infer_rsid_and_gene)
    merged["rsid"] = parsed.apply(lambda t: t[0])
    merged["gene"] = parsed.apply(lambda t: t[1])

    # Order nicely
    ordered_cols = ["global_rank", "gene_variant_id", "rsid", "gene"] + rank_cols
    merged = merged[ordered_cols].sort_values(by=["global_rank", "gene_variant_id"], na_position="last").reset_index(drop=True)
    return merged


def _build_prompt(disease: str, top_df: pd.DataFrame, topn: int) -> str:
    # Keep a compact table preview for the prompt (limit rows to topn and a few columns)
    preview_cols = ["global_rank", "gene_variant_id", "rsid", "gene"]
    preview = top_df.loc[: topn - 1, preview_cols].copy()
    table_lines = ["global_rank,gene_variant_id,rsid,gene"]
    for _, row in preview.iterrows():
        table_lines.append(",".join([str(row.get(c, "")) for c in preview_cols]))
    table_text = "\n".join(table_lines)

    sys_msg = dedent(
        f'''
        You are a precise scientific writer. Write a concise, factual narrative (≤400 words) about the top-ranked variants associated with {disease}.
        Requirements:
        - Summarize key biological themes across the genes (pathways, mechanisms) WITHOUT over-claiming causality.
        - Do NOT invent data. Only infer modestly from known functions of listed genes.
        - Prefer disciplined, manuscript-ready tone; no marketing language.
        - If rsIDs are missing for some rows, acknowledge coverage limits.
        '''
    ).strip()

    user_msg = dedent(
        f'''
        Top-ranked variants table (first {topn} rows):

        {table_text}

        Task:
        1) Provide a paragraph that contextualizes the list for {disease}.
        2) Highlight recurring functional themes (e.g., cell-cycle, DNA repair, immune regulation) if applicable.
        3) End with a 3–5 bullet list of caveats/next steps (e.g., LD expansion, fine-mapping, replication, enrichment checks).
        '''
    ).strip()

    prompt = f"<system>\n{sys_msg}\n</system>\n\n<user>\n{user_msg}\n</user>"
    return prompt


def _call_llm(prompt: str, model: str, temperature: float, max_tokens: int):
    """
    Try OpenAI Responses API first; fallback to Chat Completions. Raise on failure.
    Returns (text, which_api).
    """
    try:
        from openai import OpenAI  # type: ignore
    except Exception as e:
        raise RuntimeError(
            "OpenAI SDK not available. Install with `pip install openai` or run with --dry-run."
        ) from e

    client = OpenAI(api_key=os.getenv("OPENAI_API_KEY", None))

    # Try Responses API
    try:
        resp = client.responses.create(
            model=model,
            input=[
                {"role": "system", "content": "You are a helpful assistant."},
                {"role": "user", "content": prompt},
            ],
            temperature=temperature,
            max_output_tokens=max_tokens,
        )
        out_text = getattr(resp, "output_text", None)
        if not out_text:
            try:
                out_text = resp.output[0].content[0].text  # type: ignore[attr-defined]
            except Exception:
                out_text = None
        if out_text:
            return out_text, "responses"
    except Exception:
        pass

    # Fallback: Chat Completions
    resp = client.chat.completions.create(
        model=model,
        messages=[
            {"role": "system", "content": "You are a helpful assistant."},
            {"role": "user", "content": prompt},
        ],
        temperature=temperature,
        max_tokens=max_tokens,
    )
    out_text = resp.choices[0].message.content  # type: ignore[attr-defined]
    return out_text, "chat.completions"


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="PGSChat (patched)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="PGSChat (patched) — variant narrative from aggregated rankings",
    )
    parser.add_argument(
        "--agg-csv",
        nargs="+",
        required=True,
        help="Path(s) to aggregated rank CSV files (must contain columns: gene_variant_id, aggregated_rank).",
    )
    parser.add_argument("--disease", default="lung cancer", help="Disease label for narrative context.")
    parser.add_argument("--topn", type=int, default=25, help="Number of top rows to include in the analysis and prompt.")
    parser.add_argument("--output-dir", default="outputs/pgs", help="Directory for outputs.")
    parser.add_argument("--model", default="gpt-4o-mini", help="OpenAI model name.")
    parser.add_argument("--temperature", type=float, default=0.2, help="Sampling temperature.")
    parser.add_argument("--max-tokens", type=int, default=800, help="Max tokens for generation (approx.).")
    parser.add_argument("--dry-run", action="store_true", help="Do not call any LLM; save prompt only.")

    args = parser.parse_args(argv)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = _now_tag()

    # Read + merge
    merged = _read_and_merge(args.agg_csv)

    # Save a complete merged table as provenance
    merged_path = out_dir / f"merged_aggregated_{tag}.csv"
    merged.to_csv(merged_path, index=False)

    # Top-N slice for narrative
    top_df = merged.sort_values("global_rank", na_position="last").head(args.topn).reset_index(drop=True)
    top_path = out_dir / f"top_variants_{tag}.csv"
    top_df.to_csv(top_path, index=False)

    # Build prompt
    prompt = _build_prompt(args.disease, merged.sort_values("global_rank", na_position="last"), args.topn)

    if args.dry_run:
        out_txt = out_dir / f"analysis_out_PGS_{tag}_PROMPT_ONLY.txt"
        out_txt.write_text(prompt, encoding="utf-8")
        print(f"[DRY-RUN] Wrote prompt only:\n  {out_txt}\nTop-N table:\n  {top_path}\nMerged table:\n  {merged_path}")
        return

    # Call LLM
    try:
        text, api_used = _call_llm(prompt, model=args.model, temperature=args.temperature, max_tokens=args.max_tokens)
        out_txt = out_dir / f"analysis_out_PGS_{tag}.txt"
        out_txt.write_text(text, encoding="utf-8")
        print(f"[OK:{api_used}] Saved:\n  {out_txt}\nTop-N table:\n  {top_path}\nMerged table:\n  {merged_path}")
    except Exception as e:
        # Save prompt for manual use
        out_txt = out_dir / f"analysis_out_PGS_{tag}_PROMPT_ONLY.txt"
        out_txt.write_text(prompt, encoding="utf-8")
        print(f"[WARN] LLM call failed: {e}\nSaved prompt instead:\n  {out_txt}\nTop-N table:\n  {top_path}\nMerged table:\n  {merged_path}")


if __name__ == "__main__":
    main()
