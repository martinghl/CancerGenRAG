#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GeneChat (single CSV, Gemini, v3) - Enrichment + STRING PPI (TSV + high-res PNG/SVG) + JSON
-------------------------------------------------------------------------------
- Single CSV input with a 'gene' column.
- Validates Enrichr libraries; robust against pandas/gseapy version quirks.
- STRING PPI:
    * Saves TSV via /api/tsv/network
    * Saves a ready-made network figure via STRING API:
        - PNG: /api/image/network
        - High-res PNG: /api/highres_image/network
        - Vector SVG: /api/svg/network

Install:
  pip install --upgrade gseapy pandas requests matplotlib google-genai
  export GEMINI_API_KEY="YOUR_KEY"

Usage (network figure auto when --string is set):
  python GeneChat_single_with_gene_v3.py \
    --csv pgs/top_variants_*.csv \
    --disease "lung cancer" \
    --outdir genechat \
    --string \
    --model gemini-2.5-flash
"""

from pathlib import Path
from textwrap import dedent
import argparse
import datetime as _dt
import json
import os
import re
import difflib
from typing import List, Optional

import pandas as pd


def _now_tag() -> str:
    return _dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def _read_genes_from_csv(path: str, gene_col: str, rank_col: Optional[str]) -> List[str]:
    df = pd.read_csv(path)
    if gene_col not in df.columns:
        raise ValueError(f"Input CSV missing gene column '{gene_col}'. Available: {list(df.columns)}")
    if rank_col and rank_col in df.columns:
        order_cols = [rank_col, gene_col]
    elif "global_rank" in df.columns:
        order_cols = ["global_rank", gene_col]
    elif "aggregated_rank" in df.columns:
        order_cols = ["aggregated_rank", gene_col]
    else:
        order_cols = [gene_col]
    sdf = df.sort_values(order_cols, na_position="last")
    genes, seen = [], set()
    for g in sdf[gene_col].astype(str):
        g2 = re.sub(r"[^A-Za-z0-9_-]", "", g.strip())
        if len(g2) < 2:
            continue
        k = g2.upper()
        if k in seen:
            continue
        seen.add(k)
        genes.append(g2)
    return genes


def _resolve_libraries(requested: List[str]) -> List[str]:
    try:
        import gseapy as gp
    except Exception as e:
        raise RuntimeError("gseapy is required. Install with: pip install --upgrade gseapy") from e
    try:
        avail = gp.get_library_name()
    except Exception as e:
        print(f"[WARN] Could not fetch Enrichr library list: {e}")
        return requested
    avail_set = set(avail)
    valid = [lib for lib in requested if lib in avail_set]
    missing = [lib for lib in requested if lib not in avail_set]
    if missing:
        print(f"[WARN] These libraries are not recognized by Enrichr: {missing}")
        for lib in missing:
            sug = difflib.get_close_matches(lib, avail, n=3, cutoff=0.5)
            if sug:
                print(f"       Suggestions for '{lib}': {sug}")
    if not valid:
        picks = []
        for key in ["KEGG", "Reactome", "GO_Biological_Process", "GO_Molecular_Function", "GO_Cellular_Component"]:
            cand = [n for n in avail if key in n]
            if cand:
                cand.sort(reverse=True)
                picks.append(cand[0])
        valid = [p for p in picks if p]
        print(f"[INFO] Falling back to libraries: {valid}")
    return valid


def run_enrichment(genes: List[str], libraries: List[str], out_dir: Path, tag: str, top_k: int = 12) -> Path:
    try:
        import gseapy as gp
    except Exception as e:
        raise RuntimeError("gseapy is required. Install with: pip install --upgrade gseapy") from e

    libs = _resolve_libraries(libraries)

    all_rows = []
    for lib in libs:
        try:
            enr = gp.enrichr(gene_list=genes, gene_sets=lib, outdir=None, cutoff=1.0, no_plot=True, verbose=False)
            res = enr.results if hasattr(enr, "results") else None
            if res is None or res.empty:
                print(f"[WARN] No results for library: {lib}")
                continue
            res = res.copy(); res["gene_set"] = lib
            all_rows.append(res)
            try:
                _make_plots_for_library(res, lib, out_dir, tag, top_k=top_k)
            except Exception as pe:
                print(f"[WARN] Plotting failed for {lib}: {pe}")
        except AttributeError as ae:
            msg = str(ae)
            if "append" in msg:
                raise RuntimeError(
                    "Your gseapy/pandas combo is incompatible (pandas>=2 removed DataFrame.append). "
                    "Fix: `pip install --upgrade gseapy` (recommended) or `pip install 'pandas<2.0'`."
                )
            else:
                print(f"[WARN] Enrichr failed for {lib}: {ae}")
        except Exception as ie:
            print(f"[WARN] Enrichr failed for {lib}: {ie}")
            continue

    if not all_rows:
        raise RuntimeError("No enrichment results across libraries. Try upgrading gseapy and/or changing libraries via --libraries.")

    comb = pd.concat(all_rows, ignore_index=True)

    if "Adjusted P-value" in comb.columns:
        comb["adj_p"] = comb["Adjusted P-value"]
    else:
        cand = [c for c in comb.columns if "adjust" in c.lower() and "p" in c.lower()]
        comb["adj_p"] = comb[cand[0]] if cand else None
    if "Combined Score" in comb.columns:
        comb["combined_score"] = comb["Combined Score"]

    out_csv = out_dir / f"pathway_enrichment_{tag}.csv"
    comb.to_csv(out_csv, index=False)
    return out_csv


def _parse_overlap_to_ratio(overlap: str) -> float:
    try:
        if isinstance(overlap, str) and "/" in overlap:
            a, b = overlap.split("/", 1)
            a = float(a.strip()); b = float(b.strip())
            return a / b if b > 0 else 0.0
    except Exception:
        pass
    return 0.0


def _make_plots_for_library(res_df: "pd.DataFrame", lib: str, out_dir: Path, tag: str, top_k: int = 12, ext: str = "png"):
    import matplotlib.pyplot as plt
    import math

    df = res_df.copy()
    if "adj_p" not in df.columns:
        if "Adjusted P-value" in df.columns:
            df["adj_p"] = df["Adjusted P-value"]
        else:
            cols = [c for c in df.columns if "adjust" in c.lower() and "p" in c.lower()]
            if cols:
                df["adj_p"] = df[cols[0]]
            else:
                raise ValueError("Adjusted p-value column not found in enrichment results.")

    df = df.sort_values("adj_p", na_position="last").head(top_k)
    terms = df["Term"].astype(str) if "Term" in df.columns else df.iloc[:, 0].astype(str)

    # Bar: -log10(adj p)
    x = -df["adj_p"].astype(float).clip(lower=1e-300).apply(lambda v: math.log10(v))
    fig, ax = plt.subplots(figsize=(8, max(3, 0.4 * len(df))))
    ax.barh(range(len(df)), x)
    ax.set_yticks(range(len(df))); ax.set_yticklabels(list(terms))
    ax.set_xlabel("-log10(adj p)"); ax.set_title(f"Top terms - {lib}")
    plt.tight_layout()
    bar_path = out_dir / f"enrich_bar_{lib}_{tag}.{ext}"
    fig.savefig(bar_path, dpi=300); plt.close(fig)

    # Dot: Combined Score vs -log10(adj p), size by Overlap ratio
    cs = df["Combined Score"].astype(float) if "Combined Score" in df.columns else df.get("combined_score", pd.Series([0.0] * len(df))).astype(float)
    y2 = -df["adj_p"].astype(float).clip(lower=1e-300).apply(lambda v: math.log10(v))
    sizes = (df["Overlap"].astype(str).apply(_parse_overlap_to_ratio) if "Overlap" in df.columns else pd.Series([0.5] * len(df)))
    sizes = (sizes.fillna(0.0) + 0.05) * 800.0
    fig2, ax2 = plt.subplots(figsize=(7, 5))
    ax2.scatter(cs.fillna(0.0), y2, s=sizes)
    for i, term in enumerate(list(terms)):
        ax2.text(cs.iloc[i] if pd.notna(cs.iloc[i]) else 0.0, y2.iloc[i], term, fontsize=8, ha="left", va="bottom")
    ax2.set_xlabel("Combined Score"); ax2.set_ylabel("-log10(adj p)"); ax2.set_title(f"Dot plot - {lib}")
    plt.tight_layout()
    dot_path = out_dir / f"enrich_dot_{lib}_{tag}.{ext}"
    fig2.savefig(dot_path, dpi=300); plt.close(fig2)


def _build_prompt(disease: str, enrich_csv_path: Path, top_terms: int = 12) -> str:
    df = pd.read_csv(enrich_csv_path)
    df2 = df.copy()
    if "adj_p" not in df2.columns:
        if "Adjusted P-value" in df2.columns:
            df2["adj_p"] = df2["Adjusted P-value"]
        else:
            cand = [c for c in df2.columns if "adjust" in c.lower() and "p" in c.lower()]
            if cand:
                df2["adj_p"] = df2[cand[0]]
            else:
                df2["adj_p"] = None

    keep_cols = [c for c in ["gene_set", "Term", "adj_p", "Combined Score", "Genes"] if c in df2.columns]
    prev = df2.sort_values("adj_p", na_position="last").head(top_terms)[keep_cols]
    lines = [",".join(keep_cols)]
    for _, row in prev.iterrows():
        lines.append(",".join([str(row.get(c, "")) for c in keep_cols]))
    table_text = "\n".join(lines)

    sys_msg = dedent(f"""
You are a precise genetics analyst. Summarize recurring biological themes in the following enrichment results for {disease}.
Be concise and conservative; avoid speculative mechanisms. Return valid JSON only (no markdown).
""").strip()

    task = dedent("""
Required JSON keys:
- themes: array of {term, library, adj_p, summary} items for top terms with very short 1-sentence summaries
- summary: one short paragraph weaving the themes together
- verification: 2-4 bullet points suggesting how to verify (replication, LD expansion, cross-dataset consistency, etc.)
""").strip()

    data = f"Top enriched preview (CSV):\n{table_text}"
    return f"{sys_msg}\n\n{task}\n\n{data}"


def _run_gemini(prompt: str, model: str, out_dir: Path, tag: str, api_key: Optional[str]):
    try:
        from google import genai
        from google.genai import types
    except Exception as e:
        txt_path = out_dir / f"analysis_out_GENEENRICH_{tag}_PROMPT_ONLY.txt"
        txt_path.write_text(prompt, encoding="utf-8")
        print(f"[WARN] google-genai SDK not available: {e}. Saved PROMPT_ONLY: {txt_path}")
        return None

    key = api_key or os.getenv("GEMINI_API_KEY") or os.getenv("GOOGLE_API_KEY")
    if not key:
        txt_path = out_dir / f"analysis_out_GENEENRICH_{tag}_PROMPT_ONLY.txt"
        txt_path.write_text(prompt, encoding="utf-8")
        print(f"[WARN] No API key. Saved PROMPT_ONLY: {txt_path}")
        return None

    client = genai.Client(api_key=key)

    try:
        schema = types.Schema(
            type=types.Type.OBJECT,
            properties={
                "themes": types.Schema(
                    type=types.Type.ARRAY,
                    items=types.Schema(
                        type=types.Type.OBJECT,
                        properties={
                            "term": types.Schema(type=types.Type.STRING),
                            "library": types.Schema(type=types.Type.STRING),
                            "adj_p": types.Schema(type=types.Type.NUMBER),
                            "summary": types.Schema(type=types.Type.STRING),
                        },
                        required=["term", "library"]
                    )
                ),
                "summary": types.Schema(type=types.Type.STRING),
                "verification": types.Schema(type=types.Type.STRING),
            },
            required=["summary"]
        )
    except Exception:
        schema = None

    try:
        if schema is not None:
            resp = client.models.generate_content(
                model=model,
                contents=prompt,
                config=types.GenerateContentConfig(
                    response_mime_type="application/json",
                    response_schema=schema,
                ),
            )
        else:
            resp = client.models.generate_content(
                model=model,
                contents=prompt,
                config={"response_mime_type": "application/json"},
            )
        text = resp.text or ""
    except Exception as e:
        txt_path = out_dir / f"analysis_out_GENEENRICH_{tag}_PROMPT_ONLY.txt"
        txt_path.write_text(prompt, encoding="utf-8")
        print(f"[WARN] Gemini call failed: {e}. Saved PROMPT_ONLY: {txt_path}")
        return None

    try:
        obj = json.loads(text)
    except Exception:
        import re as _re
        m = _re.search(r"\{.*\}\s*$", text, flags=_re.DOTALL)
        if not m:
            raw = out_dir / f"analysis_out_GENEENRICH_{tag}_RAW.txt"
            raw.write_text(text, encoding="utf-8")
            print(f"[OK] Saved RAW text (non-JSON): {raw}")
            return raw
        obj = json.loads(m.group(0))

    out_json = out_dir / f"analysis_out_GENEENRICH_{tag}.json"
    out_json.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[OK] Saved JSON: {out_json}")
    return out_json


def fetch_string_network_tsv(genes: List[str], out_dir: Path, tag: str, species: int = 9606, required_score: int = 400):
    try:
        import requests
    except Exception as e:
        print(f"[WARN] requests not installed: {e}")
        return None
    if not genes:
        print("[WARN] empty gene list for STRING"); return None
    base = "https://string-db.org/api/tsv/network"
    identifiers = "%0d".join(genes)
    params = {"identifiers": identifiers, "species": species, "required_score": required_score, "add_nodes": 0}
    try:
        r = requests.get(base, params=params, timeout=30); r.raise_for_status()
    except Exception as e:
        print(f"[WARN] STRING TSV request failed: {e}"); return None
    out_tsv = out_dir / f"string_network_{tag}.tsv"
    out_tsv.write_text(r.text, encoding="utf-8")
    return out_tsv


def fetch_string_network_figure(
    genes: List[str],
    out_dir: Path,
    tag: str,
    species: int = 9606,
    required_score: int = 400,
    add_white_nodes: int = 0,
    network_flavor: str = "confidence",
    hide_disconnected: int = 0,
    figure_format: str = "highres_png",
    hide_node_labels: int = 0,
    label_font_size: Optional[int] = None,
    figure_name: Optional[str] = None,
):
    """Fetch a ready-made STRING network figure.

    Supported figure_format values:
      - "png"         -> output_format=image (standard PNG)
      - "highres_png" -> output_format=highres_image (high-res PNG)
      - "svg"         -> output_format=svg (vector)
    """
    try:
        import requests
    except Exception as e:
        print(f"[WARN] requests not installed: {e}")
        return None
    if not genes:
        print("[WARN] empty gene list for STRING")
        return None

    if network_flavor not in ("confidence", "evidence", "actions", "action"):
        network_flavor = "confidence"

    ff = (figure_format or "highres_png").lower().strip()
    if ff == "svg":
        output_format, ext = "svg", "svg"
    elif ff in ("highres", "highres_png", "highres-image", "highres_image"):
        output_format, ext = "highres_image", "png"
    else:
        output_format, ext = "image", "png"

    # STRING API URL pattern: /api/[output-format]/network
    base = f"https://string-db.org/api/{output_format}/network"
    identifiers = "%0d".join(genes)
    params = {
        "identifiers": identifiers,
        "species": species,
        "required_score": required_score,
        "network_flavor": network_flavor,
        "hide_disconnected_nodes": hide_disconnected,
        # For multi-protein queries, STRING shows interactions among the inputs.
        # To extend the neighborhood, use add_white_nodes / add_color_nodes.
        "add_white_nodes": add_white_nodes,
        "hide_node_labels": hide_node_labels,
    }
    if label_font_size is not None:
        params["custom_label_font_size"] = int(label_font_size)

    try:
        r = requests.get(base, params=params, timeout=60)
        r.raise_for_status()
    except Exception as e:
        print(f"[WARN] STRING figure request failed: {e}")
        return None

    fname = figure_name if figure_name else f"string_network_{network_flavor}_{tag}.{ext}"
    out_fig = out_dir / fname
    out_fig.write_bytes(r.content)
    return out_fig


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="GeneChat (single CSV, Gemini, v3)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Run Enrichr on a single CSV, export STRING PPI (TSV + PNG), and produce a Gemini JSON narrative",
    )
    parser.add_argument("--csv", required=True, help="Input CSV (must contain a gene column).")
    parser.add_argument("--gene-col", default="gene", help="Name of the gene column in the CSV.")
    parser.add_argument("--rank-col", default=None, help="Optional rank column to preserve order (default tries global_rank then aggregated_rank).")
    parser.add_argument("--max-genes", type=int, default=200, help="Max number of genes to use (after dedupe, in order).")
    parser.add_argument("--libraries", default="GO_Biological_Process_2021,GO_Molecular_Function_2021,GO_Cellular_Component_2021,KEGG_2021_Human,Reactome_2016", help="Comma-separated Enrichr libraries.")
    parser.add_argument("--top-k", type=int, default=12, help="Top terms per library for plots and prompt.")
    parser.add_argument("--outdir", default="outputs/genechat", help="Output directory.")
    parser.add_argument("--plots-ext", default="png", choices=["png", "svg"], help="Plot file extension.")
    parser.add_argument("--string", action="store_true", help="Also fetch STRING PPI TSV and PNG.")
    parser.add_argument("--string-species", type=int, default=9606, help="STRING species id (9606=human).")
    parser.add_argument("--string-score", type=int, default=400, help="STRING required_score 0-1000.")
    parser.add_argument("--string-add-nodes", type=int, default=0, help="STRING add_nodes parameter (0..N).")
    parser.add_argument("--string-flavor", default="confidence", choices=["confidence", "evidence", "actions", "action"], help="STRING network_flavor for PNG.")
    parser.add_argument("--string-hide-disconnected", type=int, default=0, help="Hide disconnected nodes (1=yes, 0=no).")
    parser.add_argument(
        "--string-fig-format",
        default="svg",
        choices=["png", "highres_png", "svg"],
        help="STRING network figure output format (highres_png gives a larger, sharper PNG; svg is vector).",
    )
    parser.add_argument(
        "--string-hide-labels",
        type=int,
        default=0,
        help="Hide node labels in STRING network figure (1=yes, 0=no).",
    )
    parser.add_argument(
        "--string-label-font-size",
        type=int,
        default=None,
        help="Custom label font size (5-50) for STRING network figure; default uses STRING's built-in size.",
    )
    parser.add_argument("--disease", default="lung cancer", help="Disease label for narrative.")
    parser.add_argument("--model", default="gemini-2.5-flash", help="Gemini model id.")
    parser.add_argument("--api-key", default=None, help="Gemini API key; or set GEMINI_API_KEY/GOOGLE_API_KEY.")
    parser.add_argument("--dry-run", action="store_true", help="Save prompt only; do not call model.")
    args = parser.parse_args(argv)

    out_dir = Path(args.outdir); out_dir.mkdir(parents=True, exist_ok=True)
    tag = _now_tag()

    genes_all = _read_genes_from_csv(args.csv, args.gene_col, args.rank_col)
    if not genes_all:
        raise RuntimeError("No genes parsed from CSV.")
    genes = genes_all[: args.max_genes]
    print(f"[INFO] Using {len(genes)} genes (deduped, ordered).")

    req_libs = [s.strip() for s in args.libraries.split(",") if s.strip()]
    enrich_csv = run_enrichment(genes, req_libs, out_dir, tag, top_k=args.top_k)
    print(f"[OK] Enrichment table: {enrich_csv}")

    if args.string:
        tsv = fetch_string_network_tsv(genes, out_dir, tag, species=args.string_species, required_score=args.string_score)
        if tsv:
            print(f"[OK] STRING network TSV: {tsv}")
        fig = fetch_string_network_figure(
            genes,
            out_dir,
            tag,
            species=args.string_species,
            required_score=args.string_score,
            add_white_nodes=args.string_add_nodes,
            network_flavor=args.string_flavor,
            hide_disconnected=args.string_hide_disconnected,
            figure_format=args.string_fig_format,
            hide_node_labels=args.string_hide_labels,
            label_font_size=args.string_label_font_size,
        )
        if fig:
            print(f"[OK] STRING network figure: {fig}")

    prompt = _build_prompt(args.disease, enrich_csv, top_terms=args.top_k)
    if args.dry_run:
        txt = out_dir / f"analysis_out_GENEENRICH_{tag}_PROMPT_ONLY.txt"
        txt.write_text(prompt, encoding="utf-8")
        print(f"[DRY-RUN] Saved PROMPT_ONLY: {txt}")
        return

    _ = _run_gemini(prompt, args.model, out_dir, tag, api_key=args.api_key)


if __name__ == "__main__":
    main()
