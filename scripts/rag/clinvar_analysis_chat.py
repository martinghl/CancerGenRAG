#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AnalysisChat (single CSV, Gemini) - ClinVar evidence + JSON narrative
---------------------------------------------------------------------
Default: SINGLE CSV input that already contains a 'gene' column (e.g., top_variants_*.csv or *_with_gene.csv).

What it does
- Reads ONE CSV and extracts (rsid, gene) pairs:
  * rsid parsed from an ID column (default: gene_variant_id); if an explicit rsid column exists, it will use it.
  * gene taken from a gene column (default: gene).
- Queries NCBI ClinVar via Entrez (snp -> clinvar links) for each pair (rate-limited; email required).
- Flattens key fields to a tidy CSV.
- Builds a compact prompt and asks Gemini for a structured JSON narrative (highlights, summary, caveats).
- If --dry-run or missing SDK/API key, saves PROMPT_ONLY instead of calling the model.

Install
  pip install biopython pandas google-genai
  export GEMINI_API_KEY="YOUR_KEY"   # or pass --gemini-api-key

Examples
  python AnalysisChat_single_with_gene.py \
    --csv pgs/top_variants_20251111_141049.csv \
    --email you@domain \
    --disease "lung cancer" \
    --outdir analysis \
    --model gemini-2.5-flash

  # Dry-run (no API calls)
  python AnalysisChat_single_with_gene.py \
    --csv geo_aggregated_ranks_with_gene.csv \
    --email you@domain \
    --disease "lung cancer" \
    --outdir analysis \
    --dry-run
"""

from pathlib import Path
from textwrap import dedent
import argparse
import datetime as _dt
import json
import os
import re
import time
from typing import Optional, Tuple, Dict, Any

import pandas as pd


def _now_tag() -> str:
    return _dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def _infer_rsid(s: str) -> Optional[str]:
    s = "" if s is None else str(s)
    m = re.search(r"(rs\d+)", s, flags=re.IGNORECASE)
    return m.group(1) if m else None


def _normalize_gene(g: str) -> Optional[str]:
    if g is None:
        return None
    s = str(g).strip()
    s = re.sub(r"[^A-Za-z0-9_-]", "", s)
    return s or None


def _read_pairs_single_csv(path: str, id_col: str, gene_col: str, rsid_col: Optional[str], max_rows: Optional[int]) -> pd.DataFrame:
    df = pd.read_csv(path)
    if gene_col not in df.columns:
        raise ValueError(f"Input CSV missing gene column '{gene_col}'. Available: {list(df.columns)}")
    if rsid_col and rsid_col in df.columns:
        rs = df[rsid_col].astype(str).apply(_infer_rsid)
    else:
        if id_col not in df.columns:
            # Try to auto-detect a column with 'gene_variant_id' or 'id' in the name
            alt = None
            for c in df.columns:
                if c.lower() in ("gene_variant_id", "variant_id", "id"):
                    alt = c; break
            if alt is None:
                raise ValueError(f"Input CSV missing id_col '{id_col}'. Available: {list(df.columns)}")
            id_col = alt
        rs = df[id_col].astype(str).apply(_infer_rsid)
    genes = df[gene_col].apply(_normalize_gene)
    out = pd.DataFrame({"rsid": rs, "gene": genes}).dropna().drop_duplicates()
    if max_rows:
        out = out.head(max_rows)
    return out.reset_index(drop=True)


# ------------------------- ClinVar via Entrez -------------------------

def _entrez_setup(email: str, ncbi_api_key: Optional[str] = None):
    try:
        from Bio import Entrez
    except Exception as e:
        raise RuntimeError("biopython is required. Install with: pip install biopython") from e
    Entrez.email = email
    key = ncbi_api_key or os.getenv("NCBI_API_KEY")
    if key:
        Entrez.api_key = key
    return Entrez


def _rate_delay(ncbi_api_key: Optional[str] = None) -> float:
    # NCBI policy: ~3 req/s without key (~0.34s), ~10 req/s with key (~0.11s)
    return 0.11 if (ncbi_api_key or os.getenv("NCBI_API_KEY")) else 0.34


def _snp_to_clinvar_ids(Entrez, rsid: str) -> list[str]:
    ids = []
    try:
        h = Entrez.esearch(db="snp", term=rsid, retmode="xml")
        rec = Entrez.read(h)
        snp_ids = rec.get("IdList", [])
        time.sleep(0.01)
        if not snp_ids:
            return []
        h2 = Entrez.elink(dbfrom="snp", db="clinvar", id=",".join(snp_ids), retmode="xml")
        lnk = Entrez.read(h2)
        for ls in lnk:
            for db in ls.get("LinkSetDb", []):
                if db.get("DbTo") == "clinvar":
                    ids.extend([l["Id"] for l in db.get("Link", [])])
        return list(dict.fromkeys(ids))
    except Exception:
        return []


def _pyify(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {k: _pyify(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_pyify(v) for v in obj]
    try:
        json.dumps(obj); return obj
    except Exception:
        return str(obj)


def _flatten_clinvar_docsum(doc: Dict[str, Any]) -> Dict[str, Any]:
    d = _pyify(doc)
    out = {
        "clinvar_id": d.get("uid") or d.get("Id") or d.get("id"),
        "title": d.get("title") or d.get("Title"),
        "review_status": d.get("review_status") or d.get("ReviewStatus") or d.get("review_status_summary"),
        "last_update": d.get("last_update_date") or d.get("LastUpdated"),
        "accession": d.get("accession") or d.get("rcvaccession") or d.get("accession_version"),
        "clinical_significance": None,
        "conditions": None,
        "raw_json": json.dumps(d, ensure_ascii=False),
    }
    cs = d.get("clinical_significance") or d.get("ClinicalSignificance")
    if isinstance(cs, dict):
        out["clinical_significance"] = cs.get("description") or cs.get("Description") or str(cs)
    elif isinstance(cs, str):
        out["clinical_significance"] = cs

    conds = []
    traitset = d.get("traitset") or d.get("TraitSet") or d.get("traits") or d.get("TraitList")
    if isinstance(traitset, list):
        for t in traitset:
            if isinstance(t, dict):
                name = t.get("name") or t.get("trait_name") or t.get("preferred_name")
                if name:
                    conds.append(str(name))
    elif isinstance(traitset, dict):
        name = traitset.get("name") or traitset.get("trait_name") or traitset.get("preferred_name")
        if name:
            conds.append(str(name))
    out["conditions"] = "; ".join(sorted(set(conds))) if conds else None
    return out


def query_clinvar_pairs(pairs_df: pd.DataFrame, email: str, ncbi_api_key: Optional[str], verbose: bool = True) -> pd.DataFrame:
    Entrez = _entrez_setup(email, ncbi_api_key)
    delay = _rate_delay(ncbi_api_key)
    rows = []

    for idx, row in pairs_df.iterrows():
        rsid = row["rsid"]; gene = row["gene"]
        if verbose:
            print(f"[{idx+1}/{len(pairs_df)}] Query rsid={rsid} gene={gene}")
        time.sleep(delay)

        clin_ids = _snp_to_clinvar_ids(Entrez, rsid)
        if not clin_ids:
            rows.append({"rsid": rsid, "gene": gene, "clinvar_id": None, "title": None, "review_status": None,
                         "last_update": None, "accession": None, "clinical_significance": None, "conditions": None, "raw_json": None})
            continue

        try:
            h = Entrez.esummary(db="clinvar", id=",".join(clin_ids), retmode="xml")
            summ = Entrez.read(h)
            docs = summ.get("DocumentSummarySet", {}).get("DocumentSummary", [])
            if isinstance(docs, dict):
                docs = [docs]
        except Exception:
            docs = []

        if not docs:
            rows.append({"rsid": rsid, "gene": gene, "clinvar_id": None, "title": None, "review_status": None,
                         "last_update": None, "accession": None, "clinical_significance": None, "conditions": None, "raw_json": None})
        else:
            for d in docs:
                flat = _flatten_clinvar_docsum(d)
                flat["rsid"] = rsid; flat["gene"] = gene
                rows.append(flat)

    return pd.DataFrame(rows)


# ------------------------- Gemini narrative -------------------------

def _build_prompt_for_gemini(disease: str, clinvar_csv: Path, max_items: int = 25) -> str:
    df = pd.read_csv(clinvar_csv)
    df2 = df.dropna(subset=["clinvar_id", "clinical_significance"]).copy()
    keep = [c for c in ["rsid", "gene", "clinical_significance", "review_status", "conditions", "last_update"] if c in df2.columns]
    df2 = df2[keep].head(max_items)
    lines = [",".join(keep)]
    for _, r in df2.iterrows():
        lines.append(",".join([str(r.get(c, "")) for c in keep]))
    table_text = "\n".join(lines)

    sys_msg = dedent(f"""
You are a conservative genetics analyst. Summarize ClinVar evidence for top {disease} variants.
Avoid speculation; synthesize what the table implies about clinical significance and disease traits.
Return valid JSON only (no markdown).
""").strip()

    task = dedent("""
Required JSON keys:
- highlights: array of {rsid, gene, significance, condition, review} objects for 5-10 notable items
- summary: short paragraph synthesizing patterns in significance and review levels
- caveats: 2-4 bullets on limitations (coverage gaps, conflicting interpretations, etc.)
""").strip()

    data = f"ClinVar preview (CSV):\n{table_text}"
    return f"{sys_msg}\n\n{task}\n\n{data}"


def _run_gemini(prompt: str, model: str, out_dir: Path, tag: str, api_key: Optional[str]) -> Optional[Path]:
    try:
        from google import genai
        from google.genai import types
    except Exception as e:
        txt = out_dir / f"analysis_out_ANALY_{tag}_PROMPT_ONLY.txt"
        txt.write_text(prompt, encoding="utf-8")
        print(f"[WARN] google-genai SDK not available: {e}. Saved PROMPT_ONLY: {txt}")
        return None

    key = api_key or os.getenv("GEMINI_API_KEY") or os.getenv("GOOGLE_API_KEY")
    if not key:
        txt = out_dir / f"analysis_out_ANALY_{tag}_PROMPT_ONLY.txt"
        txt.write_text(prompt, encoding="utf-8")
        print(f"[WARN] No API key. Saved PROMPT_ONLY: {txt}")
        return None

    client = genai.Client(api_key=key)

    try:
        schema = types.Schema(
            type=types.Type.OBJECT,
            properties={
                "highlights": types.Schema(
                    type=types.Type.ARRAY,
                    items=types.Schema(
                        type=types.Type.OBJECT,
                        properties={
                            "rsid": types.Schema(type=types.Type.STRING),
                            "gene": types.Schema(type=types.Type.STRING),
                            "significance": types.Schema(type=types.Type.STRING),
                            "condition": types.Schema(type=types.Type.STRING),
                            "review": types.Schema(type=types.Type.STRING),
                        },
                        required=["rsid", "gene"]
                    )
                ),
                "summary": types.Schema(type=types.Type.STRING),
                "caveats": types.Schema(type=types.Type.STRING),
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
        txt = out_dir / f"analysis_out_ANALY_{tag}_PROMPT_ONLY.txt"
        txt.write_text(prompt, encoding="utf-8")
        print(f"[WARN] Gemini call failed: {e}. Saved PROMPT_ONLY: {txt}")
        return None

    # Validate JSON
    try:
        obj = json.loads(text)
    except Exception:
        import re as _re
        m = _re.search(r"\{.*\}\s*$", text, flags=_re.DOTALL)
        if not m:
            raw = out_dir / f"analysis_out_ANALY_{tag}_RAW.txt"
            raw.write_text(text, encoding="utf-8")
            print(f"[OK] Saved RAW text (non-JSON): {raw}")
            return raw
        obj = json.loads(m.group(0))

    out_json = out_dir / f"analysis_out_ANALY_{tag}.json"
    out_json.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[OK] Saved JSON narrative: {out_json}")
    return out_json


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="AnalysisChat (single CSV, Gemini)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="ClinVar evidence retrieval for a SINGLE CSV with (rsid,gene) derivable columns",
    )
    parser.add_argument("--csv", required=True, help="Single CSV path (must contain a gene column).")
    parser.add_argument("--id-col", default="gene_variant_id", help="Column to parse rsid from if --rsid-col not given (default gene_variant_id).")
    parser.add_argument("--gene-col", default="gene", help="Gene column (default gene).")
    parser.add_argument("--rsid-col", default=None, help="Optional explicit rsid column (e.g., 'rsid').")
    parser.add_argument("--max-rows", type=int, default=60, help="Max pairs to query ClinVar for (after dedupe).")
    parser.add_argument("--email", required=True, help="Entrez requires your email for API usage.")
    parser.add_argument("--ncbi-api-key", default=None, help="Optional NCBI API key (or set env NCBI_API_KEY).")
    parser.add_argument("--disease", default="lung cancer", help="Disease label for the narrative.")
    parser.add_argument("--outdir", default="outputs/analysis", help="Output directory.")
    parser.add_argument("--model", default="gemini-2.5-flash", help="Gemini model for narrative.")
    parser.add_argument("--gemini-api-key", default=None, help="Gemini API key; or set GEMINI_API_KEY/GOOGLE_API_KEY env.")
    parser.add_argument("--dry-run", action="store_true", help="Build pairs and prompt only; skip NCBI & Gemini calls.")
    args = parser.parse_args(argv)

    out_dir = Path(args.outdir); out_dir.mkdir(parents=True, exist_ok=True)
    tag = _now_tag()

    pairs = _read_pairs_single_csv(args.csv, args.id_col, args.gene_col, args.rsid_col, args.max_rows)
    if pairs.empty:
        raise RuntimeError("No (rsid,gene) pairs after parsing. Check column names or file content.")
    pairs_path = out_dir / f"pairs_{tag}.csv"
    pairs.to_csv(pairs_path, index=False)
    print(f"[OK] Parsed pairs: {pairs_path}  (n={len(pairs)})")

    if args.dry_run:
        preview_csv = out_dir / f"clinvar_results_{tag}_EMPTY.csv"
        pairs.assign(clinvar_id=None, title=None, review_status=None, last_update=None,
                     accession=None, clinical_significance=None, conditions=None).to_csv(preview_csv, index=False)
        prompt = _build_prompt_for_gemini(args.disease, preview_csv, max_items=min(25, len(pairs)))
        txt = out_dir / f"analysis_out_ANALY_{tag}_PROMPT_ONLY.txt"
        txt.write_text(prompt, encoding="utf-8")
        print(f"[DRY-RUN] Saved PROMPT_ONLY: {txt}")
        return

    res_df = query_clinvar_pairs(pairs, email=args.email, ncbi_api_key=args.ncbi_api_key, verbose=True)
    res_csv = out_dir / f"clinvar_results_{tag}.csv"
    res_df.to_csv(res_csv, index=False)
    print(f"[OK] ClinVar results: {res_csv}")

    nf = res_df[res_df["clinvar_id"].isna()][["rsid", "gene"]].drop_duplicates()
    if not nf.empty:
        nf_csv = out_dir / f"clinvar_not_found_{tag}.csv"
        nf.to_csv(nf_csv, index=False)
        print(f"[INFO] Not found pairs: {nf_csv}")

    prompt = _build_prompt_for_gemini(args.disease, res_csv, max_items=25)
    _ = _run_gemini(prompt, args.model, out_dir, tag, api_key=args.gemini_api_key)


if __name__ == "__main__":
    main()
