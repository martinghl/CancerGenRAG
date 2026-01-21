# CancerGen-RAG

CancerGen-RAG is a modular workflow for **lung cancer variant / gene prioritization and interpretation**.

Core idea (see `docs/CancerGen-RAG_paper.pdf`): **decouple prioritization from interpretation**
1) **Prioritization:** consolidate heterogeneous genomic evidence into a transparent consensus ranking (rank aggregation).
2) **Interpretation:** after ranking, run evidence synthesis (LLM + curated sources) for the top candidates.

This repository is a cleaned, GitHub-ready re-organization of the project code.

---

## What is included

- `scripts/preprocess/` – optional dbSNP SQLite builder + rsID imputation utilities
- `scripts/aggregation/` – ID alignment, label coverage checks, postprocessing helpers, evaluation utilities, plots
- `scripts/rag/` – LLM-driven evidence synthesis utilities (OpenAI / Gemini), Enrichr enrichment, STRING PPI, ClinVar lookups
- `r/` – R utilities for PGS data preparation/annotation (Bioconductor/Quincunx-based)
- `docs/` – paper PDF + pipeline notes
- `data/sample/` – tiny templates / sample files for sanity checks

## What is intentionally NOT included (you will add later)

- **GWAS priority score table** (e.g., `LC_GWAS_Priority_Scores.csv`).
  - A template is provided at `data/sample/LC_GWAS_Priority_Scores_TEMPLATE.csv`.
- Large/private raw data (full PGS scoring files, big intermediate PKLs, dbSNP VCF, etc.).
- (Optional) The *rank-aggregation engines* themselves, if you run them elsewhere.
  - This repo assumes you can provide method-specific `aggregated_ranks.csv` files and offers tooling around them.

---

## Repository layout

```
CancerGen-RAG/
  scripts/
    preprocess/   # dbSNP SQLite + rsID imputation (optional)
    aggregation/  # alignment, evaluation, plots
    rag/          # LLM + evidence sources (ClinVar / Enrichr / STRING)
  r/              # R-based PGS parsing/annotation utilities
  data/sample/    # templates + tiny samples
  docs/           # paper + notes
  notebooks/      # exploratory notebooks (kept for provenance)
  legacy_original/# original tree (for auditability)
```

---

## Quickstart (Python)

### 1) Create env & install dependencies

```bash
python -m venv .venv
source .venv/bin/activate   # (Windows: .venv\\Scripts\\activate)

pip install -r requirements/core.txt
# Optional modules:
pip install -r requirements/preprocess.txt
pip install -r requirements/rag.txt
```

### 2) Optional: build a dbSNP SQLite DB (for rsID imputation)

```bash
python scripts/preprocess/build_dbsnp_sqlite.py \
  --vcf path/to/All_20180418.filtered.vcf.gz \
  --db  data/variants.db \
  --batch-size 100000
```

### 3) Optional: impute rsIDs in annotated CSVs

```bash
python scripts/preprocess/impute_rsid_from_sqlite.py \
  --db data/variants.db \
  --input-dir path/to/Annotated \
  --output-dir path/to/Annotated_with_rsid
```

---

## Prioritization utilities

### A) Build crosswalk from an intermediate `data_dict.pkl`

```bash
python scripts/aggregation/build_crosswalk_from_pkl.py \
  --data-pkl path/to/data_dict.pkl \
  --out-csv outputs/crosswalk.csv
```

### B) Align your priority-score table to canonical IDs

```bash
python scripts/aggregation/align_priority_ids.py \
  --priority-csv data/LC_GWAS_Priority_Scores.csv \
  --priority-id-col variant_id \
  --priority-score-col priority_score \
  --crosswalk-csv outputs/crosswalk.csv \
  --target-key rsID \
  --out-csv outputs/LC_GWAS_Priority_Scores_aligned.csv
```

### C) Check label coverage / overlap

```bash
python scripts/aggregation/check_labels_coverage.py \
  --priority-csv outputs/LC_GWAS_Priority_Scores_aligned.csv \
  --priority-id-col rsID \
  --data-pkl path/to/data_dict.pkl \
  --out-dir outputs/label_coverage
```

### D) Postprocess an aggregated ranking (Top-50)

If you already have a method-specific `aggregated_ranks.csv` (see template at `data/sample/aggregated_ranks_TEMPLATE.csv`):

```bash
python scripts/aggregation/unsup_postprocess_top50.py \
  --output_path outputs/unsup \
  --prefix GEO \
  --topn 50
```

### E) Plot comparison vs GWAS priority

```bash
python scripts/aggregation/plot_unsup_comparison.py \
  --unsup-dir outputs/unsup \
  --methods GEO,MC-3,THURSTONE,RRA,BARD,CEMC_SPEARMAN,CEMC_KENDALL \
  --priority-csv data/LC_GWAS_Priority_Scores.csv \
  --priority-id-col variant_id \
  --priority-col priority_score \
  --topn 12 \
  --ext pdf
```

---

## RAG / LLM utilities

### Environment variables

Copy `config/example.env` to `.env` and fill in the keys you need.

### OpenAI example

```bash
python scripts/rag/pgs_chat.py \
  --input outputs/unsup/GEO/aggregated_ranks.csv \
  --disease "lung cancer" \
  --topn 25 \
  --model gpt-4.1-mini \
  --output-dir outputs/pgs
```

### Gemini + ClinVar example

ClinVar script requires an email for NCBI Entrez.

```bash
python scripts/rag/clinvar_analysis_chat.py \
  --csv outputs/unsup/GEO/aggregated_ranks.csv \
  --disease "lung cancer" \
  --outdir outputs/clinvar \
  --email you@example.com \
  --model gemini-2.5-flash
```

---

