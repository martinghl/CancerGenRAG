#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_dbsnp_sqlite.py

Build a lightweight SQLite index for rsID lookup from a dbSNP VCF (bgzip + tabix).

Why this exists
---------------
Some PGS/annotation tables only contain (chr, pos) but not rsID. This script
creates a local (chrom, pos) -> rsID lookup DB to support fast rsID imputation.

Output schema
-------------
Table: variants
  - chrom TEXT
  - pos   INTEGER
  - rsid  TEXT
Primary index: (chrom, pos)

Table: processed_chromosomes
  - chrom TEXT PRIMARY KEY
  - last_pos INTEGER

The `processed_chromosomes` table enables resume/restart: if a chromosome was
partially processed, the next run will continue after `last_pos`.

Example
-------
python scripts/preprocess/build_dbsnp_sqlite.py \
  --vcf path/to/All_20180418.filtered.vcf.gz \
  --db  data/variants.db \
  --batch-size 100000
"""
from __future__ import annotations

import argparse
import os
import sqlite3
from pathlib import Path
from typing import Iterable, List, Tuple, Optional

import pysam
from tqdm import tqdm


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _init_db(conn: sqlite3.Connection) -> None:
    c = conn.cursor()
    c.execute(
        """
        CREATE TABLE IF NOT EXISTS variants (
            chrom TEXT NOT NULL,
            pos   INTEGER NOT NULL,
            rsid  TEXT NOT NULL,
            PRIMARY KEY (chrom, pos)
        )
        """
    )
    c.execute(
        """
        CREATE TABLE IF NOT EXISTS processed_chromosomes (
            chrom TEXT PRIMARY KEY,
            last_pos INTEGER
        )
        """
    )
    c.execute("CREATE INDEX IF NOT EXISTS idx_variants_chrom_pos ON variants(chrom, pos)")
    conn.commit()


def _get_resume_pos(conn: sqlite3.Connection, chrom: str) -> Optional[int]:
    c = conn.cursor()
    c.execute("SELECT last_pos FROM processed_chromosomes WHERE chrom = ?", (chrom,))
    row = c.fetchone()
    return int(row[0]) if row and row[0] is not None else None


def _set_resume_pos(conn: sqlite3.Connection, chrom: str, last_pos: int) -> None:
    c = conn.cursor()
    c.execute(
        """
        INSERT INTO processed_chromosomes (chrom, last_pos)
        VALUES (?, ?)
        ON CONFLICT(chrom) DO UPDATE SET last_pos = excluded.last_pos
        """,
        (chrom, int(last_pos)),
    )
    conn.commit()


def _iter_records(
    vcf: pysam.VariantFile, chrom: str, start_pos: int | None = None
) -> Iterable[Tuple[str, int, str]]:
    """
    Yield (chrom, pos, rsid) tuples. Skips records without rsID.
    """
    # pysam fetch uses 0-based start; positions in VCF are 1-based.
    if start_pos is None:
        it = vcf.fetch(chrom)
    else:
        it = vcf.fetch(chrom, start_pos)  # start_pos is 0-based here
    for rec in it:
        # VCF: rec.pos is 1-based int
        rid = rec.id
        if not rid:
            continue
        # Keep only rsIDs; dbSNP VCF may have multiple IDs separated by ';'
        rs = None
        for tok in str(rid).split(";"):
            tok = tok.strip()
            if tok.startswith("rs"):
                rs = tok
                break
        if not rs:
            continue
        yield chrom, int(rec.pos), rs


def _batched(iterable: Iterable[Tuple[str, int, str]], batch_size: int) -> Iterable[List[Tuple[str, int, str]]]:
    batch: List[Tuple[str, int, str]] = []
    for item in iterable:
        batch.append(item)
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def build_db(vcf_path: Path, db_path: Path, batch_size: int, chroms: List[str] | None, force: bool) -> None:
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF not found: {vcf_path}")
    # Tabix index expected as .tbi (bgzip)
    if not (vcf_path.as_posix() + ".tbi") and not (vcf_path.with_suffix(vcf_path.suffix + ".tbi").exists()):
        # best-effort check
        if not (Path(str(vcf_path) + ".tbi").exists()):
            raise FileNotFoundError(f"Tabix index (.tbi) not found for: {vcf_path}")

    _ensure_parent(db_path)

    if force and db_path.exists():
        db_path.unlink()

    conn = sqlite3.connect(str(db_path))
    try:
        _init_db(conn)
        vcf = pysam.VariantFile(str(vcf_path))

        available = list(vcf.header.contigs)
        if chroms:
            targets = [c for c in chroms if c in available]
            missing = [c for c in chroms if c not in available]
            if missing:
                print(f"[warn] requested chroms not in VCF header and will be skipped: {missing}")
        else:
            targets = available

        c = conn.cursor()
        for chrom in targets:
            resume_pos = _get_resume_pos(conn, chrom)
            start0 = max(0, (resume_pos or 0) - 1) if resume_pos else None  # fetch uses 0-based
            label = f"{chrom} (resume after pos={resume_pos})" if resume_pos else chrom
            print(f"[chrom] {label}")

            last_seen = resume_pos or 0
            inserted = 0

            rec_iter = _iter_records(vcf, chrom, start_pos=start0)
            for batch in tqdm(_batched(rec_iter, batch_size), desc=f"ingest:{chrom}", unit="batch"):
                # executemany with PRIMARY KEY: ignore duplicates
                c.executemany(
                    "INSERT OR IGNORE INTO variants (chrom, pos, rsid) VALUES (?, ?, ?)",
                    batch,
                )
                conn.commit()
                inserted += c.rowcount if c.rowcount != -1 else 0
                last_seen = max(last_seen, max(p for _, p, _ in batch))
                _set_resume_pos(conn, chrom, last_seen)

            print(f"[ok] {chrom}: last_pos={last_seen} inserted~={inserted}")

        print(f"[done] DB written: {db_path}")
    finally:
        conn.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="Path to bgzipped dbSNP VCF (.vcf.gz) with a .tbi index.")
    ap.add_argument("--db", required=True, help="Output SQLite DB path.")
    ap.add_argument("--batch-size", type=int, default=100000, help="Insert batch size.")
    ap.add_argument("--chroms", default=None, help="Comma-separated chromosome list (default: all contigs).")
    ap.add_argument("--force", action="store_true", help="Overwrite DB if it already exists.")
    args = ap.parse_args()

    chroms = [c.strip() for c in args.chroms.split(",")] if args.chroms else None
    build_db(Path(args.vcf), Path(args.db), int(args.batch_size), chroms, bool(args.force))


if __name__ == "__main__":
    main()
