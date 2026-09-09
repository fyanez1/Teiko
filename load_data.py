#!/usr/bin/env python3
"""Part 1 - initialise the SQLite schema and load ``cell-count.csv``.

Run from the repository root with no arguments::

    python load_data.py

It (re)creates ``cell_counts.db`` next to this file. The schema lives in
``cellcount/schema.sql``; this script only depends on the Python standard
library so it works before any third-party packages are installed.
"""
from __future__ import annotations

import csv
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"
SCHEMA_PATH = ROOT / "cellcount" / "schema.sql"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
METADATA_COLUMNS = [
    "project", "subject", "condition", "age", "sex", "treatment", "response",
    "sample", "sample_type", "time_from_treatment_start",
]
EXPECTED_COLUMNS = METADATA_COLUMNS + POPULATIONS
TABLES = ["projects", "subjects", "samples", "cell_populations", "cell_counts"]


def _nullable(value: str | None) -> str | None:
    """Map empty CSV cells (e.g. response for healthy subjects) to SQL NULL."""
    if value is None:
        return None
    value = value.strip()
    return value if value else None


def read_rows(csv_path: Path = CSV_PATH) -> list[dict]:
    """Read and validate the CSV, returning one dict per row."""
    if not csv_path.exists():
        raise FileNotFoundError(f"Input file not found: {csv_path}")
    with open(csv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in EXPECTED_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{csv_path.name} is missing expected columns: {missing}")
        rows = [{k: (v.strip() if isinstance(v, str) else v) for k, v in row.items()} for row in reader]
    if not rows:
        raise ValueError(f"{csv_path.name} contains no data rows")
    return rows


def build_database(csv_path: Path = CSV_PATH, db_path: Path = DB_PATH) -> dict[str, int]:
    """Create a fresh database at ``db_path`` from ``csv_path``.

    Returns a mapping of table name -> number of rows loaded.
    """
    rows = read_rows(Path(csv_path))
    db_path = Path(db_path)
    if db_path.exists():
        db_path.unlink()  # always rebuild from scratch so the load is reproducible

    conn = sqlite3.connect(str(db_path))
    try:
        conn.execute("PRAGMA foreign_keys = ON")
        conn.executescript(SCHEMA_PATH.read_text(encoding="utf-8"))

        # --- lookup table of populations ------------------------------------
        conn.executemany(
            "INSERT INTO cell_populations (name, display_order) VALUES (?, ?)",
            [(name, i) for i, name in enumerate(POPULATIONS)],
        )

        # --- projects ---------------------------------------------------------
        projects = sorted({r["project"] for r in rows})
        conn.executemany("INSERT INTO projects (project_id) VALUES (?)", [(p,) for p in projects])

        # --- subjects (validate that subject-level metadata is consistent) ----
        subjects: dict[str, tuple] = {}
        for r in rows:
            attrs = (
                r["project"], r["condition"], int(r["age"]), r["sex"],
                r["treatment"], _nullable(r["response"]),
            )
            previous = subjects.setdefault(r["subject"], attrs)
            if previous != attrs:
                raise ValueError(
                    f"Subject {r['subject']} has conflicting metadata across rows: "
                    f"{previous} vs {attrs}"
                )
        conn.executemany(
            "INSERT INTO subjects (subject_id, project_id, condition, age, sex, treatment, response)"
            " VALUES (?, ?, ?, ?, ?, ?, ?)",
            [(subject_id, *attrs) for subject_id, attrs in subjects.items()],
        )

        # --- samples ----------------------------------------------------------
        conn.executemany(
            "INSERT INTO samples (sample_id, subject_id, sample_type, time_from_treatment_start)"
            " VALUES (?, ?, ?, ?)",
            [(r["sample"], r["subject"], r["sample_type"], int(r["time_from_treatment_start"])) for r in rows],
        )

        # --- cell counts in long form ------------------------------------------
        population_ids = dict(conn.execute("SELECT name, population_id FROM cell_populations"))
        conn.executemany(
            "INSERT INTO cell_counts (sample_id, population_id, count) VALUES (?, ?, ?)",
            [(r["sample"], population_ids[p], int(r[p])) for r in rows for p in POPULATIONS],
        )
        conn.commit()

        counts = {t: conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0] for t in TABLES}
    finally:
        conn.close()
    return counts


def main() -> int:
    try:
        counts = build_database()
    except (FileNotFoundError, ValueError, sqlite3.DatabaseError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"Loaded {CSV_PATH.name} into {DB_PATH.name}")
    for table, n in counts.items():
        print(f"  {table:<17} {n:>7,} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
