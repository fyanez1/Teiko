"""Database helpers shared by the analysis modules, tests and the dashboard."""
from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Iterable

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"
SCHEMA_PATH = Path(__file__).resolve().parent / "schema.sql"
OUTPUT_DIR = ROOT / "outputs"

# Canonical population order (matches the CSV column order and `display_order`).
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def connect(db_path: Path | str = DB_PATH) -> sqlite3.Connection:
    """Open a connection with foreign-key enforcement switched on."""
    db_path = Path(db_path)
    if not db_path.exists():
        raise FileNotFoundError(
            f"Database not found at {db_path}. Run `python load_data.py` (or `make pipeline`) first."
        )
    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def query(conn: sqlite3.Connection, sql: str, params: Iterable = ()) -> pd.DataFrame:
    """Run a parameterised SELECT and return the result as a DataFrame."""
    return pd.read_sql_query(sql, conn, params=tuple(params))
