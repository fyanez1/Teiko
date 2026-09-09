"""Part 2 - relative frequency of each cell population in each sample."""
from __future__ import annotations

import sqlite3

import pandas as pd

from .db import query

SUMMARY_SQL = """
SELECT sample, total_count, population, count, percentage
FROM sample_frequencies
ORDER BY sample, display_order
"""


def summary_table(conn: sqlite3.Connection) -> pd.DataFrame:
    """One row per (sample, population) with the population's share of the sample.

    Columns: sample, total_count, population, count, percentage.
    ``total_count`` is the sum of the five population counts for the sample and
    ``percentage`` is ``100 * count / total_count``. Both are computed by the
    ``sample_frequencies`` SQL view so every consumer uses the same definition.
    """
    return query(conn, SUMMARY_SQL)


def summary_for_sample(conn: sqlite3.Connection, sample_id: str) -> pd.DataFrame:
    """Summary rows for a single sample (used by the dashboard)."""
    sql = SUMMARY_SQL.replace("FROM sample_frequencies", "FROM sample_frequencies WHERE sample = ?")
    return query(conn, sql, (sample_id,))
