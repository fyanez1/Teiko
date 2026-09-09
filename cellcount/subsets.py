"""Part 4 - subset queries for early-treatment (baseline) samples.

All functions issue parameterised SQL against the database so that the same
code paths serve the pipeline, the tests and the dashboard.
"""
from __future__ import annotations

import sqlite3

import pandas as pd

from .db import query

# Filters describing the baseline cohort in the brief.
BASELINE_FILTERS = {
    "condition": "melanoma",
    "treatment": "miraclib",
    "sample_type": "PBMC",
    "time_from_treatment_start": 0,
}

_FILTERABLE = (
    "project_id", "condition", "treatment", "sample_type",
    "time_from_treatment_start", "response", "sex",
)


def _conditions(filters: dict) -> tuple[list[str], list]:
    """Translate a {column: value} mapping into SQL predicates on alias ``m``."""
    clauses, params = [], []
    for column, value in filters.items():
        if column not in _FILTERABLE:
            raise ValueError(f"Unsupported filter column: {column}")
        if value is None:  # None means "do not filter on this column"
            continue
        clauses.append(f"m.{column} = ?")
        params.append(value)
    return clauses, params


def _where(filters: dict) -> tuple[str, list]:
    clauses, params = _conditions(filters)
    return (" WHERE " + " AND ".join(clauses)) if clauses else "", params


def baseline_samples(conn: sqlite3.Connection, **overrides) -> pd.DataFrame:
    """All samples in the baseline cohort (melanoma, miraclib, PBMC, day 0 by default)."""
    filters = {**BASELINE_FILTERS, **overrides}
    where, params = _where(filters)
    sql = f"""
    SELECT m.sample_id AS sample, m.subject_id AS subject, m.project_id AS project,
           m.condition, m.treatment, m.sample_type, m.time_from_treatment_start,
           m.response, m.sex, m.age
    FROM sample_metadata AS m{where}
    ORDER BY m.project_id, m.subject_id, m.sample_id
    """
    return query(conn, sql, params)


def samples_per_project(conn: sqlite3.Connection, **overrides) -> pd.DataFrame:
    """Number of cohort samples contributed by each project (projects with none show 0)."""
    filters = {**BASELINE_FILTERS, **overrides}
    clauses, params = _conditions(filters)
    join_on = " AND ".join(["m.project_id = p.project_id", *clauses])
    sql = f"""
    SELECT p.project_id AS project, COUNT(m.sample_id) AS n_samples
    FROM projects AS p
    LEFT JOIN sample_metadata AS m ON {join_on}
    GROUP BY p.project_id
    ORDER BY p.project_id
    """
    return query(conn, sql, params)


def subjects_by_response(conn: sqlite3.Connection, **overrides) -> pd.DataFrame:
    """Number of distinct cohort subjects who responded / did not respond."""
    filters = {**BASELINE_FILTERS, **overrides}
    where, params = _where(filters)
    sql = f"""
    SELECT CASE m.response WHEN 'yes' THEN 'responder'
                           WHEN 'no'  THEN 'non-responder'
                           ELSE 'unknown' END AS response,
           COUNT(DISTINCT m.subject_id) AS n_subjects
    FROM sample_metadata AS m{where}
    GROUP BY m.response
    ORDER BY m.response DESC
    """
    return query(conn, sql, params)


def subjects_by_sex(conn: sqlite3.Connection, **overrides) -> pd.DataFrame:
    """Number of distinct cohort subjects by sex."""
    filters = {**BASELINE_FILTERS, **overrides}
    where, params = _where(filters)
    sql = f"""
    SELECT CASE m.sex WHEN 'M' THEN 'male' WHEN 'F' THEN 'female' ELSE 'unknown' END AS sex,
           COUNT(DISTINCT m.subject_id) AS n_subjects
    FROM sample_metadata AS m{where}
    GROUP BY m.sex
    ORDER BY m.sex
    """
    return query(conn, sql, params)


def average_population_count(
    conn: sqlite3.Connection,
    population: str = "b_cell",
    condition: str | None = "melanoma",
    sex: str | None = "M",
    response: str | None = "yes",
    time_from_treatment_start: int | None = 0,
    treatment: str | None = None,
    sample_type: str | None = None,
) -> tuple[float | None, int]:
    """Average raw count of ``population`` over the samples matching the filters.

    The defaults answer the brief's question: melanoma males, responders, at
    time 0, across *all* sample types and treatments (``None`` = no filter).
    Returns ``(average, n_samples)``; the average is ``None`` if nothing matches.
    """
    filters = {
        "condition": condition, "sex": sex, "response": response,
        "time_from_treatment_start": time_from_treatment_start,
        "treatment": treatment, "sample_type": sample_type,
    }
    clauses, params = _conditions(filters)
    clauses.append("p.name = ?")
    params.append(population)
    sql = f"""
    SELECT AVG(c.count) AS avg_count, COUNT(*) AS n_samples
    FROM cell_counts AS c
    JOIN cell_populations AS p ON p.population_id = c.population_id
    JOIN sample_metadata  AS m ON m.sample_id     = c.sample_id
    WHERE {" AND ".join(clauses)}
    """
    row = query(conn, sql, params).iloc[0]
    avg = None if pd.isna(row["avg_count"]) else float(row["avg_count"])
    return avg, int(row["n_samples"])


def melanoma_male_responder_baseline_b_cells(conn: sqlite3.Connection) -> tuple[float | None, int]:
    """Average B-cell count for melanoma male responders at time 0 (all sample & treatment types)."""
    return average_population_count(conn)
