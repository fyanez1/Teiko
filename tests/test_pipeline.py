"""End-to-end checks: loader, schema views and the Part 2-4 analysis functions.

Every expectation is recomputed independently with pandas from the raw CSV so
the SQL logic is verified against a second implementation.
"""
import numpy as np
import pandas as pd
import pytest

from cellcount.db import CSV_PATH, POPULATIONS, connect
from cellcount.stats import bh_adjust, compare_groups, responder_dataset, subject_level_dataset
from cellcount.subsets import (
    average_population_count,
    baseline_samples,
    samples_per_project,
    subjects_by_response,
    subjects_by_sex,
)
from cellcount.summary import summary_table
from load_data import build_database


@pytest.fixture(scope="module")
def raw() -> pd.DataFrame:
    return pd.read_csv(CSV_PATH)


@pytest.fixture(scope="module")
def conn(tmp_path_factory):
    db_path = tmp_path_factory.mktemp("db") / "test.db"
    build_database(CSV_PATH, db_path)
    connection = connect(db_path)
    yield connection
    connection.close()


def test_loader_row_counts(conn, raw):
    counts = {t: conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
              for t in ["projects", "subjects", "samples", "cell_populations", "cell_counts"]}
    assert counts["projects"] == raw["project"].nunique()
    assert counts["subjects"] == raw["subject"].nunique()
    assert counts["samples"] == len(raw)
    assert counts["cell_populations"] == len(POPULATIONS)
    assert counts["cell_counts"] == len(raw) * len(POPULATIONS)


def test_loader_maps_blank_response_to_null(conn, raw):
    n_null = conn.execute("SELECT COUNT(*) FROM subjects WHERE response IS NULL").fetchone()[0]
    expected = raw.drop_duplicates("subject")["response"].isna().sum()
    assert n_null == expected


def test_loader_rejects_inconsistent_subject_metadata(tmp_path):
    bad = pd.read_csv(CSV_PATH).head(2).copy()
    bad.loc[1, "sex"] = "F"  # same subject, conflicting sex
    csv_path = tmp_path / "bad.csv"
    bad.to_csv(csv_path, index=False)
    with pytest.raises(ValueError, match="conflicting metadata"):
        build_database(csv_path, tmp_path / "bad.db")


def test_summary_table_matches_pandas(conn, raw):
    summary = summary_table(conn)
    assert list(summary.columns) == ["sample", "total_count", "population", "count", "percentage"]
    assert len(summary) == len(raw) * len(POPULATIONS)

    totals = summary.groupby("sample")["percentage"].sum()
    assert np.allclose(totals, 100.0)

    expected_total = raw.set_index("sample")[POPULATIONS].sum(axis=1)
    got_total = summary.drop_duplicates("sample").set_index("sample")["total_count"]
    pd.testing.assert_series_equal(got_total.sort_index(), expected_total.sort_index(),
                                   check_names=False, check_dtype=False)

    row = raw.iloc[123]
    got = summary[(summary["sample"] == row["sample"]) & (summary["population"] == "nk_cell")].iloc[0]
    assert got["count"] == row["nk_cell"]
    assert got["percentage"] == pytest.approx(100 * row["nk_cell"] / row[POPULATIONS].sum())


def test_responder_cohort_and_stats(conn, raw):
    df = responder_dataset(conn)
    expected = raw[(raw["condition"] == "melanoma") & (raw["treatment"] == "miraclib")
                   & (raw["sample_type"] == "PBMC") & raw["response"].isin(["yes", "no"])]
    assert df["sample"].nunique() == len(expected)
    assert set(df["response"]) == {"yes", "no"}

    stats = compare_groups(df)
    assert list(stats["population"]) == POPULATIONS
    assert ((stats["mannwhitney_p"] >= 0) & (stats["mannwhitney_p"] <= 1)).all()
    assert (stats["mannwhitney_p_adj"] >= stats["mannwhitney_p"] - 1e-12).all()
    assert stats["cliffs_delta"].between(-1, 1).all()

    # Cross-check one Mann-Whitney p-value with scipy on pandas-filtered data.
    from scipy.stats import mannwhitneyu
    pct = expected.assign(total=expected[POPULATIONS].sum(axis=1))
    pct["cd4_pct"] = 100 * pct["cd4_t_cell"] / pct["total"]
    p = mannwhitneyu(pct.loc[pct["response"] == "yes", "cd4_pct"],
                     pct.loc[pct["response"] == "no", "cd4_pct"], alternative="two-sided").pvalue
    assert stats.set_index("population").loc["cd4_t_cell", "mannwhitney_p"] == pytest.approx(p)

    subject_level = subject_level_dataset(df)
    assert len(subject_level) == df["subject"].nunique() * len(POPULATIONS)


def test_bh_adjust():
    p = np.array([0.01, 0.04, 0.03, 0.20, 0.50])
    adj = bh_adjust(p)
    # Step-up procedure: sorted p * n / rank = [0.05, 0.075, 0.0667, 0.25, 0.5],
    # then enforce monotonicity from the largest p downwards.
    assert np.allclose(adj, [0.05, 0.2 / 3, 0.2 / 3, 0.25, 0.50])
    assert (adj >= p).all() and (adj <= 1).all()


def test_baseline_subsets(conn, raw):
    expected = raw[(raw["condition"] == "melanoma") & (raw["treatment"] == "miraclib")
                   & (raw["sample_type"] == "PBMC") & (raw["time_from_treatment_start"] == 0)]
    cohort = baseline_samples(conn)
    assert len(cohort) == len(expected)

    per_project = samples_per_project(conn).set_index("project")["n_samples"]
    assert set(per_project.index) == set(raw["project"].unique())  # every project listed, even with 0
    for project, n in expected.groupby("project").size().items():
        assert per_project[project] == n

    by_response = subjects_by_response(conn).set_index("response")["n_subjects"]
    assert by_response["responder"] == expected.loc[expected["response"] == "yes", "subject"].nunique()
    assert by_response["non-responder"] == expected.loc[expected["response"] == "no", "subject"].nunique()

    by_sex = subjects_by_sex(conn).set_index("sex")["n_subjects"]
    assert by_sex["male"] == expected.loc[expected["sex"] == "M", "subject"].nunique()
    assert by_sex["female"] == expected.loc[expected["sex"] == "F", "subject"].nunique()


def test_average_b_cells_melanoma_male_responders_baseline(conn, raw):
    expected = raw[(raw["condition"] == "melanoma") & (raw["sex"] == "M")
                   & (raw["response"] == "yes") & (raw["time_from_treatment_start"] == 0)]
    avg, n = average_population_count(conn)
    assert n == len(expected)
    assert avg == pytest.approx(expected["b_cell"].mean())
    assert round(avg, 2) == round(expected["b_cell"].mean(), 2)


def test_average_population_count_handles_empty_selection(conn):
    avg, n = average_population_count(conn, condition="no-such-condition")
    assert avg is None and n == 0
