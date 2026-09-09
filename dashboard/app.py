"""Interactive dashboard presenting Parts 2-4 from ``cell_counts.db`` (Streamlit).

Start it with ``make dashboard`` (equivalent to
``python -m streamlit run dashboard/app.py``). Every number shown here is
computed live from the SQLite database created by the pipeline, using the same
``cellcount`` library functions as ``run_analysis.py``.
"""
from __future__ import annotations

import sys
from contextlib import closing
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:  # `streamlit run dashboard/app.py` puts dashboard/ first on sys.path
    sys.path.insert(0, str(ROOT))

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from cellcount import DB_PATH, POPULATIONS, connect, query
from cellcount.stats import (
    ALPHA,
    GROUP_COLOURS,
    compare_by_timepoint,
    compare_groups,
    responder_dataset,
    subject_level_dataset,
)
from cellcount.subsets import (
    BASELINE_FILTERS,
    average_population_count,
    baseline_samples,
    samples_per_project,
    subjects_by_response,
    subjects_by_sex,
)
from cellcount.summary import summary_table

st.set_page_config(page_title="Loblaw Bio - Immune Cell Dashboard", page_icon="🧬", layout="wide")

STATS_DISPLAY = [
    "population", "n_responders", "n_non_responders", "median_responders", "median_non_responders",
    "median_difference", "cliffs_delta", "mannwhitney_p", "mannwhitney_p_adj", "welch_p", "welch_p_adj",
    "significant", "direction",
]


# --------------------------------------------------------------------------- data access
def ensure_database() -> None:
    """Build the database on first use so the app also works on a fresh clone/deploy."""
    if DB_PATH.exists():
        return
    import load_data  # root-level loader (stdlib only)

    with st.spinner("Database not found - building it from cell-count.csv ..."):
        load_data.build_database()
    st.cache_data.clear()


@st.cache_data(show_spinner=False)
def db_overview() -> dict:
    with closing(connect()) as conn:
        return {
            "projects": conn.execute("SELECT COUNT(*) FROM projects").fetchone()[0],
            "subjects": conn.execute("SELECT COUNT(*) FROM subjects").fetchone()[0],
            "samples": conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0],
            "cell_counts": conn.execute("SELECT COUNT(*) FROM cell_counts").fetchone()[0],
        }


@st.cache_data(show_spinner=False)
def load_summary_with_metadata() -> pd.DataFrame:
    with closing(connect()) as conn:
        summary = summary_table(conn)
        meta = query(conn, "SELECT sample_id AS sample, subject_id AS subject, project_id AS project, "
                           "condition, treatment, response, sex, age, sample_type, "
                           "time_from_treatment_start FROM sample_metadata")
    return summary.merge(meta, on="sample", how="left")


@st.cache_data(show_spinner=False)
def cohort_options() -> dict:
    with closing(connect()) as conn:
        return {
            "condition": query(conn, "SELECT DISTINCT condition FROM subjects WHERE response IS NOT NULL "
                                     "ORDER BY condition")["condition"].tolist(),
            "treatment": query(conn, "SELECT DISTINCT treatment FROM subjects WHERE response IS NOT NULL "
                                     "ORDER BY treatment")["treatment"].tolist(),
            "sample_type": query(conn, "SELECT DISTINCT sample_type FROM samples ORDER BY sample_type")["sample_type"].tolist(),
            "time": query(conn, "SELECT DISTINCT time_from_treatment_start AS t FROM samples ORDER BY t")["t"].tolist(),
            "all_conditions": query(conn, "SELECT DISTINCT condition FROM subjects ORDER BY condition")["condition"].tolist(),
            "all_treatments": query(conn, "SELECT DISTINCT treatment FROM subjects ORDER BY treatment")["treatment"].tolist(),
            "projects": query(conn, "SELECT project_id FROM projects ORDER BY project_id")["project_id"].tolist(),
        }


@st.cache_data(show_spinner=False)
def load_responders(condition: str, treatment: str, sample_type: str) -> pd.DataFrame:
    with closing(connect()) as conn:
        return responder_dataset(conn, condition, treatment, sample_type)


@st.cache_data(show_spinner=False)
def load_baseline() -> dict:
    with closing(connect()) as conn:
        avg, n = average_population_count(conn)
        return {
            "samples": baseline_samples(conn),
            "per_project": samples_per_project(conn),
            "by_response": subjects_by_response(conn),
            "by_sex": subjects_by_sex(conn),
            "avg_b_cells": avg,
            "avg_b_cells_n": n,
        }


@st.cache_data(show_spinner=False)
def explore_average(population, condition, sex, response, time_point, treatment, sample_type):
    with closing(connect()) as conn:
        return average_population_count(conn, population, condition, sex, response, time_point,
                                        treatment, sample_type)


def quantile_boxplot(df: pd.DataFrame, group_col: str) -> go.Figure:
    """Boxplot from precomputed quantiles (whiskers at the 5th/95th percentiles).

    Sends five numbers per box to the browser instead of every sample.
    """
    q = (df.groupby([group_col, "population"])["percentage"]
           .quantile([0.05, 0.25, 0.5, 0.75, 0.95]).unstack())
    fig = go.Figure()
    for group, sub in q.groupby(level=0):
        sub = sub.droplevel(0)  # index is now the population name
        sub = sub.reindex([pop for pop in POPULATIONS if pop in sub.index])
        fig.add_trace(go.Box(name=str(group), x=sub.index.tolist(), lowerfence=sub[0.05], q1=sub[0.25],
                             median=sub[0.5], q3=sub[0.75], upperfence=sub[0.95]))
    fig.update_layout(boxmode="group", height=380, margin=dict(t=20, b=20),
                      yaxis_title="% of total cells", legend_title_text=group_col)
    return fig


def fmt_p(p: float) -> str:
    return "n/a" if pd.isna(p) else (f"{p:.2e}" if p < 1e-3 else f"{p:.3f}")


def option_or_all(label: str, options: list, key: str):
    """Selectbox with an explicit 'all' entry that maps to None (no filter)."""
    choice = st.selectbox(label, ["(all)"] + list(options), key=key)
    return None if choice == "(all)" else choice


# --------------------------------------------------------------------------- layout
ensure_database()
overview = db_overview()
options = cohort_options()

st.title("🧬 Loblaw Bio - immune cell population dashboard")
st.caption(
    f"Source: `{DB_PATH.name}` - {overview['projects']} projects, {overview['subjects']:,} subjects, "
    f"{overview['samples']:,} samples, {overview['cell_counts']:,} cell-count records. "
    "All figures are computed live from the database built by `make pipeline`."
)

with st.sidebar:
    st.header("About")
    st.markdown(
        """
        **Part 2** - relative frequency of each immune cell population per sample.

        **Part 3** - responders vs non-responders (melanoma, miraclib, PBMC):
        boxplots, Mann-Whitney U / Welch tests, effect sizes, FDR correction.

        **Part 4** - baseline (day 0) melanoma PBMC samples from miraclib-treated
        patients, plus the B-cell question.
        """
    )
    st.markdown("---")
    st.markdown("**Populations**: " + ", ".join(f"`{p}`" for p in POPULATIONS))
    st.markdown("Reproduce everything with `make pipeline`; outputs are written to `outputs/`.")

tab2, tab3, tab4 = st.tabs([
    "Part 2 · Cell frequencies", "Part 3 · Responders vs non-responders", "Part 4 · Baseline subsets",
])

# ---------------------------------------------------------------- Part 2
with tab2:
    st.subheader("Relative frequency of each cell population in each sample")
    st.markdown(
        "For every sample, `total_count` is the sum of the five population counts and "
        "`percentage = 100 × count / total_count`. Use the filters to narrow the table; "
        "the download always contains the filtered rows."
    )
    data = load_summary_with_metadata()

    f1, f2, f3, f4, f5, f6 = st.columns(6)
    with f1:
        sel_project = st.multiselect("Project", options["projects"], default=options["projects"])
    with f2:
        sel_condition = st.multiselect("Condition", options["all_conditions"], default=options["all_conditions"])
    with f3:
        sel_treatment = st.multiselect("Treatment", options["all_treatments"], default=options["all_treatments"])
    with f4:
        sel_type = st.multiselect("Sample type", options["sample_type"], default=options["sample_type"])
    with f5:
        sel_time = st.multiselect("Time point (days)", options["time"], default=options["time"])
    with f6:
        sel_pop = st.multiselect("Population", POPULATIONS, default=POPULATIONS)

    filtered = data[
        data["project"].isin(sel_project) & data["condition"].isin(sel_condition)
        & data["treatment"].isin(sel_treatment) & data["sample_type"].isin(sel_type)
        & data["time_from_treatment_start"].isin(sel_time) & data["population"].isin(sel_pop)
    ]
    search = st.text_input("Search by sample or subject id (e.g. `sample00012` or `sbj004`)", "")
    if search:
        filtered = filtered[filtered["sample"].str.contains(search, case=False)
                            | filtered["subject"].str.contains(search, case=False)]

    m1, m2, m3 = st.columns(3)
    m1.metric("Samples shown", f"{filtered['sample'].nunique():,}")
    m2.metric("Rows shown", f"{len(filtered):,}")
    m3.metric("Mean total cells / sample", f"{filtered.drop_duplicates('sample')['total_count'].mean():,.0f}"
              if len(filtered) else "-")

    summary_cols = ["sample", "total_count", "population", "count", "percentage"]
    show_meta = st.checkbox("Show sample metadata columns", value=False)
    table = filtered if show_meta else filtered[summary_cols]
    st.dataframe(table, width="stretch", height=420, hide_index=True,
                 column_config={"percentage": st.column_config.NumberColumn(format="%.2f")})
    st.download_button("Download filtered summary table (CSV)", table.to_csv(index=False).encode(),
                       "summary_table_filtered.csv", "text/csv")

    if len(filtered):
        c1, c2 = st.columns(2)
        with c1:
            st.markdown("**Average composition by condition** (filtered samples)")
            comp = (filtered.groupby(["condition", "population"], as_index=False)["percentage"].mean())
            fig = px.bar(comp, x="condition", y="percentage", color="population", barmode="stack",
                         category_orders={"population": POPULATIONS},
                         labels={"percentage": "Mean % of total cells"})
            fig.update_layout(height=380, margin=dict(t=20, b=20))
            st.plotly_chart(fig)
        with c2:
            st.markdown("**Distribution of relative frequencies** (filtered samples)")
            fig = quantile_boxplot(filtered, "sample_type")
            st.plotly_chart(fig)

# ---------------------------------------------------------------- Part 3
with tab3:
    st.subheader("Do responders and non-responders differ in cell population frequencies?")
    st.markdown(
        "Default cohort (per the brief): **melanoma** patients treated with **miraclib**, **PBMC** samples "
        "only. Responders (`response = yes`) are compared with non-responders (`response = no`) for each "
        "population using the per-sample relative frequencies from Part 2."
    )
    c1, c2, c3, c4, c5 = st.columns(5)
    with c1:
        condition = st.selectbox("Condition", options["condition"],
                                 index=options["condition"].index("melanoma") if "melanoma" in options["condition"] else 0)
    with c2:
        treatment = st.selectbox("Treatment", options["treatment"],
                                 index=options["treatment"].index("miraclib") if "miraclib" in options["treatment"] else 0)
    with c3:
        sample_type = st.selectbox("Sample type", options["sample_type"],
                                   index=options["sample_type"].index("PBMC") if "PBMC" in options["sample_type"] else 0)
    with c4:
        unit = st.radio("Unit of analysis", ["Samples (all time points)", "Subjects (mean per subject)"],
                        help="Each subject contributes several samples. The subject-level option averages them "
                             "first so repeated measures from one person are not treated as independent.")
    with c5:
        alpha = st.select_slider("Significance level α", options=[0.01, 0.05, 0.10], value=ALPHA)

    df = load_responders(condition, treatment, sample_type)
    if df.empty:
        st.warning("No samples with a recorded response match this cohort.")
    else:
        time_points = sorted(df["time_from_treatment_start"].unique())
        sel_times = st.multiselect("Time points included (days from treatment start)", time_points,
                                   default=time_points)
        df = df[df["time_from_treatment_start"].isin(sel_times)]
        analysis_df = subject_level_dataset(df) if unit.startswith("Subjects") else df
        stats = compare_groups(analysis_df, alpha=alpha)

        n_r = df.loc[df["response"] == "yes", "sample"].nunique()
        n_nr = df.loc[df["response"] == "no", "sample"].nunique()
        s_r = df.loc[df["response"] == "yes", "subject"].nunique()
        s_nr = df.loc[df["response"] == "no", "subject"].nunique()
        k1, k2, k3, k4 = st.columns(4)
        k1.metric("Responder samples", n_r, f"{s_r} subjects", delta_color="off")
        k2.metric("Non-responder samples", n_nr, f"{s_nr} subjects", delta_color="off")
        sig = stats.loc[stats["significant"], "population"].tolist()
        nominal = stats.loc[stats["mannwhitney_p"] < alpha, "population"].tolist()
        k3.metric("Significant after FDR (BH)", len(sig), ", ".join(sig) if sig else "none", delta_color="off")
        k4.metric(f"Nominal p < {alpha}", len(nominal), ", ".join(nominal) if nominal else "none",
                  delta_color="off")

        st.markdown("#### Boxplots - one panel per population")
        plot_df = analysis_df.copy()
        fig = px.box(
            plot_df, x="group", y="percentage", color="group", facet_col="population",
            points="all", color_discrete_map=GROUP_COLOURS,
            category_orders={"population": POPULATIONS, "group": ["responder", "non-responder"]},
            labels={"percentage": "% of total cells", "group": ""},
            hover_data=[c for c in ["sample", "subject", "time_from_treatment_start"] if c in plot_df.columns],
        )
        fig.update_traces(marker=dict(size=3, opacity=0.35), jitter=0.4)
        fig.update_yaxes(matches=None, showticklabels=True)
        fig.update_xaxes(showticklabels=False)
        lookup = stats.set_index("population")
        fig.for_each_annotation(lambda a: a.update(text=(
            f"{a.text.split('=')[-1]}{' *' if bool(lookup.loc[a.text.split('=')[-1], 'significant']) else ''}"
            f"<br><sup>p={fmt_p(lookup.loc[a.text.split('=')[-1], 'mannwhitney_p'])}, "
            f"q={fmt_p(lookup.loc[a.text.split('=')[-1], 'mannwhitney_p_adj'])}, "
            f"δ={lookup.loc[a.text.split('=')[-1], 'cliffs_delta']:+.2f}</sup>")))
        fig.update_layout(height=480, legend_title_text="", margin=dict(t=70, b=20))
        st.plotly_chart(fig)
        st.caption("p = two-sided Mann-Whitney U; q = Benjamini-Hochberg adjusted p across the five populations; "
                   "δ = Cliff's delta (> 0: higher in responders). * marks q < α.")

        st.markdown("#### Statistics per population")
        st.dataframe(
            stats[STATS_DISPLAY].style.format({
                "median_responders": "{:.2f}", "median_non_responders": "{:.2f}", "median_difference": "{:+.2f}",
                "cliffs_delta": "{:+.3f}", "mannwhitney_p": "{:.4f}", "mannwhitney_p_adj": "{:.4f}",
                "welch_p": "{:.4f}", "welch_p_adj": "{:.4f}",
            }),
            width="stretch", hide_index=True,
        )

        st.markdown("#### Conclusion")
        strongest = stats.loc[stats["mannwhitney_p"].idxmin()]
        if sig:
            st.success(
                f"After Benjamini-Hochberg correction (q < {alpha}) the following populations differ between "
                f"responders and non-responders: **{', '.join(sig)}**. "
                + " ".join(
                    f"{r['population']} is {r['direction']} (median {r['median_responders']:.2f}% vs "
                    f"{r['median_non_responders']:.2f}%, q = {fmt_p(r['mannwhitney_p_adj'])}, "
                    f"Cliff's δ = {r['cliffs_delta']:+.3f})."
                    for _, r in stats[stats["significant"]].iterrows())
            )
        else:
            st.info(
                f"No population differs significantly once the five tests are corrected for multiple "
                f"comparisons (all q ≥ {alpha}). "
                + (f"Nominally significant before correction: **{', '.join(nominal)}**. " if nominal else "")
                + f"The strongest candidate is **{strongest['population']}** "
                f"({'higher' if strongest['cliffs_delta'] > 0 else 'lower'} in responders; "
                f"p = {fmt_p(strongest['mannwhitney_p'])}, q = {fmt_p(strongest['mannwhitney_p_adj'])}, "
                f"Cliff's δ = {strongest['cliffs_delta']:+.3f}), a small effect that would need confirmation "
                "in an independent cohort before being used to predict response."
            )

        with st.expander("Sensitivity analysis: each time point tested separately"):
            by_time = compare_by_timepoint(df, alpha=alpha)
            st.dataframe(
                by_time[["time_from_treatment_start", *STATS_DISPLAY]].style.format({
                    "median_responders": "{:.2f}", "median_non_responders": "{:.2f}",
                    "median_difference": "{:+.2f}", "cliffs_delta": "{:+.3f}", "mannwhitney_p": "{:.4f}",
                    "mannwhitney_p_adj": "{:.4f}", "welch_p": "{:.4f}", "welch_p_adj": "{:.4f}",
                }),
                width="stretch", hide_index=True,
            )
        with st.expander("Methods"):
            st.markdown(
                """
                * **Data**: per-sample relative frequencies (Part 2) for the selected cohort; only samples
                  with a recorded response are used.
                * **Primary test**: two-sided Mann-Whitney U on percentages (no normality assumption).
                  **Companion**: Welch's t-test. **Effect size**: Cliff's delta, derived from U.
                * **Multiple testing**: Benjamini-Hochberg across the five populations; a population is
                  called significant when the adjusted Mann-Whitney p is below α.
                * **Repeated measures**: subjects contribute three samples (days 0, 7, 14). Switch the unit
                  of analysis to *Subjects* to average them first, or restrict the time points.
                """
            )

# ---------------------------------------------------------------- Part 4
with tab4:
    st.subheader("Baseline melanoma PBMC samples from miraclib-treated patients")
    st.markdown(
        f"Filters: condition = **{BASELINE_FILTERS['condition']}**, treatment = **{BASELINE_FILTERS['treatment']}**, "
        f"sample type = **{BASELINE_FILTERS['sample_type']}**, time from treatment start = "
        f"**{BASELINE_FILTERS['time_from_treatment_start']}** days."
    )
    base = load_baseline()
    cohort = base["samples"]
    k1, k2, k3 = st.columns(3)
    k1.metric("Baseline samples", len(cohort))
    k2.metric("Subjects", cohort["subject"].nunique())
    k3.metric("Projects contributing", int((base["per_project"]["n_samples"] > 0).sum()))

    c1, c2, c3 = st.columns(3)
    with c1:
        st.markdown("**Samples per project**")
        st.dataframe(base["per_project"], width="stretch", hide_index=True)
        fig = px.bar(base["per_project"], x="project", y="n_samples", text="n_samples",
                     labels={"n_samples": "samples"})
        fig.update_layout(height=280, margin=dict(t=10, b=10), showlegend=False)
        st.plotly_chart(fig)
    with c2:
        st.markdown("**Subjects by response**")
        st.dataframe(base["by_response"], width="stretch", hide_index=True)
        fig = px.bar(base["by_response"], x="response", y="n_subjects", text="n_subjects", color="response",
                     color_discrete_map=GROUP_COLOURS, labels={"n_subjects": "subjects"})
        fig.update_layout(height=280, margin=dict(t=10, b=10), showlegend=False)
        st.plotly_chart(fig)
    with c3:
        st.markdown("**Subjects by sex**")
        st.dataframe(base["by_sex"], width="stretch", hide_index=True)
        fig = px.bar(base["by_sex"], x="sex", y="n_subjects", text="n_subjects", color="sex",
                     labels={"n_subjects": "subjects"})
        fig.update_layout(height=280, margin=dict(t=10, b=10), showlegend=False)
        st.plotly_chart(fig)

    with st.expander(f"Show the {len(cohort)} baseline samples"):
        st.dataframe(cohort, width="stretch", hide_index=True, height=350)
        st.download_button("Download baseline samples (CSV)", cohort.to_csv(index=False).encode(),
                           "baseline_melanoma_miraclib_pbmc.csv", "text/csv")

    st.markdown("---")
    st.markdown("#### Average B-cell count - melanoma males, responders, time 0")
    st.markdown("All sample types and all treatment types are included, as specified in the brief.")
    avg = base["avg_b_cells"]
    st.metric("Mean b_cell count", "n/a" if avg is None else f"{avg:,.2f}",
              f"n = {base['avg_b_cells_n']} samples", delta_color="off")

    with st.expander("Explore: average count for any subset"):
        e1, e2, e3, e4 = st.columns(4)
        with e1:
            population = st.selectbox("Population", POPULATIONS, index=0, key="exp_pop")
            condition_e = option_or_all("Condition", options["all_conditions"], "exp_cond")
        with e2:
            sex_e = option_or_all("Sex", ["M", "F"], "exp_sex")
            response_e = option_or_all("Response", ["yes", "no"], "exp_resp")
        with e3:
            time_e = option_or_all("Time from treatment start", options["time"], "exp_time")
            treatment_e = option_or_all("Treatment", options["all_treatments"], "exp_treat")
        with e4:
            sample_type_e = option_or_all("Sample type", options["sample_type"], "exp_type")
        avg_e, n_e = explore_average(population, condition_e, sex_e, response_e, time_e, treatment_e, sample_type_e)
        st.metric(f"Mean {population} count", "n/a" if avg_e is None else f"{avg_e:,.2f}",
                  f"n = {n_e} samples", delta_color="off")
