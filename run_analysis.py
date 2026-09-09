#!/usr/bin/env python3
"""Parts 2-4 - run the complete analysis against ``cell_counts.db``.

Usage (after ``python load_data.py``)::

    python run_analysis.py

Every result is written to ``outputs/`` (CSV tables, PNG boxplots, a JSON
summary and a human-readable ``REPORT.md``) and a short summary is printed.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from cellcount import DB_PATH, OUTPUT_DIR, connect
from cellcount.stats import (
    ALPHA,
    compare_by_timepoint,
    compare_groups,
    plot_boxplots,
    plot_single_boxplot,
    responder_dataset,
    significant_populations,
    subject_level_dataset,
)
from cellcount.subsets import (
    BASELINE_FILTERS,
    baseline_samples,
    melanoma_male_responder_baseline_b_cells,
    samples_per_project,
    subjects_by_response,
    subjects_by_sex,
)
from cellcount.summary import summary_table

COHORT = {"condition": "melanoma", "treatment": "miraclib", "sample_type": "PBMC"}

STATS_COLUMNS = [
    "population", "n_responders", "n_non_responders",
    "mean_responders", "mean_non_responders", "median_responders", "median_non_responders",
    "median_difference", "cliffs_delta",
    "mannwhitney_u", "mannwhitney_p", "mannwhitney_p_adj",
    "welch_t", "welch_p", "welch_p_adj",
    "significant", "direction",
]


# --------------------------------------------------------------------------- helpers
def md_table(df: pd.DataFrame, float_fmt: str = "{:.3f}") -> str:
    """Render a DataFrame as a GitHub-flavoured markdown table (no extra deps)."""
    def fmt(value):
        if isinstance(value, (float, np.floating)):
            if np.isnan(value):
                return "n/a"
            return f"{value:.2e}" if 0 < abs(value) < 1e-3 else float_fmt.format(value)
        if isinstance(value, (bool, np.bool_)):
            return "yes" if value else "no"
        return str(value)

    header = "| " + " | ".join(df.columns) + " |"
    divider = "|" + "|".join(" --- " for _ in df.columns) + "|"
    body = ["| " + " | ".join(fmt(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join([header, divider, *body])


def fmt_p(p: float) -> str:
    return "n/a" if pd.isna(p) else (f"{p:.2e}" if p < 1e-3 else f"{p:.3f}")


def interpret_part3(sample_stats: pd.DataFrame, subject_stats: pd.DataFrame, alpha: float) -> list[str]:
    """Data-driven interpretation of the responder comparison."""
    lines = []
    sig = significant_populations(sample_stats)
    nominal = sample_stats.loc[sample_stats["mannwhitney_p"] < alpha, "population"].tolist()
    if sig:
        lines.append(
            f"After Benjamini-Hochberg correction (q < {alpha}), the following populations differ "
            f"significantly between responders and non-responders: {', '.join(sig)}."
        )
        for _, r in sample_stats[sample_stats["significant"]].iterrows():
            lines.append(
                f"- **{r['population']}** is {r['direction']} "
                f"(median {r['median_responders']:.2f}% vs {r['median_non_responders']:.2f}%, "
                f"Mann-Whitney p = {fmt_p(r['mannwhitney_p'])}, q = {fmt_p(r['mannwhitney_p_adj'])}, "
                f"Cliff's delta = {r['cliffs_delta']:+.3f})."
            )
    else:
        lines.append(
            f"No population differs significantly between responders and non-responders once the "
            f"five tests are corrected for multiple comparisons (all Benjamini-Hochberg q >= {alpha})."
        )
    if nominal:
        parts = []
        for _, r in sample_stats[sample_stats["mannwhitney_p"] < alpha].iterrows():
            direction = "higher" if r["cliffs_delta"] > 0 else "lower"
            parts.append(
                f"{r['population']} ({direction} in responders; median {r['median_responders']:.2f}% vs "
                f"{r['median_non_responders']:.2f}%, p = {fmt_p(r['mannwhitney_p'])}, "
                f"q = {fmt_p(r['mannwhitney_p_adj'])}, Cliff's delta = {r['cliffs_delta']:+.3f}, "
                f"Welch q = {fmt_p(r['welch_p_adj'])})"
            )
        lines.append(
            f"Nominally significant before correction (uncorrected p < {alpha}): " + "; ".join(parts) + "."
        )
    else:
        lines.append(f"No population is even nominally significant (uncorrected p < {alpha}).")

    # Agreement with the subject-level (pseudo-replication-free) analysis.
    subj_nominal = subject_stats.loc[subject_stats["mannwhitney_p"] < alpha, "population"].tolist()
    subj_sig = significant_populations(subject_stats)
    lines.append(
        "Averaging each subject's samples first (one value per subject, which removes the repeated-measures "
        f"dependence) gives the same picture: nominal p < {alpha} for "
        f"{', '.join(subj_nominal) if subj_nominal else 'no population'}; "
        f"significant after correction: {', '.join(subj_sig) if subj_sig else 'none'}."
    )
    strongest = sample_stats.loc[sample_stats["mannwhitney_p"].idxmin()]
    lines.append(
        f"The strongest candidate predictor of miraclib response is **{strongest['population']}** "
        f"({'higher' if strongest['cliffs_delta'] > 0 else 'lower'} in responders), but the effect is small "
        f"(Cliff's delta = {strongest['cliffs_delta']:+.3f}, i.e. a {abs(strongest['cliffs_delta'])*50 + 50:.0f}% "
        f"chance that a random responder exceeds a random non-responder) and should be confirmed in an "
        f"independent cohort before being used to predict response."
    )
    return lines


# --------------------------------------------------------------------------- parts
def run_part2(conn, out: Path) -> pd.DataFrame:
    df = summary_table(conn)
    path = out / "part2_summary_table.csv"
    df.to_csv(path, index=False, float_format="%.4f")
    print(f"[Part 2] Summary table: {len(df):,} rows ({df['sample'].nunique():,} samples x "
          f"{df['population'].nunique()} populations) -> {path.relative_to(out.parent)}")
    print(df.head(10).to_string(index=False))
    return df


def run_part3(conn, out: Path) -> dict:
    df = responder_dataset(conn, **COHORT)
    df.to_csv(out / "part3_responder_dataset.csv", index=False, float_format="%.4f")

    sample_stats = compare_groups(df)[STATS_COLUMNS]
    subject_stats = compare_groups(subject_level_dataset(df))[STATS_COLUMNS]
    time_stats = compare_by_timepoint(df)[["time_from_treatment_start", *STATS_COLUMNS]]
    sample_stats.to_csv(out / "part3_stats_sample_level.csv", index=False)
    subject_stats.to_csv(out / "part3_stats_subject_level.csv", index=False)
    time_stats.to_csv(out / "part3_stats_by_timepoint.csv", index=False)

    figure = plot_boxplots(df, out / "part3_boxplots.png", sample_stats)
    per_population = {}
    for population in sample_stats["population"]:
        per_population[population] = str(
            plot_single_boxplot(df, population, out / f"part3_boxplot_{population}.png", sample_stats).name
        )

    groups = df.groupby("group")
    sizes = {
        "responder_samples": int(groups["sample"].nunique().get("responder", 0)),
        "non_responder_samples": int(groups["sample"].nunique().get("non-responder", 0)),
        "responder_subjects": int(groups["subject"].nunique().get("responder", 0)),
        "non_responder_subjects": int(groups["subject"].nunique().get("non-responder", 0)),
    }
    interpretation = interpret_part3(sample_stats, subject_stats, ALPHA)

    print(f"\n[Part 3] Cohort: {COHORT} -> {sizes['responder_samples']} responder samples "
          f"({sizes['responder_subjects']} subjects) vs {sizes['non_responder_samples']} non-responder samples "
          f"({sizes['non_responder_subjects']} subjects)")
    show = sample_stats[["population", "median_responders", "median_non_responders", "cliffs_delta",
                         "mannwhitney_p", "mannwhitney_p_adj", "welch_p", "welch_p_adj", "significant"]]
    print(show.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print("Boxplots ->", figure.relative_to(out.parent))
    for line in interpretation:
        print(" ", line.replace("**", ""))

    return {
        "cohort": COHORT,
        "alpha": ALPHA,
        "sizes": sizes,
        "sample_level": sample_stats,
        "subject_level": subject_stats,
        "by_timepoint": time_stats,
        "significant_populations": significant_populations(sample_stats),
        "nominally_significant_populations": sample_stats.loc[sample_stats["mannwhitney_p"] < ALPHA, "population"].tolist(),
        "interpretation": interpretation,
        "figure": figure.name,
        "figures_per_population": per_population,
    }


def run_part4(conn, out: Path) -> dict:
    cohort = baseline_samples(conn)
    per_project = samples_per_project(conn)
    by_response = subjects_by_response(conn)
    by_sex = subjects_by_sex(conn)
    avg_b, n_b = melanoma_male_responder_baseline_b_cells(conn)

    cohort.to_csv(out / "part4_baseline_samples.csv", index=False)
    per_project.to_csv(out / "part4_samples_per_project.csv", index=False)
    by_response.to_csv(out / "part4_subjects_by_response.csv", index=False)
    by_sex.to_csv(out / "part4_subjects_by_sex.csv", index=False)
    avg_text = "n/a" if avg_b is None else f"{avg_b:.2f}"
    (out / "part4_avg_b_cells_melanoma_male_responders_baseline.txt").write_text(
        f"{avg_text}\n# average b_cell count, melanoma, male, response=yes, time_from_treatment_start=0, "
        f"all sample types and treatments (n={n_b} samples)\n", encoding="utf-8"
    )

    print(f"\n[Part 4] Baseline cohort {BASELINE_FILTERS}: {len(cohort)} samples from "
          f"{cohort['subject'].nunique()} subjects")
    print("Samples per project:\n" + per_project.to_string(index=False))
    print("Subjects by response:\n" + by_response.to_string(index=False))
    print("Subjects by sex:\n" + by_sex.to_string(index=False))
    print(f"Average B-cell count, melanoma males, responders, time 0 (all sample & treatment types): "
          f"{avg_text} (n={n_b} samples)")

    return {
        "filters": BASELINE_FILTERS,
        "n_samples": int(len(cohort)),
        "n_subjects": int(cohort["subject"].nunique()),
        "samples_per_project": per_project,
        "subjects_by_response": by_response,
        "subjects_by_sex": by_sex,
        "avg_b_cells_melanoma_male_responders_baseline": None if avg_b is None else round(avg_b, 2),
        "avg_b_cells_n_samples": n_b,
    }


# --------------------------------------------------------------------------- reporting
def write_results_json(out: Path, part2: pd.DataFrame, part3: dict, part4: dict) -> Path:
    def records(df: pd.DataFrame) -> list[dict]:
        return json.loads(df.to_json(orient="records"))

    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "database": DB_PATH.name,
        "part2": {
            "rows": int(len(part2)),
            "samples": int(part2["sample"].nunique()),
            "populations": part2["population"].unique().tolist(),
            "file": "part2_summary_table.csv",
        },
        "part3": {
            **{k: v for k, v in part3.items() if k not in ("sample_level", "subject_level", "by_timepoint")},
            "sample_level": records(part3["sample_level"]),
            "subject_level": records(part3["subject_level"]),
            "by_timepoint": records(part3["by_timepoint"]),
        },
        "part4": {
            **{k: v for k, v in part4.items()
               if k not in ("samples_per_project", "subjects_by_response", "subjects_by_sex")},
            "samples_per_project": records(part4["samples_per_project"]),
            "subjects_by_response": records(part4["subjects_by_response"]),
            "subjects_by_sex": records(part4["subjects_by_sex"]),
        },
    }
    path = out / "results.json"
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def write_report(out: Path, part2: pd.DataFrame, part3: dict, part4: dict) -> Path:
    s = part3["sizes"]
    stats_view = part3["sample_level"][[
        "population", "n_responders", "n_non_responders", "median_responders", "median_non_responders",
        "cliffs_delta", "mannwhitney_p", "mannwhitney_p_adj", "welch_p", "welch_p_adj", "significant",
    ]]
    subject_view = part3["subject_level"][[
        "population", "n_responders", "n_non_responders", "median_responders", "median_non_responders",
        "cliffs_delta", "mannwhitney_p", "mannwhitney_p_adj", "significant",
    ]]
    time_view = part3["by_timepoint"][[
        "time_from_treatment_start", "population", "n_responders", "n_non_responders",
        "cliffs_delta", "mannwhitney_p", "mannwhitney_p_adj", "significant",
    ]]
    avg_b = part4["avg_b_cells_melanoma_male_responders_baseline"]
    avg_text = "n/a" if avg_b is None else f"{avg_b:.2f}"

    md = f"""# Analysis report

Generated by `run_analysis.py` on {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')} from `{DB_PATH.name}`.

## Part 2 - Relative frequency of each cell population per sample

`part2_summary_table.csv` holds {len(part2):,} rows: one per (sample, population) for
{part2['sample'].nunique():,} samples. `total_count` is the sum of the five population counts of the
sample and `percentage = 100 * count / total_count`. First rows:

{md_table(part2.head(10), "{:.2f}")}

## Part 3 - Responders vs non-responders (melanoma, miraclib, PBMC)

Cohort: condition = melanoma, treatment = miraclib, sample type = PBMC, response recorded.
{s['responder_samples']} responder samples from {s['responder_subjects']} subjects vs
{s['non_responder_samples']} non-responder samples from {s['non_responder_subjects']} subjects
(each subject contributes samples at days 0, 7 and 14).

Per-population tests on the per-sample relative frequencies (primary analysis). `mannwhitney_p_adj`
and `welch_p_adj` are Benjamini-Hochberg adjusted across the five populations; `significant` means
adjusted Mann-Whitney p < {part3['alpha']}. Cliff's delta > 0 means responders tend to have higher values.

{md_table(stats_view, "{:.3f}")}

![Boxplots]({part3['figure']})

### Interpretation

{chr(10).join(part3['interpretation'])}

### Sensitivity analysis 1 - one value per subject (mean of the subject's samples)

{md_table(subject_view, "{:.3f}")}

### Sensitivity analysis 2 - each time point separately

{md_table(time_view, "{:.3f}")}

## Part 4 - Baseline melanoma PBMC samples from miraclib-treated patients

Filters: {part4['filters']} -> **{part4['n_samples']} samples** from **{part4['n_subjects']} subjects**
(`part4_baseline_samples.csv`).

Samples per project:

{md_table(part4['samples_per_project'])}

Subjects by response:

{md_table(part4['subjects_by_response'])}

Subjects by sex:

{md_table(part4['subjects_by_sex'])}

### Average B-cell count - melanoma males, responders, time 0 (all sample and treatment types)

**{avg_text}** cells (mean of the raw `b_cell` count over {part4['avg_b_cells_n_samples']} samples).
"""
    path = out / "REPORT.md"
    path.write_text(md, encoding="utf-8")
    return path


def main() -> int:
    if not DB_PATH.exists():
        print(f"ERROR: {DB_PATH.name} not found. Run `python load_data.py` first.", file=sys.stderr)
        return 1
    OUTPUT_DIR.mkdir(exist_ok=True)
    conn = connect(DB_PATH)
    try:
        part2 = run_part2(conn, OUTPUT_DIR)
        part3 = run_part3(conn, OUTPUT_DIR)
        part4 = run_part4(conn, OUTPUT_DIR)
    finally:
        conn.close()
    results = write_results_json(OUTPUT_DIR, part2, part3, part4)
    report = write_report(OUTPUT_DIR, part2, part3, part4)
    print(f"\nWrote {results.relative_to(OUTPUT_DIR.parent)} and {report.relative_to(OUTPUT_DIR.parent)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
