"""Part 3 - compare cell-population frequencies between responders and non-responders.

The primary analysis follows the brief literally: it uses the per-sample
relative frequencies from the Part 2 summary table for melanoma patients
treated with miraclib (PBMC samples only) and compares responders with
non-responders for each population.

Because every subject contributes several samples (days 0, 7 and 14), the
sample-level test treats repeated measures from one person as independent.
Two complementary views are therefore also provided:

* a *subject-level* comparison that first averages each subject's samples, and
* a *per-time-point* comparison that tests each collection day separately.

Statistics reported per population
----------------------------------
* Mann-Whitney U (two-sided) - primary, distribution-free test on percentages.
* Welch's t-test - parametric companion (unequal variances).
* Cliff's delta - effect size in [-1, 1]; > 0 means responders tend to be higher.
* Benjamini-Hochberg adjusted p-values across the five populations.
"""
from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless rendering for the pipeline
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps

from .db import POPULATIONS, query

ALPHA = 0.05
RESPONSE_LABELS = {"yes": "responder", "no": "non-responder"}
GROUP_COLOURS = {"responder": "#2a9d8f", "non-responder": "#e76f51"}

RESPONDER_SQL = """
SELECT
    f.sample,
    m.subject_id                 AS subject,
    m.project_id                 AS project,
    m.sex,
    m.age,
    m.time_from_treatment_start  AS time_from_treatment_start,
    m.response,
    f.population,
    f.count,
    f.total_count,
    f.percentage
FROM sample_frequencies AS f
JOIN sample_metadata    AS m ON m.sample_id = f.sample
WHERE m.condition   = ?
  AND m.treatment   = ?
  AND m.sample_type = ?
  AND m.response IN ('yes', 'no')
ORDER BY f.display_order, f.sample
"""


def responder_dataset(
    conn: sqlite3.Connection,
    condition: str = "melanoma",
    treatment: str = "miraclib",
    sample_type: str = "PBMC",
) -> pd.DataFrame:
    """Per-sample population frequencies for the cohort under comparison.

    Defaults to the cohort in the brief: melanoma patients on miraclib, PBMC
    samples only, with a recorded response (yes/no).
    """
    df = query(conn, RESPONDER_SQL, (condition, treatment, sample_type))
    df["group"] = df["response"].map(RESPONSE_LABELS)
    return df


def subject_level_dataset(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse repeated samples to one value per subject and population (mean %)."""
    keys = ["subject", "project", "sex", "response", "group", "population"]
    return df.groupby(keys, as_index=False, sort=False)["percentage"].mean()


def bh_adjust(pvalues) -> np.ndarray:
    """Benjamini-Hochberg false-discovery-rate adjustment."""
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adjusted, 0.0, 1.0)
    return out


def cliffs_delta_from_u(u_stat: float, n1: int, n2: int) -> float:
    """Cliff's delta derived from the Mann-Whitney U statistic of group 1."""
    return 2.0 * u_stat / (n1 * n2) - 1.0


def _ordered_populations(df: pd.DataFrame) -> list[str]:
    present = list(dict.fromkeys(df["population"]))
    return [p for p in POPULATIONS if p in present] + [p for p in present if p not in POPULATIONS]


def compare_groups(df: pd.DataFrame, value_col: str = "percentage", alpha: float = ALPHA) -> pd.DataFrame:
    """Test responders vs non-responders for every population.

    ``df`` must contain ``population``, ``response`` ('yes'/'no') and ``value_col``.
    Returns one row per population with descriptive statistics, test results,
    effect sizes and BH-adjusted p-values.
    """
    rows = []
    for population in _ordered_populations(df):
        sub = df[df["population"] == population]
        resp = sub.loc[sub["response"] == "yes", value_col].to_numpy(dtype=float)
        non = sub.loc[sub["response"] == "no", value_col].to_numpy(dtype=float)
        row = {
            "population": population,
            "n_responders": int(resp.size),
            "n_non_responders": int(non.size),
            "mean_responders": float(np.mean(resp)) if resp.size else np.nan,
            "mean_non_responders": float(np.mean(non)) if non.size else np.nan,
            "median_responders": float(np.median(resp)) if resp.size else np.nan,
            "median_non_responders": float(np.median(non)) if non.size else np.nan,
        }
        row["median_difference"] = row["median_responders"] - row["median_non_responders"]
        if resp.size >= 2 and non.size >= 2:
            mwu = sps.mannwhitneyu(resp, non, alternative="two-sided")
            welch = sps.ttest_ind(resp, non, equal_var=False)
            row.update(
                mannwhitney_u=float(mwu.statistic),
                mannwhitney_p=float(mwu.pvalue),
                welch_t=float(welch.statistic),
                welch_p=float(welch.pvalue),
                cliffs_delta=cliffs_delta_from_u(float(mwu.statistic), resp.size, non.size),
            )
        else:
            row.update(mannwhitney_u=np.nan, mannwhitney_p=np.nan, welch_t=np.nan, welch_p=np.nan, cliffs_delta=np.nan)
        rows.append(row)

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["mannwhitney_p_adj"] = bh_adjust(out["mannwhitney_p"].fillna(1.0))
    out["welch_p_adj"] = bh_adjust(out["welch_p"].fillna(1.0))
    out["significant"] = out["mannwhitney_p_adj"] < alpha
    out["direction"] = np.where(
        out["significant"],
        np.where(out["cliffs_delta"] > 0, "higher in responders", "lower in responders"),
        "no significant difference",
    )
    return out


def compare_by_timepoint(df: pd.DataFrame, alpha: float = ALPHA) -> pd.DataFrame:
    """Run ``compare_groups`` separately for every time point."""
    frames = []
    for time_point, sub in df.groupby("time_from_treatment_start", sort=True):
        res = compare_groups(sub, alpha=alpha)
        res.insert(0, "time_from_treatment_start", time_point)
        frames.append(res)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def significant_populations(stats_df: pd.DataFrame) -> list[str]:
    return stats_df.loc[stats_df["significant"], "population"].tolist()


def _format_p(p: float) -> str:
    if np.isnan(p):
        return "n/a"
    return f"{p:.2e}" if p < 1e-3 else f"{p:.3f}"


def plot_boxplots(
    df: pd.DataFrame,
    path: Path | str,
    stats_df: pd.DataFrame | None = None,
    title: str = "Melanoma / miraclib / PBMC: responders vs non-responders",
) -> Path:
    """Save one boxplot per population (responders vs non-responders) as a PNG."""
    populations = _ordered_populations(df)
    fig, axes = plt.subplots(1, len(populations), figsize=(3.6 * len(populations), 5.2))
    axes = np.atleast_1d(axes)
    rng = np.random.default_rng(0)
    stats_lookup = stats_df.set_index("population") if stats_df is not None else None

    for ax, population in zip(axes, populations):
        sub = df[df["population"] == population]
        groups = ["responder", "non-responder"]
        data = [sub.loc[sub["group"] == g, "percentage"].to_numpy() for g in groups]
        box = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False)
        for patch, g in zip(box["boxes"], groups):
            patch.set_facecolor(GROUP_COLOURS[g])
            patch.set_alpha(0.55)
        for median in box["medians"]:
            median.set_color("black")
        for i, (values, g) in enumerate(zip(data, groups), start=1):
            jitter = rng.uniform(-0.16, 0.16, size=values.size)
            ax.scatter(np.full(values.size, i) + jitter, values, s=7, alpha=0.35,
                       color=GROUP_COLOURS[g], linewidths=0)
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"{g}\n(n={len(v)})" for g, v in zip(groups, data)], fontsize=9)
        ax.set_title(population, fontsize=11, fontweight="bold")
        ax.grid(axis="y", alpha=0.3)
        if stats_lookup is not None and population in stats_lookup.index:
            s = stats_lookup.loc[population]
            star = " *" if bool(s["significant"]) else ""
            ax.text(0.5, 0.97,
                    f"MWU p={_format_p(s['mannwhitney_p'])}, q={_format_p(s['mannwhitney_p_adj'])}{star}\n"
                    f"Cliff's δ={s['cliffs_delta']:+.2f}",
                    transform=ax.transAxes, ha="center", va="top", fontsize=8.5,
                    bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.85, linewidth=0))
    axes[0].set_ylabel("Relative frequency (% of total cells)")
    fig.suptitle(f"{title}\n(* = significant after Benjamini-Hochberg correction, q < {ALPHA})", fontsize=12)
    fig.tight_layout()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_single_boxplot(df: pd.DataFrame, population: str, path: Path | str,
                        stats_df: pd.DataFrame | None = None) -> Path:
    """Boxplot for one population only (one file per population)."""
    sub = df[df["population"] == population]
    single_stats = stats_df[stats_df["population"] == population] if stats_df is not None else None
    return plot_boxplots(sub, path, single_stats, title=f"{population}: responders vs non-responders")
