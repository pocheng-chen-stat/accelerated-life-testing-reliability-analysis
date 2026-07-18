"""
R-style visualization functions for the accelerated life testing project.

These functions intentionally reproduce the plotting logic of the original R
course project: transformed plotting coordinates with readable original-scale
axes on the bottom/left and transformed scales on the top/right.

The goal of this module is reproducibility of the original analysis figures,
not modern chart styling.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

from src.inverse_power_law import (
    failure_probability_confidence_band,
    predict_failure_probability,
)


R_COLORS = {
    100.3: "purple",
    122.4: "blue",
    157.1: "lime",
    219.0: "brown",
    361.4: "red",
    50.0: "black",
}


def _savefig(path: str | Path, fig=None, *, dpi: int = 150) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fig is None:
        fig = plt.gcf()
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def _fmt(value: float) -> str:
    if np.isclose(value, 0):
        return "0"
    if abs(value) >= 1_000_000:
        return f"{value:.0e}"
    if abs(value) >= 10_000:
        return f"{int(round(value)):d}"
    if abs(value - round(value)) < 1e-10:
        return f"{int(round(value))}"
    return f"{value:g}"


def _valid_positive(values: Iterable[float]) -> list[float]:
    return [float(v) for v in values if float(v) > 0]


def _set_four_sided_log_axes(
    ax,
    *,
    bottom_values: Iterable[float],
    left_values: Iterable[float],
    top_values: Iterable[float] | None,
    right_values: Iterable[float] | None,
    bottom_label: str,
    left_label: str,
    top_label: str,
    right_label: str,
) -> None:
    """Axes for plots whose coordinates are x=log(original x), y=log(original y)."""
    bottom_values = _valid_positive(bottom_values)
    left_values = _valid_positive(left_values)

    ax.set_xticks(np.log(bottom_values))
    ax.set_xticklabels([_fmt(v) for v in bottom_values])
    ax.set_yticks(np.log(left_values))
    ax.set_yticklabels([_fmt(v) for v in left_values])
    ax.set_xlabel(bottom_label)
    ax.set_ylabel(left_label)

    top = ax.secondary_xaxis("top")
    if top_values is not None:
        top_values = list(top_values)
        top.set_xticks(top_values)
        top.set_xticklabels([_fmt(v) for v in top_values])
    top.set_xlabel(top_label)

    right = ax.secondary_yaxis("right")
    if right_values is not None:
        right_values = list(right_values)
        right.set_yticks(right_values)
        right.set_yticklabels([_fmt(v) for v in right_values])
    right.set_ylabel(right_label)


def _set_probability_axes(
    ax,
    *,
    time_values: Iterable[float],
    probability_values: Iterable[float],
    top_log_time_values: Iterable[float] | None,
    right_quantile_values: Iterable[float] | None,
    bottom_label: str = "Time",
    left_label: str = "Proportion Failing",
    top_label: str = "log Time",
    right_label: str = "Standard quantile",
) -> None:
    """Axes for probability plots whose coordinates are x=log(time), y=qnorm(F)."""
    time_values = _valid_positive(time_values)
    probability_values = [float(p) for p in probability_values]

    ax.set_xticks(np.log(time_values))
    ax.set_xticklabels([_fmt(v) for v in time_values])
    ax.set_yticks(stats.norm.ppf(probability_values))
    ax.set_yticklabels([f"{p:g}" for p in probability_values])
    ax.set_xlabel(bottom_label)
    ax.set_ylabel(left_label)

    top = ax.secondary_xaxis("top")
    if top_log_time_values is not None:
        top_log_time_values = list(top_log_time_values)
        top.set_xticks(top_log_time_values)
        top.set_xticklabels([_fmt(v) for v in top_log_time_values])
    top.set_xlabel(top_label)

    right = ax.secondary_yaxis("right")
    if right_quantile_values is not None:
        right_quantile_values = list(right_quantile_values)
        right.set_yticks(right_quantile_values)
        right.set_yticklabels([_fmt(v) for v in right_quantile_values])
    right.set_ylabel(right_label)


def _mu_map(parameter_table: pd.DataFrame) -> dict[float, float]:
    rows = parameter_table[parameter_table["parameter"] == "mu"]
    return {float(v): float(m) for v, m in zip(rows["voltage_kv_mm"], rows["mle"])}


def _sigma_map(parameter_table: pd.DataFrame) -> dict[float, float]:
    rows = parameter_table[parameter_table["parameter"] == "sigma"]
    return {float(v): float(m) for v, m in zip(rows["voltage_kv_mm"], rows["mle"])}


def _equal_sigma(parameter_table: pd.DataFrame) -> float:
    return float(parameter_table.loc[parameter_table["parameter"] == "sigma", "mle"].iloc[0])


def _empirical_probability_data(df: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for voltage, group in df.groupby("voltage_kv_mm", sort=True):
        g = group.sort_values("failure_time").copy()
        n = len(g)
        g["p_failing"] = (np.arange(1, n + 1) - 0.5) / n
        frames.append(g)
    return pd.concat(frames, ignore_index=True)


# -----------------------------------------------------------------------------
# Original R-style figures
# -----------------------------------------------------------------------------


def plot_01_raw_scatter(df: pd.DataFrame, output_path: str | Path) -> None:
    """R block: plot(data_original2$V, data_original2$time)."""
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.scatter(df["voltage_kv_mm"], df["failure_time"], facecolors="none", edgecolors="black", s=22)
    ax.set_xlabel("kV / mm")
    ax.set_ylabel("Time")
    _savefig(output_path, fig)


def plot_02_log_scatter_with_regression(df: pd.DataFrame, output_path: str | Path) -> None:
    """R transformed scatter with included/excluded 361.4 kV/mm regression lines."""
    x = np.log(df["voltage_kv_mm"].to_numpy(dtype=float))
    y = np.log(df["failure_time"].to_numpy(dtype=float))

    all_beta = np.linalg.lstsq(np.column_stack([np.ones(len(x)), x]), y, rcond=None)[0]
    mask = df["voltage_kv_mm"].to_numpy(dtype=float) != 361.4
    excl_beta = np.linalg.lstsq(np.column_stack([np.ones(mask.sum()), x[mask]]), y[mask], rcond=None)[0]

    fig, ax = plt.subplots(figsize=(9.3, 6.5))
    ax.scatter(x, y, facecolors="none", edgecolors="black", s=20)
    x_grid = np.linspace(np.log(100), np.log(370), 250)
    ax.plot(x_grid, all_beta[0] + all_beta[1] * x_grid, color="red", linewidth=0.9, label="included 361.4kV/mm")
    ax.plot(x_grid, excl_beta[0] + excl_beta[1] * x_grid, color="blue", linewidth=0.9, label="excluded 361.4kV/mm")

    ax.set_xlim(np.log(100), np.log(370))
    ax.set_ylim(np.log(0.09), np.log(7200))
    _set_four_sided_log_axes(
        ax,
        bottom_values=[100.3, 122.4, 157.1, 219, 361.4],
        left_values=[0.1, 0.5, 1, 5, 10, 30, 100, 500, 1000, 3000, 6000, 7000],
        top_values=[4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8],
        right_values=[-2, 0, 2, 4, 6, 8],
        bottom_label="kV/mm",
        left_label="time",
        top_label="log kV/mm",
        right_label="log time",
    )
    ax.legend(loc="upper right", frameon=True, fancybox=False, edgecolor="black", fontsize=8)
    _savefig(output_path, fig)


def _draw_probability_points(ax, df: pd.DataFrame, *, include_legend: bool = True) -> None:
    empirical = _empirical_probability_data(df)
    for voltage, group in empirical.groupby("voltage_kv_mm", sort=True):
        color = R_COLORS.get(float(voltage), "black")
        label = f"{voltage:g} kV / mm" if include_legend else None
        ax.scatter(
            np.log(group["failure_time"].to_numpy(dtype=float)),
            stats.norm.ppf(group["p_failing"].to_numpy(dtype=float)),
            color=color,
            s=18,
            label=label,
        )


def plot_03_probability_plot_separate_sigma(
    df: pd.DataFrame,
    separate_params: pd.DataFrame,
    output_path: str | Path,
) -> None:
    """Multiple probability plot with sigma estimated separately for each stress level."""
    empirical = _empirical_probability_data(df)
    mu = _mu_map(separate_params)
    sigma = _sigma_map(separate_params)

    fig, ax = plt.subplots(figsize=(9.2, 6.2))
    _draw_probability_points(ax, empirical)
    xlim = (np.log(empirical["failure_time"].min()), np.log(empirical["failure_time"].max()))
    x_line = np.linspace(*xlim, 300)
    for voltage in sorted(empirical["voltage_kv_mm"].unique()):
        color = R_COLORS.get(float(voltage), "black")
        ax.plot(x_line, (x_line - mu[float(voltage)]) / sigma[float(voltage)], color=color, linewidth=0.8)

    ax.set_xlim(*xlim)
    ax.set_ylim(stats.norm.ppf(0.04), stats.norm.ppf(0.98))
    _set_probability_axes(
        ax,
        time_values=[0.1, 0.5, 1, 3, 10, 30, 100, 500, 1000, 3000, 6000, 7000],
        probability_values=[0.05, 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99],
        top_log_time_values=None,
        right_quantile_values=None,
    )
    ax.legend(loc="upper left", frameon=True, fancybox=False, edgecolor="black", fontsize=8)
    _savefig(output_path, fig)


def plot_04_probability_plot_equal_sigma(
    df: pd.DataFrame,
    equal_params: pd.DataFrame,
    output_path: str | Path,
) -> None:
    """Multiple probability plot with a common sigma across stress levels."""
    empirical = _empirical_probability_data(df)
    mu = _mu_map(equal_params)
    sigma = _equal_sigma(equal_params)

    fig, ax = plt.subplots(figsize=(9.2, 6.2))
    _draw_probability_points(ax, empirical)
    xlim = (np.log(empirical["failure_time"].min()), np.log(empirical["failure_time"].max()))
    x_line = np.linspace(*xlim, 300)
    for voltage in sorted(empirical["voltage_kv_mm"].unique()):
        color = R_COLORS.get(float(voltage), "black")
        ax.plot(x_line, (x_line - mu[float(voltage)]) / sigma, color=color, linewidth=0.8)

    ax.set_xlim(*xlim)
    ax.set_ylim(stats.norm.ppf(0.04), stats.norm.ppf(0.98))
    _set_probability_axes(
        ax,
        time_values=[0.1, 0.5, 1, 3, 10, 30, 100, 500, 1000, 3000, 6000, 7000],
        probability_values=[0.05, 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99],
        top_log_time_values=None,
        right_quantile_values=None,
    )
    ax.legend(loc="upper left", frameon=True, fancybox=False, edgecolor="black", fontsize=8)
    _savefig(output_path, fig)


def plot_05_probability_plot_overlay_sigma_assumptions(
    df: pd.DataFrame,
    separate_params: pd.DataFrame,
    equal_params: pd.DataFrame,
    output_path: str | Path,
) -> None:
    """Faithful to the R script when separate-sigma and equal-sigma lines are drawn on one canvas."""
    empirical = _empirical_probability_data(df)
    sep_mu, sep_sigma = _mu_map(separate_params), _sigma_map(separate_params)
    eq_mu, eq_sigma = _mu_map(equal_params), _equal_sigma(equal_params)

    fig, ax = plt.subplots(figsize=(9.2, 6.2))
    _draw_probability_points(ax, empirical)
    xlim = (np.log(empirical["failure_time"].min()), np.log(empirical["failure_time"].max()))
    x_line = np.linspace(*xlim, 300)
    for voltage in sorted(empirical["voltage_kv_mm"].unique()):
        color = R_COLORS.get(float(voltage), "black")
        ax.plot(x_line, (x_line - sep_mu[float(voltage)]) / sep_sigma[float(voltage)], color=color, linewidth=0.6, alpha=0.45)
        ax.plot(x_line, (x_line - eq_mu[float(voltage)]) / eq_sigma, color=color, linewidth=0.9)

    ax.set_xlim(*xlim)
    ax.set_ylim(stats.norm.ppf(0.04), stats.norm.ppf(0.98))
    _set_probability_axes(
        ax,
        time_values=[0.1, 0.5, 1, 3, 10, 30, 100, 500, 1000, 3000, 6000, 7000],
        probability_values=[0.05, 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99],
        top_log_time_values=None,
        right_quantile_values=None,
    )
    ax.legend(loc="upper left", frameon=True, fancybox=False, edgecolor="black", fontsize=8)
    _savefig(output_path, fig)


def plot_06_probability_plot_inverse_power_with_50kv_ci(
    df: pd.DataFrame,
    inverse_summary: dict,
    output_path: str | Path,
    *,
    prediction_voltage: float = 50.0,
    marker_time: float = 10000.0,
    time_grid_start: float = 1_000.0,
    time_grid_end: float = 10_000_000.0,
    x_min: float = 0.001,
    x_max: float = 150_000.0,
) -> None:
    """R probability plot: inverse power law lines + 50 kV/mm extrapolation + confidence band."""
    beta0 = inverse_summary["beta0"]
    beta1 = inverse_summary["beta1"]
    sigma = inverse_summary["sigma"]

    empirical = _empirical_probability_data(df)
    fig, ax = plt.subplots(figsize=(9.5, 6.5))
    _draw_probability_points(ax, empirical)

    # x_line = np.linspace(np.log(x_min), np.log(x_max), 400)
    # for voltage in sorted(empirical["voltage_kv_mm"].unique()):
    #     mu = beta0 + beta1 * np.log(float(voltage))
    #     ax.plot(x_line, (x_line - mu) / sigma, color=R_COLORS.get(float(voltage), "black"), linewidth=0.8)

    # mu_pred = beta0 + beta1 * np.log(prediction_voltage)
    # ax.plot(x_line, (x_line - mu_pred) / sigma, color="black", linewidth=1.1, label=f"{prediction_voltage:g} kV / mm")
    y_min = stats.norm.ppf(0.001)
    y_max = 3
    z_line = np.linspace(y_min, y_max, 400)

    for voltage in sorted(empirical["voltage_kv_mm"].unique()):
        mu = beta0 + beta1 * np.log(float(voltage))
        ax.plot(
            mu + sigma * z_line,
            z_line,
            color=R_COLORS.get(float(voltage), "black"),
            linewidth=0.8,
        )

    mu_pred = beta0 + beta1 * np.log(prediction_voltage)
    ax.plot(
        mu_pred + sigma * z_line,
        z_line,
        color="black",
        linewidth=1.1,
        label=f"{prediction_voltage:g} kV / mm",
    )

    # Use a geometric grid because the x-axis is log(time).  A fixed increment
    # in raw time produces visibly sparse points in the lower-time region and
    # unnecessarily dense points in the upper-time region.
    time_grid = np.geomspace(time_grid_start, time_grid_end, 3000, dtype=float)
    band = failure_probability_confidence_band(inverse_summary, prediction_voltage, time_grid)

    # Continuous lines communicate a confidence band more clearly than a cloud
    # of tiny points and remain smooth across the full probability scale.
    ax.plot(
        np.log(band["time"]),
        band["lower_95_quantile"],
        color="black",
        linewidth=1.25,
        alpha=0.9,
    )
    ax.plot(
        np.log(band["time"]),
        band["upper_95_quantile"],
        color="black",
        linewidth=1.25,
        alpha=0.9,
    )

    f_marker = predict_failure_probability(inverse_summary, prediction_voltage, marker_time)
    q_marker = stats.norm.ppf(f_marker)
    ax.scatter([np.log(marker_time)], [q_marker], marker="x", color="black", s=65, zorder=5)
    ax.axvline(np.log(marker_time), linestyle="--", color="gray", linewidth=0.9)
    ax.axhline(q_marker, linestyle="--", color="gray", linewidth=0.9)

    ax.set_xlim(np.log(x_min), np.log(x_max))
    ax.set_ylim(stats.norm.ppf(0.001), 3)
    _set_probability_axes(
        ax,
        time_values=[0.1, 0.5, 3, 10, 30, 100, 500, 1500, 3000, 7000, 20000, 70000, 150000],
        probability_values=[0.002, 0.05, 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99],
        top_log_time_values=None,
        right_quantile_values=None,
    )
    handles, labels = ax.get_legend_handles_labels()
    # The data points already include observed stress legends; add 50 kV/mm.
    ax.legend(loc="upper left", frameon=True, fancybox=False, edgecolor="black", fontsize=8)
    _savefig(output_path, fig)


def plot_07_stress_life_quantiles(
    df: pd.DataFrame,
    inverse_summary: dict,
    output_path: str | Path,
    *,
    prediction_voltage: float = 50.0,
) -> None:
    """R stress-life figure: quantile lines and vertical lognormal density traces."""
    beta0 = inverse_summary["beta0"]
    beta1 = inverse_summary["beta1"]
    sigma = inverse_summary["sigma"]

    fig, ax = plt.subplots(figsize=(9.5, 6.8))
    ax.scatter(
        np.log(df["voltage_kv_mm"].to_numpy(dtype=float)),
        np.log(df["failure_time"].to_numpy(dtype=float)),
        facecolors="none",
        edgecolors="black",
        s=20,
        linewidths=0.7,
    )

    x_grid = np.linspace(np.log(40), np.log(370), 400)
    for prob, label_y in [(0.9, 2.2), (0.5, 0.8), (0.1, -0.9)]:
        y_grid = beta0 + beta1 * x_grid + stats.norm.ppf(prob) * sigma
        ax.plot(x_grid, y_grid, color="black", linewidth=0.7)
        ax.text(5.96, label_y, f"{int(prob * 100)} %", fontsize=8)

    for voltage in [prediction_voltage] + sorted(df["voltage_kv_mm"].unique().tolist()):
        mu = beta0 + beta1 * np.log(float(voltage))
        y_seq = np.linspace(mu + stats.norm.ppf(0.001) * sigma, mu + stats.norm.ppf(0.999) * sigma, 300)
        density = stats.norm.pdf(y_seq, loc=mu, scale=sigma)
        ax.plot(np.repeat(np.log(voltage), len(y_seq)), y_seq, color="black", linewidth=0.55)
        ax.plot(np.log(voltage) - density * 0.5, y_seq, color="black", linewidth=0.55)

    ax.set_xlim(np.log(40), np.log(370))
    ax.set_ylim(np.log(0.09), np.log(1e7))
    bottom_ticks = [prediction_voltage] + sorted(df["voltage_kv_mm"].unique().tolist())
    _set_four_sided_log_axes(
        ax,
        bottom_values=bottom_ticks,
        left_values=[0.1, 1, 5, 30, 100, 500, 3000, 20000, 1e5, 1e6, 1e7],
        top_values=[4.0, 4.5, 5.0, 5.5, 6.0],
        right_values=[0, 5, 10, 15],
        bottom_label="kV/mm",
        left_label="time",
        top_label="log kV/mm",
        right_label="log time",
    )
    _savefig(output_path, fig)



def plot_12_stress_life_quantile_comparison(
    df: pd.DataFrame,
    included_summary: dict,
    excluded_summary: dict,
    output_path: str | Path,
    *,
    prediction_voltage: float = 50.0,
) -> None:
    """
    Compare stress-life quantile curves from models fitted with and without
    the 361.4 kV/mm stress group.

    Each model uses one color for all three quantiles. The 10%, 50%, and 90%
    labels are written directly on the corresponding curves so the legend only
    needs to distinguish the two model fits.
    """
    included_color = "#D97706"  # orange: model including 361.4 kV/mm
    excluded_color = "#1F4E79"  # blue: model excluding 361.4 kV/mm

    fig, ax = plt.subplots(figsize=(10.2, 6.8))

    high_stress = np.isclose(df["voltage_kv_mm"].to_numpy(dtype=float), 361.4)
    regular = ~high_stress

    # Common observations are kept neutral. The 361.4 kV/mm observations are
    # highlighted with the same color as the model that includes them.
    ax.scatter(
        np.log(df.loc[regular, "voltage_kv_mm"].to_numpy(dtype=float)),
        np.log(df.loc[regular, "failure_time"].to_numpy(dtype=float)),
        facecolors="none",
        edgecolors="black",
        s=25,
        linewidths=0.8,
        zorder=3,
        label="Observed data",
    )
    ax.scatter(
        np.log(df.loc[high_stress, "voltage_kv_mm"].to_numpy(dtype=float)),
        np.log(df.loc[high_stress, "failure_time"].to_numpy(dtype=float)),
        facecolors="none",
        edgecolors=included_color,
        s=30,
        linewidths=1.1,
        zorder=4,
        label="361.4 kV/mm observations",
    )

    x_grid = np.linspace(np.log(40), np.log(370), 500)
    quantiles = [0.9, 0.5, 0.1]

    def draw_model(
        summary: dict,
        *,
        color: str,
        legend_label: str,
        label_voltage: float,
    ) -> None:
        beta0 = float(summary["beta0"])
        beta1 = float(summary["beta1"])
        sigma = float(summary["sigma"])

        for index, probability in enumerate(quantiles):
            y_grid = beta0 + beta1 * x_grid + stats.norm.ppf(probability) * sigma
            ax.plot(
                x_grid,
                y_grid,
                color=color,
                linewidth=1.55,
                label=legend_label if index == 0 else None,
                zorder=2,
            )

            x_text = np.log(label_voltage)
            y_text = (
                beta0
                + beta1 * x_text
                + stats.norm.ppf(probability) * sigma
            )
            ax.text(
                x_text,
                y_text,
                f"{int(probability * 100)}%",
                color=color,
                fontsize=9,
                fontweight="bold",
                ha="left",
                va="center",
                bbox={
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.78,
                    "pad": 0.35,
                },
                zorder=5,
            )

    # Put the labels at different x positions to keep the two model families
    # visually separate while preserving the actual fitted curves.
    draw_model(
        included_summary,
        color=included_color,
        legend_label="Including 361.4 kV/mm",
        label_voltage=315.0,
    )
    draw_model(
        excluded_summary,
        color=excluded_color,
        legend_label="Excluding 361.4 kV/mm",
        label_voltage=235.0,
    )

    ax.set_xlim(np.log(40), np.log(370))
    ax.set_ylim(np.log(0.09), np.log(1e7))
    bottom_ticks = [prediction_voltage, 100.3, 122.4, 157.1, 219.0, 361.4]
    _set_four_sided_log_axes(
        ax,
        bottom_values=bottom_ticks,
        left_values=[0.1, 1, 5, 30, 100, 500, 3000, 20000, 1e5, 1e6, 1e7],
        top_values=[4.0, 4.5, 5.0, 5.5, 6.0],
        right_values=[0, 5, 10, 15],
        bottom_label="kV/mm",
        left_label="time",
        top_label="log kV/mm",
        right_label="log time",
    )

    # Reorder the legend so the two fitted models appear first.
    handles, labels = ax.get_legend_handles_labels()
    preferred = [
        "Including 361.4 kV/mm",
        "Excluding 361.4 kV/mm",
        "Observed data",
        "361.4 kV/mm observations",
    ]
    lookup = dict(zip(labels, handles))
    ax.legend(
        [lookup[label] for label in preferred if label in lookup],
        [label for label in preferred if label in lookup],
        loc="upper right",
        frameon=True,
        fancybox=False,
        edgecolor="black",
        fontsize=8,
    )
    _savefig(output_path, fig, dpi=200)

def inverse_power_residuals(
    df: pd.DataFrame,
    inverse_summary: dict,
) -> pd.DataFrame:
    """
    Standardized residuals for the inverse power law model.

    This uses sigma estimated from the inverse-power-law regression itself:
        z_i = (log(t_i) - beta0 - beta1 log(v_i)) / sigma_IPL

    The visual design still follows the original R figure: residuals are plotted
    on transformed coordinates, while the bottom/left axes show readable original
    scales and the top/right axes show transformed scales.
    """
    beta0 = inverse_summary["beta0"]
    beta1 = inverse_summary["beta1"]
    sigma_ipl = inverse_summary["sigma"]

    v = df["voltage_kv_mm"].to_numpy(dtype=float)
    t = df["failure_time"].to_numpy(dtype=float)
    fitted_mu = beta0 + beta1 * np.log(v)
    log_residual = (np.log(t) - fitted_mu) / sigma_ipl

    return pd.DataFrame(
        {
            "voltage_kv_mm": v,
            "failure_time": t,
            "fitted_value": np.exp(fitted_mu),
            "log_standardized_residual": log_residual,
            "exp_standardized_residual": np.exp(log_residual),
            "is_outlier_abs_log_residual_gt_2": np.abs(log_residual) > 2,
            "sigma_used": sigma_ipl,
        }
    )


def equal_sigma_residuals(df: pd.DataFrame, equal_params: pd.DataFrame) -> pd.DataFrame:
    mu = _mu_map(equal_params)
    sigma = _equal_sigma(equal_params)
    rows = []
    for _, row in df.iterrows():
        z = (np.log(row["failure_time"]) - mu[float(row["voltage_kv_mm"])]) / sigma
        rows.append(
            {
                "voltage_kv_mm": float(row["voltage_kv_mm"]),
                "failure_time": float(row["failure_time"]),
                "fitted_value": float(np.exp(mu[float(row["voltage_kv_mm"])])),
                "log_standardized_residual": float(z),
                "exp_standardized_residual": float(np.exp(z)),
            }
        )
    return pd.DataFrame(rows)


def plot_08_standardized_residuals_inverse_power(
    df: pd.DataFrame,
    inverse_summary: dict,
    output_path: str | Path,
) -> pd.DataFrame:
    """R-style standardized residuals versus fitted values using sigma_IPL."""
    residuals = inverse_power_residuals(df, inverse_summary)

    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    ax.scatter(
        np.log(residuals["fitted_value"]),
        np.log(residuals["exp_standardized_residual"]),
        color="black",
        s=18,
    )
    out = residuals["is_outlier_abs_log_residual_gt_2"].to_numpy(dtype=bool)
    if out.any():
        ax.scatter(
            np.log(residuals.loc[out, "fitted_value"]),
            np.log(residuals.loc[out, "exp_standardized_residual"]),
            color="red",
            s=20,
            label="standardized residuals > 2",
        )
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlim(-0.2, 10)
    ax.set_ylim(np.log(0.08), np.log(22))

    _set_four_sided_log_axes(
        ax,
        bottom_values=[0.8, 10, 30, 80, 330, 795, 2630, 8000],
        left_values=[0.1, 0.2, 0.4, 0.6, 1, 1.6, 3, 4, 7, 10, 14, 19],
        top_values=[0, 2, 4, 6, 8, 10],
        right_values=[-2, -1, 0, 1, 2, 3],
        bottom_label="Time",
        left_label="exponential Standardized Residuals",
        top_label="Log Time",
        right_label="Standardized Residuals",
    )
    ax.legend(loc="upper left", frameon=True, fancybox=False, edgecolor="black", fontsize=8)
    _savefig(output_path, fig)
    return residuals


def plot_09_residual_probability_plot_equal_sigma(
    df: pd.DataFrame,
    equal_params: pd.DataFrame,
    output_path: str | Path,
) -> pd.DataFrame:
    """R-style Q-Q/probability plot of exponential standardized residuals under equal-sigma model."""
    residuals = equal_sigma_residuals(df, equal_params)
    n = len(residuals)
    p = (np.arange(1, n + 1) - 0.5) / n
    x = np.sort(np.log(residuals["exp_standardized_residual"].to_numpy(dtype=float)))
    y = stats.norm.ppf(p)

    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    ax.scatter(x, y, color="black", s=18)
    ax.plot([min(x.min(), y.min()), max(x.max(), y.max())], [min(x.min(), y.min()), max(x.max(), y.max())], color="gray", linewidth=0.8)

    ax.set_xlim(np.log(0.08), np.log(9))
    ax.set_ylim(stats.norm.ppf(0.015), stats.norm.ppf(0.995))
    _set_probability_axes(
        ax,
        time_values=[0.1, 0.2, 0.5, 1, 2, 3, 4, 5, 7, 8],
        probability_values=[0.002, 0.05, 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99],
        top_log_time_values=[-2, -1, 0, 1, 2],
        right_quantile_values=[-2, -1, 0, 1, 2],
        bottom_label="exponential Standardized Residuals",
        left_label="Probability",
        top_label="Standardized Residuals",
        right_label="Standard quantile",
    )
    _savefig(output_path, fig)
    return residuals


def plot_supplemental_aic_bar(model_summary: pd.DataFrame, output_path: str | Path) -> None:
    """Supplemental GitHub-friendly AIC bar chart; not in the original R script."""
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    ax.bar(model_summary["model"], model_summary["AIC"])
    ax.set_ylabel("AIC")
    ax.set_title("Model comparison by AIC")
    ax.tick_params(axis="x", rotation=18)
    _savefig(output_path, fig)
