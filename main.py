from pathlib import Path

import pandas as pd

from src.data_loader import load_lifetime_data
from src.lognormal_models import fit_lognormal_equal_sigma, fit_lognormal_separate_sigma
from src.inverse_power_law import (
    fit_inverse_power_law,
    predict_failure_probability,
    predict_lifetime_quantile,
    failure_probability_confidence_interval,
    lifetime_quantile_confidence_interval,
)
from src.model_selection import calculate_aic, likelihood_ratio_test
from src.r_style_visualizations import (
    plot_01_raw_scatter,
    plot_02_log_scatter_with_regression,
    plot_03_probability_plot_separate_sigma,
    plot_04_probability_plot_equal_sigma,
    plot_05_probability_plot_overlay_sigma_assumptions,
    plot_06_probability_plot_inverse_power_with_50kv_ci,
    plot_07_stress_life_quantiles,
    plot_08_standardized_residuals_inverse_power,
    plot_09_residual_probability_plot_equal_sigma,
    plot_12_stress_life_quantile_comparison,
    plot_supplemental_aic_bar,
)


ROOT = Path(__file__).resolve().parent
DATA_PATH = ROOT / "data" / "insulation_lifetime_data.csv"
REPORTS_DIR = ROOT / "reports"
FIGURES_DIR = REPORTS_DIR / "figures"


def _prediction_table(summary: dict, *, label: str) -> pd.DataFrame:
    f_hat = predict_failure_probability(summary, voltage_kv_mm=50, time=10000)
    f_lower, f_upper = failure_probability_confidence_interval(
        summary, voltage_kv_mm=50, time=10000
    )
    b10 = predict_lifetime_quantile(summary, voltage_kv_mm=50, probability=0.1)
    b10_lower, b10_upper = lifetime_quantile_confidence_interval(
        summary, voltage_kv_mm=50, probability=0.1
    )
    b50 = predict_lifetime_quantile(summary, voltage_kv_mm=50, probability=0.5)
    b50_lower, b50_upper = lifetime_quantile_confidence_interval(
        summary, voltage_kv_mm=50, probability=0.5
    )

    return pd.DataFrame(
        [
            {
                "case": label,
                "voltage_kv_mm": 50,
                "F_10000": f_hat,
                "F_10000_lower_95": f_lower,
                "F_10000_upper_95": f_upper,
                "t_0.1": b10,
                "t_0.1_lower_95": b10_lower,
                "t_0.1_upper_95": b10_upper,
                "t_0.5": b50,
                "t_0.5_lower_95": b50_lower,
                "t_0.5_upper_95": b50_upper,
            }
        ]
    )


def main() -> None:
    df = load_lifetime_data(DATA_PATH)

    REPORTS_DIR.mkdir(exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # Lognormal models used in the original analysis.
    separate_params, separate_summary = fit_lognormal_separate_sigma(df)
    equal_params, equal_summary = fit_lognormal_equal_sigma(df)

    # Inverse power law using all stress levels. Use the covariance matrix for
    # the actual predictor log(V), matching the fitted model specification.
    inverse_params, inverse_summary = fit_inverse_power_law(
        df,
        label="Inverse power law",
        covariance_method="standard",
    )

    # Refit after excluding 361.4 kV/mm, matching the later part of the report.
    df_without_361 = df[df["voltage_kv_mm"] != 361.4].copy()
    inverse_params_without_361, inverse_summary_without_361 = fit_inverse_power_law(
        df_without_361,
        label="Inverse power law excluding 361.4 kV/mm",
        covariance_method="standard",
    )

    # Model-selection table.
    summaries = [separate_summary, equal_summary, inverse_summary]
    for item in summaries:
        item["AIC"] = calculate_aic(item["log_likelihood"], item["num_params"])
        item["minus_2_log_likelihood"] = -2 * item["log_likelihood"]

    model_summary = pd.DataFrame(
        [
            {
                "model": item["model"],
                "minus_2_log_likelihood": item["minus_2_log_likelihood"],
                "AIC": item["AIC"],
                "num_params": item["num_params"],
                "log_likelihood": item["log_likelihood"],
            }
            for item in summaries
        ]
    )

    lrt_equal_vs_separate = likelihood_ratio_test(
        ll_reduced=equal_summary["log_likelihood"],
        ll_full=separate_summary["log_likelihood"],
        df=separate_summary["num_params"] - equal_summary["num_params"],
    )

    lrt_inverse_vs_equal = likelihood_ratio_test(
        ll_reduced=inverse_summary["log_likelihood"],
        ll_full=equal_summary["log_likelihood"],
        df=equal_summary["num_params"] - inverse_summary["num_params"],
    )

    prediction_all = _prediction_table(inverse_summary, label="all_stress_levels")
    prediction_without_361 = _prediction_table(inverse_summary_without_361, label="excluding_361_4")

    # Tables.
    separate_params.to_csv(REPORTS_DIR / "separate_sigma_parameters.csv", index=False)
    equal_params.to_csv(REPORTS_DIR / "equal_sigma_parameters.csv", index=False)
    inverse_params.to_csv(REPORTS_DIR / "inverse_power_law_parameters.csv", index=False)
    inverse_params_without_361.to_csv(REPORTS_DIR / "inverse_power_law_parameters_excluding_361_4.csv", index=False)
    model_summary.to_csv(REPORTS_DIR / "model_summary.csv", index=False)
    prediction_all.to_csv(REPORTS_DIR / "prediction_50kv.csv", index=False)
    prediction_without_361.to_csv(REPORTS_DIR / "prediction_50kv_excluding_361_4.csv", index=False)

    # Original R-style visualization set. These are the main portfolio figures.
    plot_01_raw_scatter(df, FIGURES_DIR / "01_r_raw_scatter_voltage_time.png")
    plot_02_log_scatter_with_regression(df, FIGURES_DIR / "02_r_log_scatter_with_included_excluded_regression.png")
    plot_03_probability_plot_separate_sigma(df, separate_params, FIGURES_DIR / "03_r_probability_plot_separate_sigma.png")
    plot_04_probability_plot_equal_sigma(df, equal_params, FIGURES_DIR / "04_r_probability_plot_equal_sigma.png")
    plot_05_probability_plot_overlay_sigma_assumptions(
        df,
        separate_params,
        equal_params,
        FIGURES_DIR / "05_r_probability_plot_overlay_separate_and_equal_sigma.png",
    )
    plot_06_probability_plot_inverse_power_with_50kv_ci(
        df,
        inverse_summary,
        FIGURES_DIR / "06_r_probability_plot_50kv_ci_all_stress_levels.png",
        time_grid_start=1_000,
        time_grid_end=10_000_000,
        x_min=0.1,
        x_max=150_000,
    )
    plot_07_stress_life_quantiles(
        df,
        inverse_summary,
        FIGURES_DIR / "07_r_stress_life_quantiles_all_stress_levels.png",
    )
    inverse_power_residual_table = plot_08_standardized_residuals_inverse_power(
        df,
        inverse_summary,
        FIGURES_DIR / "08_r_standardized_residuals_inverse_power.png",
    )
    equal_residuals = plot_09_residual_probability_plot_equal_sigma(
        df,
        equal_params,
        FIGURES_DIR / "09_r_residual_probability_plot_equal_sigma.png",
    )
    plot_06_probability_plot_inverse_power_with_50kv_ci(
        df_without_361,
        inverse_summary_without_361,
        FIGURES_DIR / "10_r_probability_plot_50kv_ci_excluding_361_4.png",
        time_grid_start=10,
        time_grid_end=2_000_000,
        x_min=0.1,
        x_max=150_000,
    )
    plot_07_stress_life_quantiles(
        df_without_361,
        inverse_summary_without_361,
        FIGURES_DIR / "11_r_stress_life_quantiles_excluding_361_4.png",
    )
    plot_12_stress_life_quantile_comparison(
        df,
        inverse_summary,
        inverse_summary_without_361,
        FIGURES_DIR / "12_stress_life_quantile_comparison_including_vs_excluding_361_4.png",
    )

    # Supplemental figure: not in the original R script, but useful for GitHub/README.
    plot_supplemental_aic_bar(model_summary, FIGURES_DIR / "supplemental_model_comparison_aic.png")

    inverse_power_residual_table.to_csv(REPORTS_DIR / "inverse_power_law_residuals.csv", index=False)
    equal_residuals.to_csv(REPORTS_DIR / "equal_sigma_residuals.csv", index=False)

    # Figure catalog for README/GitHub.
    figure_catalog = pd.DataFrame(
        [
            ["01_r_raw_scatter_voltage_time.png", "Original R-style raw scatter plot before log transformation."],
            ["02_r_log_scatter_with_included_excluded_regression.png", "Log-log scatter plot with four-side axes and regression lines with/without 361.4 kV/mm."],
            ["03_r_probability_plot_separate_sigma.png", "Lognormal probability plot with each voltage group having its own sigma."],
            ["04_r_probability_plot_equal_sigma.png", "Lognormal probability plot under the common-sigma assumption."],
            ["05_r_probability_plot_overlay_separate_and_equal_sigma.png", "Overlay version following the R script if both line sets are drawn on one canvas."],
            ["06_r_probability_plot_50kv_ci_all_stress_levels.png", "Inverse power law probability plot with 50 kV/mm extrapolation and confidence band using all stress levels."],
            ["07_r_stress_life_quantiles_all_stress_levels.png", "Stress-life plot with 10%, 50%, and 90% quantile lines and lognormal density traces using all stress levels."],
            ["08_r_standardized_residuals_inverse_power.png", "R-style standardized residual plot for the inverse power law model, using sigma estimated from the inverse power law fit."],
            ["09_r_residual_probability_plot_equal_sigma.png", "Probability plot of equal-sigma exponential standardized residuals."],
            ["10_r_probability_plot_50kv_ci_excluding_361_4.png", "Same as Figure 06 after excluding 361.4 kV/mm."],
            ["11_r_stress_life_quantiles_excluding_361_4.png", "Same as Figure 07 after excluding 361.4 kV/mm."],
            ["12_stress_life_quantile_comparison_including_vs_excluding_361_4.png", "Presentation-focused comparison of 10%, 50%, and 90% lifetime quantile curves with and without the 361.4 kV/mm stress group."],
            ["supplemental_model_comparison_aic.png", "Supplemental GitHub-friendly AIC bar chart, added for quick model comparison."],
        ],
        columns=["figure", "description"],
    )
    figure_catalog.to_csv(REPORTS_DIR / "figure_catalog.csv", index=False)

    print("\n=== Model comparison ===")
    print(model_summary.round(4).to_string(index=False))

    print("\n=== LRT: Equal sigma vs Separate sigma ===")
    print({key: round(value, 4) if isinstance(value, float) else value for key, value in lrt_equal_vs_separate.items()})

    print("\n=== LRT: Inverse power law vs Equal sigma ===")
    print({key: round(value, 4) if isinstance(value, float) else value for key, value in lrt_inverse_vs_equal.items()})

    print("\n=== Prediction at 50 kV/mm, all stress levels ===")
    print(prediction_all.round(4).to_string(index=False))

    print("\n=== Prediction at 50 kV/mm, excluding 361.4 kV/mm ===")
    print(prediction_without_361.round(4).to_string(index=False))

    print("\n=== Inverse-power-law residual outliers, abs(log residual) > 2 ===")
    print(
        inverse_power_residual_table.loc[
            inverse_power_residual_table["is_outlier_abs_log_residual_gt_2"],
            ["voltage_kv_mm", "failure_time", "fitted_value", "log_standardized_residual"],
        ].round(4).to_string(index=False)
    )

    print("\nAnalysis completed.")
    print(f"Reports saved to: {REPORTS_DIR}")
    print(f"Figures saved to: {FIGURES_DIR}")


if __name__ == "__main__":
    main()
