import numpy as np
import pandas as pd
from scipy import stats

from src.utils import normal_ci, sigma_bias_correction


def fit_lognormal_separate_sigma(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """
    Fit one Lognormal distribution for each voltage level.

    Each voltage group has its own mu and sigma.
    This corresponds to the 'SepDists' model in the original R analysis.
    """
    rows = []

    for voltage, group in df.groupby("voltage_kv_mm", sort=True):
        failure_times = group["failure_time"].to_numpy(dtype=float)
        log_times = np.log(failure_times)
        n = len(failure_times)

        mu_hat = float(np.mean(log_times))
        sigma_hat = float(np.sqrt(np.mean((log_times - mu_hat) ** 2)))

        log_likelihood = float(
            np.sum(stats.lognorm.logpdf(failure_times, s=sigma_hat, scale=np.exp(mu_hat)))
        )

        mu_se = float(sigma_hat / np.sqrt(n))
        sigma_se = float(np.sqrt(sigma_hat**2 / (2 * n)) * sigma_bias_correction(n))

        mu_lower, mu_upper = normal_ci(mu_hat, mu_se)
        sigma_lower, sigma_upper = normal_ci(sigma_hat, sigma_se)

        rows.append(
            {
                "model": "Separate sigma",
                "voltage_kv_mm": voltage,
                "parameter": "mu",
                "n": n,
                "mle": mu_hat,
                "se": mu_se,
                "lower_95": mu_lower,
                "upper_95": mu_upper,
                "log_likelihood": log_likelihood,
            }
        )
        rows.append(
            {
                "model": "Separate sigma",
                "voltage_kv_mm": voltage,
                "parameter": "sigma",
                "n": n,
                "mle": sigma_hat,
                "se": sigma_se,
                "lower_95": sigma_lower,
                "upper_95": sigma_upper,
                "log_likelihood": np.nan,
            }
        )

    parameter_table = pd.DataFrame(rows)
    total_log_likelihood = float(parameter_table["log_likelihood"].sum(skipna=True))
    num_groups = df["voltage_kv_mm"].nunique()
    summary = {
        "model": "Separate sigma",
        "log_likelihood": total_log_likelihood,
        "num_params": int(2 * num_groups),
    }

    return parameter_table, summary


def fit_lognormal_equal_sigma(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """
    Fit voltage-specific mu values with one common sigma across all voltage groups.

    This corresponds to the 'EqualSig' model in the original R analysis.
    """
    group_rows = []

    for voltage, group in df.groupby("voltage_kv_mm", sort=True):
        log_times = np.log(group["failure_time"].to_numpy(dtype=float))
        group_rows.append(
            {
                "voltage_kv_mm": voltage,
                "n": len(log_times),
                "mu": float(np.mean(log_times)),
            }
        )

    group_table = pd.DataFrame(group_rows)
    mu_map = dict(zip(group_table["voltage_kv_mm"], group_table["mu"]))

    residuals = []
    for _, row in df.iterrows():
        residuals.append(np.log(row["failure_time"]) - mu_map[row["voltage_kv_mm"]])

    sigma_hat = float(np.sqrt(np.mean(np.asarray(residuals) ** 2)))

    total_log_likelihood = 0.0
    for voltage, group in df.groupby("voltage_kv_mm", sort=True):
        failure_times = group["failure_time"].to_numpy(dtype=float)
        mu_hat = mu_map[voltage]
        total_log_likelihood += float(
            np.sum(stats.lognorm.logpdf(failure_times, s=sigma_hat, scale=np.exp(mu_hat)))
        )

    n_all = len(df)
    rows = []

    for _, row in group_table.iterrows():
        mu_hat = row["mu"]
        mu_se = float(sigma_hat / np.sqrt(row["n"]))
        lower, upper = normal_ci(mu_hat, mu_se)
        rows.append(
            {
                "model": "Equal sigma",
                "voltage_kv_mm": row["voltage_kv_mm"],
                "parameter": "mu",
                "n": int(row["n"]),
                "mle": mu_hat,
                "se": mu_se,
                "lower_95": lower,
                "upper_95": upper,
                "log_likelihood": np.nan,
            }
        )

    sigma_se = float(np.sqrt(sigma_hat**2 / (2 * n_all)) * sigma_bias_correction(n_all))
    lower, upper = normal_ci(sigma_hat, sigma_se)

    rows.append(
        {
            "model": "Equal sigma",
            "voltage_kv_mm": np.nan,
            "parameter": "sigma",
            "n": n_all,
            "mle": sigma_hat,
            "se": sigma_se,
            "lower_95": lower,
            "upper_95": upper,
            "log_likelihood": total_log_likelihood,
        }
    )

    parameter_table = pd.DataFrame(rows)
    summary = {
        "model": "Equal sigma",
        "log_likelihood": float(total_log_likelihood),
        "num_params": int(df["voltage_kv_mm"].nunique() + 1),
    }

    return parameter_table, summary
