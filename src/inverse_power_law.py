import numpy as np
import pandas as pd
from scipy import stats

from src.utils import normal_ci, sigma_bias_correction


def _legacy_r_covariance(voltage: np.ndarray, y: np.ndarray, beta0: float, beta1: float, sigma: float) -> np.ndarray:
    """
    Covariance matrix based on the observed information used by the original
    R project, corrected to match the fitted predictor x = log(voltage).

    The earlier implementation accidentally used raw voltage in two information
    matrix terms for beta1.  Because the model is
    log(T) = beta0 + beta1 * log(V) + error, those terms must use log(V).
    A standard regression-based covariance estimate is also available through
    `_standard_covariance`.
    """
    n = len(y)
    # The fitted regression uses x = log(voltage), so every derivative with
    # respect to beta1 must also use log(voltage).  The previous implementation
    # mixed the raw voltage into the beta1-beta1 and beta1-sigma information
    # terms, which severely understated Var(beta1).
    x = np.log(voltage)
    residual = y - beta0 - beta1 * x

    values = [
        n / sigma**2,
        np.sum(x) / sigma**2,
        2 * np.sum(residual) / sigma**3,
        np.sum(x) / sigma**2,
        np.sum(x**2) / sigma**2,
        2 * np.sum(residual * x) / sigma**3,
        2 * np.sum(residual) / sigma**3,
        2 * np.sum(residual * x) / sigma**3,
        2 * n / sigma**2,
    ]

    information = np.array(values, dtype=float).reshape((3, 3), order="F")
    information_inverse = np.linalg.inv(information)

    correction = np.diag([1.0, 1.0, sigma_bias_correction(n)])
    return correction @ information_inverse @ correction


def _standard_covariance(voltage: np.ndarray, y: np.ndarray, beta0: float, beta1: float, sigma: float) -> np.ndarray:
    """
    Standard large-sample covariance for log(T)=beta0+beta1 log(V)+error.

    This is statistically cleaner than the legacy matrix, but the legacy matrix
    is retained to reproduce the original course project figures.
    """
    x = np.log(voltage)
    n = len(y)
    X = np.column_stack([np.ones(n), x])
    beta_cov = sigma**2 * np.linalg.inv(X.T @ X)
    sigma_var = sigma**2 / (2 * n) * sigma_bias_correction(n) ** 2

    covariance = np.zeros((3, 3), dtype=float)
    covariance[:2, :2] = beta_cov
    covariance[2, 2] = sigma_var
    return covariance


def fit_inverse_power_law(
    df: pd.DataFrame,
    *,
    label: str = "Inverse power law",
    covariance_method: str = "legacy_r",
) -> tuple[pd.DataFrame, dict]:
    """
    Fit inverse power law model by MLE.

    Model:
        log(T) = beta0 + beta1 * log(V) + error
        error ~ Normal(0, sigma^2)

    Because T is Lognormal, the likelihood includes the lognormal Jacobian term.
    """
    voltage = df["voltage_kv_mm"].to_numpy(dtype=float)
    failure_time = df["failure_time"].to_numpy(dtype=float)

    x = np.log(voltage)
    y = np.log(failure_time)
    n = len(y)

    X = np.column_stack([np.ones(n), x])
    beta_hat = np.linalg.lstsq(X, y, rcond=None)[0]
    beta0_hat, beta1_hat = float(beta_hat[0]), float(beta_hat[1])
    fitted_mu = X @ beta_hat
    residuals = y - fitted_mu

    # MLE denominator is n, matching the original R script.
    sigma_hat = float(np.sqrt(np.mean(residuals**2)))

    log_likelihood = float(
        np.sum(stats.norm.logpdf(y, loc=fitted_mu, scale=sigma_hat))
        - np.sum(np.log(failure_time))
    )

    if covariance_method == "legacy_r":
        covariance = _legacy_r_covariance(voltage, y, beta0_hat, beta1_hat, sigma_hat)
    elif covariance_method == "standard":
        covariance = _standard_covariance(voltage, y, beta0_hat, beta1_hat, sigma_hat)
    else:
        raise ValueError("covariance_method must be 'legacy_r' or 'standard'.")

    parameter_names = ["beta0", "beta1", "sigma"]
    estimates = [beta0_hat, beta1_hat, sigma_hat]
    standard_errors = np.sqrt(np.diag(covariance))

    rows = []
    for name, estimate, se in zip(parameter_names, estimates, standard_errors):
        lower, upper = normal_ci(estimate, float(se))
        rows.append(
            {
                "model": label,
                "parameter": name,
                "mle": estimate,
                "se": float(se),
                "lower_95": lower,
                "upper_95": upper,
                "log_likelihood": log_likelihood if name == "beta0" else np.nan,
            }
        )

    parameter_table = pd.DataFrame(rows)
    summary = {
        "model": label,
        "beta0": beta0_hat,
        "beta1": beta1_hat,
        "sigma": sigma_hat,
        "log_likelihood": log_likelihood,
        "num_params": 3,
        "n": n,
        "covariance_method": covariance_method,
        "covariance": covariance,
    }

    return parameter_table, summary


def predict_failure_probability(summary: dict, voltage_kv_mm: float, time: float) -> float:
    """Return F(time | voltage), the predicted probability of failure by a given time."""
    beta0 = summary["beta0"]
    beta1 = summary["beta1"]
    sigma = summary["sigma"]
    mu = beta0 + beta1 * np.log(voltage_kv_mm)
    z = (np.log(time) - mu) / sigma
    return float(stats.norm.cdf(z))


def predict_lifetime_quantile(summary: dict, voltage_kv_mm: float, probability: float) -> float:
    """Return t_p such that P(T <= t_p) = probability."""
    beta0 = summary["beta0"]
    beta1 = summary["beta1"]
    sigma = summary["sigma"]
    mu = beta0 + beta1 * np.log(voltage_kv_mm)
    return float(np.exp(mu + stats.norm.ppf(probability) * sigma))


def failure_probability_confidence_band(
    summary: dict,
    voltage_kv_mm: float,
    time_grid: np.ndarray,
    confidence: float = 0.95,
) -> pd.DataFrame:
    """
    Confidence band for F(t | voltage) using the transformation in the R project.

    The calculation follows the original script:
        xi = (log(t)-mu)/sigma
        se_F = phi(xi)/sigma * sqrt(var_mu + 2 xi cov_mu_sigma + xi^2 var_sigma)
        w = exp(z * se_F / (F * (1-F)))
    """
    beta0 = summary["beta0"]
    beta1 = summary["beta1"]
    sigma = summary["sigma"]
    covariance = np.asarray(summary["covariance"], dtype=float)

    log_v = np.log(voltage_kv_mm)
    mu = beta0 + beta1 * log_v
    xi = (np.log(time_grid) - mu) / sigma
    failure_prob = stats.norm.cdf(xi)

    # Delta-method covariance of (mu_at_voltage, sigma).
    var_mu = (
        covariance[0, 0]
        + 2 * log_v * covariance[0, 1]
        + log_v**2 * covariance[1, 1]
    )
    cov_mu_sigma = covariance[0, 2] + log_v * covariance[1, 2]
    var_sigma = covariance[2, 2]

    se_f = (stats.norm.pdf(xi) / sigma) * np.sqrt(
        np.maximum(var_mu + 2 * xi * cov_mu_sigma + xi**2 * var_sigma, 0)
    )

    z = stats.norm.ppf((1 + confidence) / 2)
    eps = np.finfo(float).eps
    f_safe = np.clip(failure_prob, eps, 1 - eps)
    w = np.exp(z * se_f / (f_safe * (1 - f_safe)))

    lower = failure_prob / (failure_prob + (1 - failure_prob) * w)
    upper = failure_prob / (failure_prob + (1 - failure_prob) / w)

    return pd.DataFrame(
        {
            "time": time_grid,
            "failure_probability": failure_prob,
            "lower_95": lower,
            "upper_95": upper,
            "standard_quantile": stats.norm.ppf(failure_prob),
            "lower_95_quantile": stats.norm.ppf(lower),
            "upper_95_quantile": stats.norm.ppf(upper),
        }
    )
