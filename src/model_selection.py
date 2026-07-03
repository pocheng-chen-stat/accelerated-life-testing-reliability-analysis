from scipy import stats


def calculate_aic(log_likelihood: float, num_params: int) -> float:
    """Akaike Information Criterion."""
    return float(-2 * log_likelihood + 2 * num_params)


def likelihood_ratio_test(ll_reduced: float, ll_full: float, df: int) -> dict:
    """Likelihood ratio test for nested models."""
    statistic = float(-2 * (ll_reduced - ll_full))
    p_value = float(1 - stats.chi2.cdf(statistic, df))
    critical_value_95 = float(stats.chi2.ppf(0.95, df))
    return {
        "test_statistic": statistic,
        "df": df,
        "p_value": p_value,
        "critical_value_95": critical_value_95,
    }
