import math


def sigma_bias_correction(n: int) -> float:
    """
    Bias-correction multiplier used in the original R project for
    the standard error of the MLE sigma estimator.
    """
    return math.sqrt(2 / n) * math.gamma(n / 2) / math.gamma((n - 1) / 2)


def normal_ci(estimate: float, standard_error: float, z: float = 1.959963984540054) -> tuple[float, float]:
    """Return Wald-type confidence interval."""
    return estimate - z * standard_error, estimate + z * standard_error
