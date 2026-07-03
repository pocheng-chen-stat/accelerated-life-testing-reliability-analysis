import pandas as pd

from src.lognormal_models import fit_lognormal_equal_sigma, fit_lognormal_separate_sigma
from src.inverse_power_law import fit_inverse_power_law, predict_failure_probability, predict_lifetime_quantile
from src.model_selection import calculate_aic


def test_key_results_match_original_r_analysis():
    df = pd.read_csv("data/insulation_lifetime_data.csv")

    _, sep = fit_lognormal_separate_sigma(df)
    _, eq = fit_lognormal_equal_sigma(df)
    _, reg = fit_inverse_power_law(df)

    assert round(sep["log_likelihood"], 4) == -282.4511
    assert round(eq["log_likelihood"], 4) == -283.6053
    assert round(reg["log_likelihood"], 4) == -289.8775

    assert round(calculate_aic(sep["log_likelihood"], sep["num_params"]), 4) == 584.9022
    assert round(calculate_aic(eq["log_likelihood"], eq["num_params"]), 4) == 579.2105
    assert round(calculate_aic(reg["log_likelihood"], reg["num_params"]), 4) == 585.7550

    assert round(predict_failure_probability(reg, 50, 10000), 4) == 0.0028
    assert round(predict_lifetime_quantile(reg, 50, 0.1), 1) == 58541.8
    assert round(predict_lifetime_quantile(reg, 50, 0.5), 1) == 267657.6
