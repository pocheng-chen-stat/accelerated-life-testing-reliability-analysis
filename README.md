# Accelerated Life Testing Reliability Analysis

## Project Overview

This project analyzes accelerated life testing (ALT) data for electrical insulation lifetime under multiple voltage stress levels. The goal is to estimate lifetime behavior under accelerated stress conditions, evaluate model assumptions, and extrapolate reliability performance to a lower normal-use voltage level.

The original analysis was developed in R and has been refactored into a reproducible Python project with modular source code, automated tests, model outputs, and R-style reliability visualizations.

## Engineering Problem

In reliability engineering, products are often tested under high-stress conditions so that failures can be observed within a reasonable amount of time. However, the engineering objective is usually not only to describe high-stress failures, but also to infer reliability under normal-use conditions.

In this project, insulation lifetime data are observed at several voltage stress levels:

```text
100.3, 122.4, 157.1, 219.0, and 361.4 kV/mm
```

The analysis aims to answer three questions:

1. What lifetime distribution is appropriate at each observed voltage level?
2. Can the voltage groups share a common dispersion parameter?
3. How can the model be used to extrapolate lifetime behavior at 50 kV/mm?

## Dataset

The dataset contains failure times from an accelerated life testing experiment.

| Column | Description |
|---|---|
| `unit_id` | Unit identifier |
| `voltage_kv_mm` | Voltage stress level in kV/mm |
| `failure_time` | Observed failure time |
| `event` | Failure indicator; all observations are treated as observed failures in this project |

This project assumes no right censoring, consistent with the original analysis.

## Modeling Strategy

A key point in this project is that the modeling workflow has two layers:

1. **Lifetime distribution model**  
   Describes the lifetime distribution at each observed stress level.

2. **Acceleration model**  
   Describes how lifetime changes as the voltage stress level changes, allowing extrapolation to an unobserved voltage level.

These two parts should not be interpreted as the same modeling task.

## 1. Lifetime Distribution Model

The lifetime at each fixed voltage level is modeled using a Lognormal distribution:

```text
T | V_j ~ Lognormal(mu_j, sigma_j)
```

Two alternatives are considered.

### Separate Sigma Model

Each voltage level has its own mean and standard deviation on the log-time scale:

```text
T | V_j ~ Lognormal(mu_j, sigma_j)
```

This model is flexible, but it uses more parameters.

### Equal Sigma Model

Each voltage level has its own mean, but all voltage levels share a common standard deviation:

```text
T | V_j ~ Lognormal(mu_j, sigma)
```

This model asks whether the dispersion of log lifetime can be treated as common across stress levels.

The equal-sigma assumption is evaluated using likelihood-based model comparison, including AIC and likelihood ratio testing.

## 2. Acceleration Model

The equal-sigma model can describe the observed voltage levels, but it does not directly provide a value of `mu` at an unobserved voltage such as 50 kV/mm.

To extrapolate to normal-use voltage, an acceleration relationship is required. This project uses an inverse power law model:

```text
log(T) = beta0 + beta1 * log(V) + error
```

or equivalently:

```text
mu(V) = beta0 + beta1 * log(V)
```

This model enables prediction at unobserved voltage levels such as 50 kV/mm.

Important interpretation:

- The equal-sigma model is used to evaluate the lifetime distribution structure across observed stress levels.
- The inverse power law model is used to establish a stress-life relationship for extrapolation.
- The inverse power law model should be interpreted carefully if it is not preferred by AIC relative to the observed-data equal-sigma model.

## Model Comparison

The project compares three model forms:

| Model | Purpose | Number of Parameters |
|---|---|---:|
| Separate sigma | Flexible Lognormal fit for each voltage group | 10 |
| Equal sigma | Lognormal fit with common dispersion across voltage groups | 6 |
| Inverse power law | Stress-life relationship for extrapolation | 3 |

Based on the original analysis, the equal-sigma model has the lowest AIC among the candidate models for the observed data.

However, the equal-sigma model alone cannot extrapolate to 50 kV/mm because it estimates separate means only for the voltage levels observed in the data. Therefore, the inverse power law model is used specifically as an extrapolation model.

## Key Results

The Python implementation reproduces the main numerical results of the original R analysis.

Expected model comparison values are approximately:

| Model | -2 Log-Likelihood | AIC |
|---|---:|---:|
| Separate sigma | 564.9022 | 584.9022 |
| Equal sigma | 567.2105 | 579.2105 |
| Inverse power law | 579.7550 | 585.7550 |

The equal-sigma model has the lowest AIC for the observed accelerated life testing data.

For extrapolation to 50 kV/mm using the inverse power law model, the expected values from the full dataset are approximately:

| Quantity | Estimate |
|---|---:|
| `F(10000)` | 0.0028 |
| `t_0.1` | 58541.8 |
| `t_0.5` | 267657.6 |

The project also includes a sensitivity analysis excluding the highest voltage level, 361.4 kV/mm, because this stress level has strong influence on the inverse power law extrapolation.

## Visualizations

This project preserves the R-style reliability plots from the original analysis. Many of the plots are drawn on transformed coordinates while displaying original-scale labels for readability.

For example:

- Actual plotting coordinate: `log(kV/mm)` or `log(time)`
- Bottom / left axes: original-scale voltage, time, or failure probability
- Top / right axes: transformed-scale values such as `log(kV/mm)`, `log(time)`, or standard normal quantiles

This design is useful because reliability probability plots are often easier to fit visually on transformed scales, while engineering interpretation is easier on the original scale.

### Main Figures

| Figure | Purpose |
|---|---|
| `01_r_raw_scatter_voltage_time.png` | Raw voltage versus failure time scatter plot |
| `02_r_log_scatter_with_included_excluded_regression.png` | Log-log scatter plot with regression lines including and excluding 361.4 kV/mm |
| `03_r_probability_plot_separate_sigma.png` | Lognormal probability plot with separate sigma fits |
| `04_r_probability_plot_equal_sigma.png` | Lognormal probability plot with equal sigma fits |
| `05_r_probability_plot_overlay_separate_and_equal_sigma.png` | Comparison of separate-sigma and equal-sigma probability plot fits |
| `06_r_probability_plot_50kv_ci_all_stress_levels.png` | Probability plot with 50 kV/mm extrapolation and confidence bands |
| `07_r_stress_life_quantiles_all_stress_levels.png` | Stress-life plot with lifetime quantile lines |
| `08_r_standardized_residuals_inverse_power.png` | Standardized residual diagnostics for the inverse power law model |
| `09_r_residual_probability_plot_equal_sigma.png.png` | Probability plot of residuals based on the equal-sigma model |
| `10_r_probability_plot_50kv_ci_excluding_361_4.png` | 50 kV/mm extrapolation after excluding 361.4 kV/mm |
| `11_r_stress_life_quantiles_excluding_361_4.png` | Stress-life quantile plot after excluding 361.4 kV/mm |
| `supplemental_model_comparison_aic.png` | Supplemental AIC comparison chart for quick GitHub presentation |

The supplemental AIC chart is added for presentation clarity. The core reliability plots are based on the original R-style analysis.

## Residual Diagnostics

The standardized residual plot for the inverse power law model uses the inverse-power-law model's own estimated sigma:

```text
z_i = [log(t_i) - beta0 - beta1 * log(v_i)] / sigma_IPL
```

This is used because the residuals are generated from the inverse power law model itself. The plot keeps the original R-style visual design, including transformed axes and red markers for large residuals.

## How to Run

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the full analysis:

```bash
python main.py
```

Run tests:

```bash
python -m pytest
```

The output files will be saved to:

```text
reports/
reports/figures/
```

## Repository Structure

```text
accelerated-life-testing-reliability-analysis/
│
├── README.md
├── requirements.txt
├── main.py
│
├── data/
│   └── insulation_lifetime_data.csv
│
├── src/
│   ├── data_loader.py
│   ├── lognormal_models.py
│   ├── inverse_power_law.py
│   ├── model_selection.py
│   ├── r_style_visualizations.py
│   └── utils.py
│
├── tests/
│   └── test_reliability_results.py
│
└── reports/
    ├── model_summary.csv
    ├── prediction_50kv.csv
    ├── inverse_power_law_residuals.csv
    └── figures/
```

## Engineering Interpretation

The analysis supports the following interpretation:

1. The equal-sigma Lognormal model provides a parsimonious fit to the observed accelerated life testing data.
2. The common sigma assumption is useful because it reduces model complexity while retaining a reasonable fit across voltage groups.
3. The inverse power law model is necessary for extrapolation because the equal-sigma model does not define lifetime behavior at unobserved voltage levels.
4. The 50 kV/mm prediction should be interpreted as model-based extrapolation, not direct empirical evidence.
5. Excluding the highest stress level can materially change extrapolated predictions, indicating sensitivity to high-stress observations.

## Limitations

- The data are treated as complete failure observations with no right censoring.
- The inverse power law model may not be the best-fitting model for the observed data based on AIC.
- Extrapolation to 50 kV/mm is sensitive to model assumptions and influential stress levels.
- Confidence bands are based on asymptotic approximations.
- Additional physical knowledge would be needed before using the extrapolated result for actual engineering decisions.

## Future Improvements

Potential extensions include:

- Add right-censored data support
- Compare Lognormal and Weibull lifetime distributions
- Add Arrhenius or Eyring acceleration models
- Add bootstrap confidence intervals for extrapolated reliability estimates
- Add formal influence diagnostics for high-stress observations
- Build an interactive dashboard for reliability prediction

## Skills Demonstrated

This project demonstrates:

- Reliability engineering analysis
- Accelerated life testing modeling
- Maximum likelihood estimation
- Lognormal lifetime modeling
- AIC and likelihood ratio testing
- Inverse power law acceleration modeling
- Model-based extrapolation
- Residual diagnostics
- Python project refactoring from R
- Reproducible GitHub project organization
