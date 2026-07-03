from pathlib import Path
import pandas as pd


REQUIRED_COLUMNS = {"unit_id", "voltage_kv_mm", "failure_time", "event"}


def load_lifetime_data(path: str | Path) -> pd.DataFrame:
    """Load insulation lifetime data."""
    df = pd.read_csv(path)
    validate_lifetime_data(df)
    return df


def validate_lifetime_data(df: pd.DataFrame) -> None:
    """Validate required columns and basic data quality."""
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    if (df["voltage_kv_mm"] <= 0).any():
        raise ValueError("voltage_kv_mm must be positive.")

    if (df["failure_time"] <= 0).any():
        raise ValueError("failure_time must be positive.")

    if not set(df["event"].unique()).issubset({0, 1}):
        raise ValueError("event must be coded as 0 or 1.")
