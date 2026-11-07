"""Model definitions and leave-one-chromosome-out evaluation."""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .config import (FEATURES, LR_C, RF_MIN_SAMPLES_LEAF, RF_N_ESTIMATORS, SEED)


def make_lr() -> Pipeline:
    return Pipeline([
        ("scaler", StandardScaler()),
        ("lr", LogisticRegression(C=LR_C, max_iter=1000, solver="lbfgs")),
    ])


def make_rf() -> RandomForestClassifier:
    return RandomForestClassifier(
        n_estimators=RF_N_ESTIMATORS,
        min_samples_leaf=RF_MIN_SAMPLES_LEAF,
        max_features="sqrt",
        class_weight="balanced",
        random_state=SEED,
        n_jobs=-1,
    )


MODEL_FACTORIES = {"logistic_regression": make_lr, "random_forest": make_rf}


def design_matrix(df: pd.DataFrame, features: list[str] | None = None):
    features = features or FEATURES
    return df[features].to_numpy(dtype=float), df["label"].to_numpy(dtype=int)


def logo_oof(df: pd.DataFrame, model_name: str,
             features: list[str] | None = None):
    """Out-of-fold predictions with LeaveOneGroupOut over chromosomes.

    Returns (y_true, y_score, fold_chromosome). Every chromosome is held out
    exactly once, so no window from the test chromosome is seen in training.
    """
    X, y = design_matrix(df, features)
    groups = df["chr"].to_numpy()
    oof = np.full(len(df), np.nan)
    fold = np.empty(len(df), dtype=object)
    logo = LeaveOneGroupOut()
    for train_idx, test_idx in logo.split(X, y, groups):
        model = MODEL_FACTORIES[model_name]()
        model.fit(X[train_idx], y[train_idx])
        oof[test_idx] = model.predict_proba(X[test_idx])[:, 1]
        fold[test_idx] = groups[test_idx]
    return y, oof, fold


def fit_full(df: pd.DataFrame, model_name: str):
    """Fit on the whole dataset for coefficient / importance reporting."""
    X, y = design_matrix(df)
    model = MODEL_FACTORIES[model_name]()
    model.fit(X, y)
    return model


def lr_coefficients(model: Pipeline) -> pd.Series:
    """Logistic-regression coefficients on standardised features (log-odds)."""
    coefs = model.named_steps["lr"].coef_[0]
    return pd.Series(coefs, index=FEATURES, name="coefficient")


def rf_importances(model: RandomForestClassifier) -> pd.Series:
    return pd.Series(model.feature_importances_, index=FEATURES, name="importance")
