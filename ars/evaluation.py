"""Metrics, bootstrap confidence intervals and feature-level tests."""

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import average_precision_score, roc_auc_score

from .config import FEATURES, N_BOOTSTRAP, SEED


def point_metrics(y_true, y_score) -> dict[str, float]:
    return {
        "roc_auc": float(roc_auc_score(y_true, y_score)),
        "pr_auc": float(average_precision_score(y_true, y_score)),
        "n_pos": int(np.sum(y_true == 1)),
        "n_neg": int(np.sum(y_true == 0)),
    }


def bootstrap_ci(y_true, y_score, n_boot: int = N_BOOTSTRAP,
                 seed: int = SEED) -> dict[str, tuple[float, float]]:
    """Stratified bootstrap 95% CI for ROC-AUC and PR-AUC.

    Positives and negatives are resampled separately so every replicate keeps
    the original class ratio.
    """
    rng = np.random.default_rng(seed)
    pos_idx = np.flatnonzero(y_true == 1)
    neg_idx = np.flatnonzero(y_true == 0)
    roc_vals, pr_vals = [], []
    for _ in range(n_boot):
        pos = rng.choice(pos_idx, size=len(pos_idx), replace=True)
        neg = rng.choice(neg_idx, size=len(neg_idx), replace=True)
        idx = np.concatenate([pos, neg])
        yt, ys = y_true[idx], y_score[idx]
        if yt.sum() == 0 or (yt == 0).sum() == 0:
            continue
        roc_vals.append(roc_auc_score(yt, ys))
        pr_vals.append(average_precision_score(yt, ys))

    def ci(values):
        values = np.asarray(values)
        return float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))

    return {"roc_auc_ci": ci(roc_vals), "pr_auc_ci": ci(pr_vals)}


def bh_adjust(p_values: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR adjustment."""
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(n)
    adjusted[order] = np.clip(ranked, 0, 1)
    return adjusted.tolist()


def feature_tests(df: pd.DataFrame, negative_type: str) -> pd.DataFrame:
    """Mann-Whitney U tests for each feature: positives vs one negative set."""
    positives = df[df["negative_type"] == "positive"]
    negatives = df[df["negative_type"] == negative_type]
    rows = []
    for feature in FEATURES:
        x = positives[feature].to_numpy()
        y = negatives[feature].to_numpy()
        stat, p_value = mannwhitneyu(x, y, alternative="two-sided")
        # Cliff's delta: +1 means positives are entirely above negatives.
        delta = 2 * stat / (len(x) * len(y)) - 1
        rows.append({
            "feature": feature,
            "negative_type": negative_type,
            "pos_median": float(np.median(x)),
            "neg_median": float(np.median(y)),
            "cliffs_delta": float(delta),
            "p_value": float(p_value),
            "n_pos": len(x),
            "n_neg": len(y),
        })
    out = pd.DataFrame(rows)
    out["p_fdr_bh"] = bh_adjust(out["p_value"].tolist())
    return out


def summarize_dataset(df: pd.DataFrame) -> dict:
    return {
        "n_windows": int(len(df)),
        "n_positives": int((df["label"] == 1).sum()),
        "n_random_negatives": int(
            ((df["label"] == 0) & (df["negative_type"] == "random")).sum()),
        "n_matched_negatives": int(
            ((df["label"] == 0) & (df["negative_type"] == "distance_matched")).sum()),
        "chromosomes": sorted(df["chr"].unique(), key=int),
    }


def ppv_table(oof_frames: list[pd.DataFrame], genome_bp: int, n_pos: int,
              window: int, prevalences: tuple[float, ...] = ()) -> pd.DataFrame:
    """Convert sensitivity/specificity at threshold 0.5 into PPV.

    The model is trained on balanced data, where precision equals PPV at 50%
    prevalence. At the genome-wide base rate the same model has a much lower
    PPV: this is the base-rate lesson from the PSA example in the textbook.
    """
    prevalence_values = {
        "genome_wide": n_pos * window / genome_bp,
        "psa_analogy_5.8pct": 0.058,
        "balanced_50pct": 0.50,
    }
    for value in prevalences:
        prevalence_values[f"custom_{value:.3f}"] = value

    rows = []
    for frame in oof_frames:
        y = frame["label"].to_numpy()
        neg_type = frame["negative_type"].iloc[0]
        for model, column in (("logistic_regression", "lr_score"),
                              ("random_forest", "rf_score")):
            pred = (frame[column].to_numpy() >= 0.5).astype(int)
            tp = int(((pred == 1) & (y == 1)).sum())
            fn = int(((pred == 0) & (y == 1)).sum())
            tn = int(((pred == 0) & (y == 0)).sum())
            fp = int(((pred == 1) & (y == 0)).sum())
            sensitivity = tp / (tp + fn)
            specificity = tn / (tn + fp)
            for label, prevalence in prevalence_values.items():
                denom = sensitivity * prevalence + (1 - specificity) * (1 - prevalence)
                ppv = sensitivity * prevalence / denom if denom else float("nan")
                npv_denom = specificity * (1 - prevalence) + (1 - sensitivity) * prevalence
                npv = specificity * (1 - prevalence) / npv_denom if npv_denom else float("nan")
                rows.append({
                    "negative_type": neg_type,
                    "model": model,
                    "threshold": 0.5,
                    "sensitivity": sensitivity,
                    "specificity": specificity,
                    "prevalence": prevalence,
                    "prevalence_label": label,
                    "ppv": ppv,
                    "npv": npv,
                    "number_needed_to_inspect": 1 / ppv if ppv else float("inf"),
                })
    return pd.DataFrame(rows)


def bootstrap_negative_set_difference(frame_a: pd.DataFrame, frame_b: pd.DataFrame,
                                     column: str, metric,
                                     n_boot: int = N_BOOTSTRAP,
                                     seed: int = SEED) -> dict:
    """Paired bootstrap of metric(A) - metric(B) across negative-set designs.

    Positives are shared by both datasets, so the same positive resample is
    used for both models; negatives are resampled independently within each
    design.
    """
    pos_a = (frame_a[frame_a["label"] == 1]
             .sort_values("window_id")[column].to_numpy())
    pos_ids = (frame_a[frame_a["label"] == 1]
               .sort_values("window_id")["window_id"].to_numpy())
    pos_b_frame = frame_b[frame_b["label"] == 1].sort_values("window_id")
    if not np.array_equal(pos_ids, pos_b_frame["window_id"].to_numpy()):
        raise ValueError("Positive windows are not aligned between datasets")
    pos_b = pos_b_frame[column].to_numpy()
    neg_a = frame_a.loc[frame_a["label"] == 0, column].to_numpy()
    neg_b = frame_b.loc[frame_b["label"] == 0, column].to_numpy()

    def score(pos, neg):
        y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
        s = np.concatenate([pos, neg])
        return metric(y, s)

    point = score(pos_a, neg_a) - score(pos_b, neg_b)
    rng = np.random.default_rng(seed)
    diffs = []
    for _ in range(n_boot):
        ip = rng.integers(0, len(pos_a), len(pos_a))
        ia = rng.integers(0, len(neg_a), len(neg_a))
        ib = rng.integers(0, len(neg_b), len(neg_b))
        diffs.append(score(pos_a[ip], neg_a[ia]) - score(pos_b[ip], neg_b[ib]))
    diffs = np.asarray(diffs)
    return {
        "difference": float(point),
        "ci_low": float(np.percentile(diffs, 2.5)),
        "ci_high": float(np.percentile(diffs, 97.5)),
        # Two-sided bootstrap "p-value": share of replicates on the other side of 0.
        "p_two_sided": float(2 * min((diffs <= 0).mean(), (diffs >= 0).mean())),
    }


def compare_negative_sets(frames: dict[str, pd.DataFrame],
                          n_boot: int = N_BOOTSTRAP, seed: int = SEED) -> dict:
    """Difference (random - distance_matched) for both models and metrics."""
    results = {}
    for model, column in (("logistic_regression", "lr_score"),
                          ("random_forest", "rf_score")):
        for metric_name, metric in (("roc_auc", roc_auc_score),
                                    ("pr_auc", average_precision_score)):
            results[f"{model}_{metric_name}"] = bootstrap_negative_set_difference(
                frames["random"], frames["distance_matched"], column, metric,
                n_boot=n_boot, seed=seed)
    return results
