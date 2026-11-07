"""Offline smoke test on the 45-window committed excerpt.

This checks that feature extraction, model fitting and plotting code paths run
without network access. It is NOT a result reproduction: 15 positives and 30
negatives are far too few for meaningful AUC estimates.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ars import config, features  # noqa: E402

SAMPLE = config.SAMPLE_DIR / "sample_dataset.csv"


def main() -> int:
    df = pd.read_csv(SAMPLE)
    print(f"Loaded {len(df)} windows from {SAMPLE.name}")

    # Recompute features from sequence to exercise the real code path.
    recomputed = np.array([
        [features.compute_features(seq, dist)[name]
         for name in config.FEATURES]
        for seq, dist in zip(df["sequence"], df["dist_bp"])
    ])
    stored = df[config.FEATURES].to_numpy(dtype=float)
    if not np.allclose(recomputed, stored, atol=1e-9):
        raise AssertionError("Recomputed features differ from the committed sample")

    X = stored
    y = df["label"].to_numpy()
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, stratify=y, random_state=config.SEED)
    model = Pipeline([
        ("scaler", StandardScaler()),
        ("lr", LogisticRegression(C=1.0, max_iter=1000)),
    ])
    model.fit(X_train, y_train)
    auc = roc_auc_score(y_test, model.predict_proba(X_test)[:, 1])
    print(f"Smoke-test held-out ROC-AUC: {auc:.3f} (n_test={len(y_test)})")
    print("SMOKE TEST ONLY - use scripts/run_all.py on the full data for results.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
