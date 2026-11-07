"""Unit tests that do not require network access."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ars import config, features  # noqa: E402


def test_sample_dataset_schema():
    df = pd.read_csv(config.SAMPLE_DIR / "sample_dataset.csv")
    assert len(df) == 45
    assert set(df["negative_type"]) == {"positive", "random", "distance_matched"}
    assert (df["sequence"].str.len() == config.WINDOW).all()
    for name in config.FEATURES:
        assert df[name].notna().all(), name
    assert df["gc_content"].between(0, 1).all()
    assert df["at_content"].between(0, 1).all()
    assert df["kmer3_entropy"].between(0, 1).all()
    # GC + AT must be exactly complementary by construction.
    assert np.allclose(df["gc_content"] + df["at_content"], 1.0)


def test_acs_pwm_prefers_consensus():
    consensus = "ATTTACATTTA"  # WTTTAYRTTTW
    shuffled = "GGGGCCCCGGG"
    assert features.acs_pwm_max(consensus) > features.acs_pwm_max(shuffled)
    assert features.acs_pwm_max(consensus) > 12


def test_entropy_range():
    seq = "ACGT" * 125
    value = features.kmer3_entropy(seq)
    assert 0 <= value <= 1
    assert features.kmer3_entropy("AAAA" * 125) < value


def test_feature_keys():
    values = features.compute_features("ACGT" * 125, 1000)
    assert set(values) == set(config.FEATURES)
    assert all(np.isfinite(v) for v in values.values())
