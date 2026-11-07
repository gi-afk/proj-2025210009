"""Sequence feature extraction (5 features, deliberately small)."""

import math
import re
from collections import Counter

import numpy as np

# ACS consensus: (A/T)TTTA(T/C)(A/G)TTT(A/T), 11 bp.
ACS_PATTERN = "WTTTAYRTTTW"

# Probability of each base under the degenerate consensus. Fixed positions use
# 0.98 for the consensus base and 0.02 total for the other three; degenerate
# positions split 0.49/0.49 between the two allowed bases.
_PWM: np.ndarray = np.zeros((11, 4), dtype=float)
_BASES = "ACGT"
for _i, _code in enumerate(ACS_PATTERN):
    if _code == "W":
        _PWM[_i, _BASES.index("A")] = 0.49
        _PWM[_i, _BASES.index("T")] = 0.49
    elif _code == "Y":
        _PWM[_i, _BASES.index("C")] = 0.49
        _PWM[_i, _BASES.index("T")] = 0.49
    elif _code == "R":
        _PWM[_i, _BASES.index("A")] = 0.49
        _PWM[_i, _BASES.index("G")] = 0.49
    else:
        _PWM[_i, _BASES.index(_code)] = 0.98
_PWM = _PWM + 0.01  # pseudocount against log(0)
_PWM_LOG2 = np.log2(_PWM / 0.25)

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
_KMER_ENTROPY_MAX = math.log2(64)


def _encode(seq: str) -> np.ndarray:
    table = np.full(256, -1, dtype=np.int8)
    for idx, base in enumerate(_BASES):
        table[ord(base)] = idx
    return table[np.array([ord(c) for c in seq], dtype=np.int32)]


def acs_pwm_max(seq: str) -> float:
    """Maximum ACS PWM log-odds score over both strands and all offsets."""
    best = -np.inf
    for strand in (seq, seq.translate(_COMPLEMENT)[::-1]):
        codes = _encode(strand)
        n = len(codes) - 11 + 1
        if n <= 0:
            return float("nan")
        windows = np.lib.stride_tricks.sliding_window_view(codes, 11)
        valid = (windows >= 0).all(axis=1)
        if valid.sum() == 0:
            return float("nan")
        scores = np.zeros(valid.sum(), dtype=float)
        sub = windows[valid]
        for j in range(11):
            scores += _PWM_LOG2[j][sub[:, j]]
        best = max(best, float(scores.max()))
    return best


def kmer3_entropy(seq: str) -> float:
    """Shannon entropy of 64 tri-nucleotides, normalised to [0, 1]."""
    if len(seq) < 3:
        return float("nan")
    counts = Counter(seq[i:i + 3] for i in range(len(seq) - 2))
    counts = {k: v for k, v in counts.items() if set(k) <= set("ACGT")}
    total = sum(counts.values())
    if total == 0:
        return float("nan")
    probs = np.array(list(counts.values()), dtype=float) / total
    entropy = float(-(probs * np.log2(probs)).sum())
    return entropy / _KMER_ENTROPY_MAX


def compute_features(seq: str, dist_bp: float) -> dict[str, float]:
    seq = seq.upper()
    length = len(seq)
    gc = (seq.count("G") + seq.count("C")) / length
    at = (seq.count("A") + seq.count("T")) / length
    return {
        "gc_content": gc,
        "at_content": at,
        "acs_pwm_max": acs_pwm_max(seq),
        "kmer3_entropy": kmer3_entropy(seq),
        "dist_cen_tel_log": math.log1p(max(dist_bp, 0.0)),
    }


def add_features(records: list[dict], window: int) -> list[dict]:
    """Add the 5 feature values to window records (records include `sequence`)."""
    enriched = []
    for rec in records:
        rec = dict(rec)
        rec.update(compute_features(rec["sequence"], rec["dist_bp"]))
        enriched.append(rec)
    return enriched
