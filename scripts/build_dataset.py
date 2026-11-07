"""Build labelled 500 bp windows and extract the five sequence features."""

import argparse
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ars import config, dataset, evaluation, features, genome  # noqa: E402

COLUMNS = ["window_id", "chr", "start", "end", "mid", "label", "negative_type",
           "origin_name", "dist_bp", "sequence"] + config.FEATURES


def distance_match_quality(positives: list[dict], matched: list[dict]) -> dict:
    pos = pd.Series([p["dist_bp"] for p in positives])
    neg = pd.Series([n["dist_bp"] for n in matched])
    rel = np.abs(np.log1p(neg.to_numpy()) - np.log1p(pos.to_numpy()))
    return {
        "pos_median_dist_bp": float(pos.median()),
        "matched_median_dist_bp": float(neg.median()),
        "median_abs_log1p_difference": float(np.median(rel)),
        "p95_abs_log1p_difference": float(np.percentile(rel, 95)),
        "pos_fraction_within_10kb_of_landmark": float((pos < 10_000).mean()),
        "matched_fraction_within_10kb_of_landmark": float((neg < 10_000).mean()),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--window", type=int, default=config.WINDOW)
    parser.add_argument("--seed", type=int, default=config.SEED)
    args = parser.parse_args(argv)

    for path in (config.GENOME_FILE, config.ORIDB_FILE, config.SGD_FEATURES_FILE):
        if not path.exists():
            raise FileNotFoundError(
                f"Missing {path}. Run scripts/download_data.py first.")

    print("Loading reference genome and annotations")
    seqs = genome.read_fasta(config.GENOME_FILE)
    centromeres = genome.read_centromeres(config.SGD_FEATURES_FILE)
    origins = dataset.load_origins(config.ORIDB_FILE)
    print(f"  {len(seqs) - 1} nuclear chromosomes, "
          f"{len(origins)} confirmed OriDB origins, "
          f"{len(centromeres)} centromeres")

    positives = dataset.positive_windows(origins, seqs, centromeres, args.window)
    print(f"  positives after merging overlapping windows: {len(positives)}")

    randoms = dataset.random_negatives(
        seqs, positives, centromeres, args.window, args.seed, config.POS_EXCLUSION)
    matched = dataset.distance_matched_negatives(
        seqs, positives, centromeres, args.window, args.seed, config.POS_EXCLUSION)
    print(f"  random negatives: {len(randoms)}, "
          f"distance-matched negatives: {len(matched)}")

    records = features.add_features(positives + randoms + matched, args.window)
    df = pd.DataFrame(records)[COLUMNS]

    config.PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    out = config.PROCESSED_DIR / "dataset.csv"
    df.to_csv(out, index=False)
    print(f"Dataset written to {out} ({len(df)} windows)")

    quality = distance_match_quality(positives, matched)
    summary = {
        "generated_on": date.today().isoformat(),
        "window_bp": args.window,
        "seed": args.seed,
        "nuclear_genome_bp": int(sum(len(seqs[c]) for c in seqs if c.isdigit())),
        "dataset": evaluation.summarize_dataset(df),
        "distance_matching_quality": quality,
        "positive_definition": "OriDB status == Confirmed, midpoint of annotated element",
        "negative_definitions": {
            "random": "uniform genome-wide window, >=1000 bp from any positive",
            "distance_matched": "one per positive, same chromosome, matched "
                                "distance to nearest centromere/telomere",
        },
        "features": config.FEATURES,
    }
    summary_path = config.PROCESSED_DIR / "dataset_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Summary written to {summary_path}")

    # Small, committed excerpt for offline smoke tests (15 + 15 + 15 windows).
    sample = pd.concat([
        df[df["negative_type"] == "positive"].head(15),
        df[df["negative_type"] == "random"].head(15),
        df[df["negative_type"] == "distance_matched"].head(15),
    ])
    config.SAMPLE_DIR.mkdir(parents=True, exist_ok=True)
    sample.to_csv(config.SAMPLE_DIR / "sample_dataset.csv", index=False)
    print(f"Sample excerpt written to {config.SAMPLE_DIR / 'sample_dataset.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
