"""One-command pipeline: download -> build dataset -> train/evaluate -> figures.

Usage:
    python scripts/run_all.py                # full run
    python scripts/run_all.py --skip-download --skip-build
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ars import config, evaluation, models, plotting  # noqa: E402

import build_dataset  # noqa: E402
import download_data  # noqa: E402


def evaluate_subset(subset: pd.DataFrame, model_name: str,
                    feature_list: list[str] | None = None,
                    n_boot: int = config.N_BOOTSTRAP) -> dict:
    y, score, _ = models.logo_oof(subset, model_name, feature_list)
    point = evaluation.point_metrics(y, score)
    ci = evaluation.bootstrap_ci(y, score, n_boot=n_boot)
    row = {
        "model": model_name,
        "n_pos": point["n_pos"],
        "n_neg": point["n_neg"],
        "roc_auc": point["roc_auc"],
        "roc_auc_low": ci["roc_auc_ci"][0],
        "roc_auc_high": ci["roc_auc_ci"][1],
        "pr_auc": point["pr_auc"],
        "pr_auc_low": ci["pr_auc_ci"][0],
        "pr_auc_high": ci["pr_auc_ci"][1],
    }
    return row, y, score


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-download", action="store_true")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--bootstrap", type=int, default=config.N_BOOTSTRAP)
    args = parser.parse_args(argv)

    started = time.time()
    if not args.skip_download:
        download_data.main([])
    if not args.skip_build:
        build_dataset.main([])
    else:
        if not (config.PROCESSED_DIR / "dataset.csv").exists():
            raise FileNotFoundError("Processed dataset missing; run without --skip-build.")

    config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    config.FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(config.PROCESSED_DIR / "dataset.csv")
    summary_path = config.PROCESSED_DIR / "dataset_summary.json"
    dataset_summary = (json.loads(summary_path.read_text(encoding="utf-8"))
                       if summary_path.exists() else {})

    metrics_rows, oof_frames, oof_payload = [], [], {}
    importance = {}
    for neg_type in config.NEG_TYPES:
        subset = df[(df["negative_type"] == "positive")
                    | (df["negative_type"] == neg_type)].reset_index(drop=True)
        oof = {}
        for model_name in ("logistic_regression", "random_forest"):
            row, y, score = evaluate_subset(subset, model_name, n_boot=args.bootstrap)
            row.update({"negative_type": neg_type, "analysis": "main"})
            metrics_rows.append(row)
            oof[model_name] = score
            print(f"  {neg_type:16s} {model_name:20s} "
                  f"ROC-AUC={row['roc_auc']:.3f} PR-AUC={row['pr_auc']:.3f}")
        oof_payload[neg_type] = {"y": y, "scores": oof}
        oof_frames.append(pd.DataFrame({
            "window_id": subset["window_id"],
            "chr": subset["chr"],
            "label": y,
            "negative_type": neg_type,
            # Round to remove last-bit floating point noise from parallel RF.
            "lr_score": np.round(oof["logistic_regression"], 12),
            "rf_score": np.round(oof["random_forest"], 12),
        }))

        # Ablations isolate the source of predictive signal.
        for label, feature_list in (
            ("lr_distance_only", ["dist_cen_tel_log"]),
            ("lr_sequence_only", [f for f in config.FEATURES
                                  if f != "dist_cen_tel_log"]),
        ):
            row, _, _ = evaluate_subset(
                subset, "logistic_regression", feature_list, n_boot=args.bootstrap)
            row.update({"negative_type": neg_type, "analysis": "ablation",
                        "model": label})
            metrics_rows.append(row)
            print(f"  ablation {neg_type:16s} {label:20s} "
                  f"ROC-AUC={row['roc_auc']:.3f} PR-AUC={row['pr_auc']:.3f}")

        importance[neg_type] = {
            "lr": models.lr_coefficients(models.fit_full(subset, "logistic_regression")),
            "rf": models.rf_importances(models.fit_full(subset, "random_forest")),
        }

    metrics_df = pd.DataFrame(metrics_rows)
    metrics_df.to_csv(config.RESULTS_DIR / "metrics_table.csv", index=False,
                      lineterminator="\n")

    tests = pd.concat([evaluation.feature_tests(df, t) for t in config.NEG_TYPES])
    tests.to_csv(config.RESULTS_DIR / "feature_tests.csv", index=False,
                 lineterminator="\n")

    for neg_type in config.NEG_TYPES:
        pd.DataFrame({
            "feature": config.FEATURES,
            "lr_coefficient": importance[neg_type]["lr"].values,
            "rf_importance": importance[neg_type]["rf"].values,
        }).round(10).to_csv(config.RESULTS_DIR / f"importance_{neg_type}.csv",
                            index=False, lineterminator="\n")

    pd.concat(oof_frames).to_csv(config.RESULTS_DIR / "oof_predictions.csv",
                                 index=False, lineterminator="\n")

    frames_by_type = {f["negative_type"].iloc[0]: f for f in oof_frames}
    comparison = evaluation.compare_negative_sets(
        frames_by_type, n_boot=args.bootstrap)
    with open(config.RESULTS_DIR / "negative_set_comparison.json", "w",
              encoding="utf-8", newline="\n") as handle:
        json.dump(comparison, handle, indent=2)

    n_pos = int((df["label"] == 1).sum())
    genome_bp = dataset_summary.get("nuclear_genome_bp", 12_071_326)
    ppv = evaluation.ppv_table(oof_frames, genome_bp, n_pos, config.WINDOW)
    ppv.round(10).to_csv(config.RESULTS_DIR / "ppv_table.csv", index=False,
                        lineterminator="\n")

    plotting.plot_feature_distributions(
        df, tests, config.FIGURE_DIR / "fig1_feature_distributions.png")
    plotting.plot_roc_pr(oof_payload, config.FIGURE_DIR / "fig2_roc_pr.png")
    plotting.plot_importance(importance, config.FIGURE_DIR / "fig3_feature_importance.png")
    plotting.plot_metrics_comparison(
        [r for r in metrics_rows if r["analysis"] == "main"],
        config.FIGURE_DIR / "fig4_metrics_comparison.png")

    metadata = {
        "dataset": dataset_summary,
        "metrics": metrics_rows,
        "settings": {
            "window_bp": config.WINDOW,
            "seed": config.SEED,
            "cv": "LeaveOneGroupOut over chromosomes (16 folds)",
            "bootstrap_replicates": args.bootstrap,
            "model_defaults": {
                "logistic_regression": {"C": config.LR_C, "penalty": "l2"},
                "random_forest": {
                    "n_estimators": config.RF_N_ESTIMATORS,
                    "min_samples_leaf": config.RF_MIN_SAMPLES_LEAF,
                    "max_features": "sqrt",
                },
            },
        },
    }
    with open(config.RESULTS_DIR / "metrics.json", "w", encoding="utf-8",
              newline="\n") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)

    print("\nMain results")
    main_rows = metrics_df[metrics_df["analysis"] == "main"]
    print(main_rows[["negative_type", "model", "roc_auc", "pr_auc"]]
          .to_string(index=False))
    runtime = time.time() - started
    print(f"\nDone in {runtime:.1f} s. Figures: {config.FIGURE_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
