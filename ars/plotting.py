"""Result figures (matplotlib only, no seaborn)."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import precision_recall_curve, roc_curve

from .config import FEATURE_LABELS, FEATURES

POS_COLOR = "#2c7fb8"
NEG_COLOR = "#d95f0e"
LR_COLOR = "#1b7837"
RF_COLOR = "#762a83"
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
})


def plot_feature_distributions(df, tests, out_path):
    """Two rows (negative set) x five columns (features)."""
    fig, axes = plt.subplots(2, len(FEATURES), figsize=(16, 6.2))
    for row, neg_type in enumerate(("random", "distance_matched")):
        sub_tests = tests[tests["negative_type"] == neg_type].set_index("feature")
        for col, feature in enumerate(FEATURES):
            ax = axes[row, col]
            pos = df.loc[df["negative_type"] == "positive", feature]
            neg = df.loc[df["negative_type"] == neg_type, feature]
            bins = np.linspace(
                min(pos.min(), neg.min()), max(pos.max(), neg.max()), 30)
            ax.hist(pos, bins=bins, density=True, alpha=0.55, color=POS_COLOR,
                    label="positive", edgecolor="none")
            ax.hist(neg, bins=bins, density=True, alpha=0.55, color=NEG_COLOR,
                    label=neg_type, edgecolor="none")
            ax.axvline(pos.median(), color=POS_COLOR, ls="--", lw=1)
            ax.axvline(neg.median(), color=NEG_COLOR, ls="--", lw=1)
            row_data = sub_tests.loc[feature]
            ax.set_title(
                f"{FEATURE_LABELS[feature]}\n"
                f"p={row_data['p_fdr_bh']:.1e}, delta={row_data['cliffs_delta']:+.2f}",
                fontsize=8)
            if row == 0:
                ax.legend(frameon=False, fontsize=7, loc="upper right")
            ax.set_yticks([])
            if col == 0:
                ax.set_ylabel(neg_type.replace("_", " "), fontsize=10)
    fig.suptitle("Feature distributions: positives vs each negative set", y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_roc_pr(oof, out_path):
    """Rows: negative set; columns: ROC and PR curves."""
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 8.2))
    for row, neg_type in enumerate(("random", "distance_matched")):
        y_true, scores = oof[neg_type]["y"], oof[neg_type]["scores"]
        prevalence = float(np.mean(y_true))
        for name, color in (("logistic_regression", LR_COLOR),
                            ("random_forest", RF_COLOR)):
            score = scores[name]
            fpr, tpr, _ = roc_curve(y_true, score)
            precision, recall, _ = precision_recall_curve(y_true, score)
            label = name.replace("_", " ")
            axes[row, 0].plot(fpr, tpr, color=color, lw=1.5, label=label)
            axes[row, 1].plot(recall, precision, color=color, lw=1.5, label=label)
        axes[row, 0].plot([0, 1], [0, 1], color="grey", ls=":", lw=1)
        axes[row, 1].axhline(prevalence, color="grey", ls=":", lw=1,
                             label=f"baseline ({prevalence:.2f})")
        axes[row, 0].set_xlabel("False positive rate")
        axes[row, 0].set_ylabel("True positive rate")
        axes[row, 0].set_title(f"ROC - {neg_type.replace('_', ' ')} negatives")
        axes[row, 1].set_xlabel("Recall")
        axes[row, 1].set_ylabel("Precision")
        axes[row, 1].set_title(f"Precision-recall - {neg_type.replace('_', ' ')}")
        for ax in axes[row]:
            ax.legend(frameon=False, fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_importance(importance, out_path):
    """Rows: negative set; columns: LR coefficients and RF importances."""
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.0))
    for row, neg_type in enumerate(("random", "distance_matched")):
        lr = importance[neg_type]["lr"].sort_values()
        rf = importance[neg_type]["rf"].sort_values()
        ax = axes[row, 0]
        ax.barh([FEATURE_LABELS[f].replace("\n", " ") for f in lr.index],
                lr.values, color=[LR_COLOR if v >= 0 else "#b2182b" for v in lr.values])
        ax.axvline(0, color="black", lw=0.8)
        ax.set_title(f"Logistic regression coefficients ({neg_type.replace('_', ' ')})")
        ax.set_xlabel("log-odds per 1 SD")
        ax = axes[row, 1]
        ax.barh([FEATURE_LABELS[f].replace("\n", " ") for f in rf.index],
                rf.values, color=RF_COLOR)
        ax.set_title(f"Random forest importance ({neg_type.replace('_', ' ')})")
        ax.set_xlabel("Mean decrease in impurity")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_metrics_comparison(metrics_rows, out_path):
    """Grouped bars with bootstrap 95% CI for ROC-AUC and PR-AUC."""
    datasets = ["random", "distance_matched"]
    models = ["logistic_regression", "random_forest"]
    colors = {"logistic_regression": LR_COLOR, "random_forest": RF_COLOR}
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
    x = np.arange(len(datasets))
    width = 0.34
    for ax, metric in zip(axes, ("roc_auc", "pr_auc")):
        for offset, model in zip((-width / 2, width / 2), models):
            values, errors = [], []
            for dataset in datasets:
                row = next(r for r in metrics_rows
                           if r["negative_type"] == dataset and r["model"] == model)
                values.append(row[metric])
                low, high = row[f"{metric}_low"], row[f"{metric}_high"]
                errors.append([row[metric] - low, high - row[metric]])
            ax.bar(x + offset, values, width, color=colors[model],
                   label=model.replace("_", " "), yerr=np.asarray(errors).T,
                   capsize=3, error_kw={"lw": 0.8})
        ax.set_xticks(x)
        ax.set_xticklabels([d.replace("_", " ") for d in datasets])
        ax.set_ylim(0.45, 1.0)
        ax.set_ylabel(metric.replace("_", " ").upper())
        ax.set_title(f"{metric.replace('_', ' ').upper()} by negative-set design")
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
