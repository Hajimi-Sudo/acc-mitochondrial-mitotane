#!/usr/bin/env python3
"""Remote-only G2M-adjusted and factorial sensitivity analysis for GSE140818."""

import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


ROOT = pathlib.Path(sys.argv[1])
INPUT = ROOT / "results/public_mitotane_geo/GSE140818_HAC15_mito_scores.tsv"
OUT = ROOT / "results/public_mitotane_geo/GSE140818_HAC15_adjusted_stats.json"
SEED = 20261021
RNG = np.random.default_rng(SEED)


def residualize(frame, feature, control):
    use = frame[[feature, control]].dropna().copy()
    x = np.column_stack([np.ones(len(use)), use[control].to_numpy(float)])
    coef = np.linalg.lstsq(x, use[feature].to_numpy(float), rcond=None)[0]
    out = pd.Series(index=frame.index, dtype=float)
    out.loc[use.index] = use[feature] - x @ coef
    return out, {"intercept": float(coef[0]), "g2m_slope": float(coef[1]), "n": int(len(use))}


def robust_contrast(a, b, label, n_boot=5000, n_perm=10000):
    a = np.asarray(a.dropna(), dtype=float)
    b = np.asarray(b.dropna(), dtype=float)
    observed = float(np.median(a) - np.median(b))
    boot = np.empty(n_boot)
    for i in range(n_boot):
        boot[i] = np.median(RNG.choice(a, len(a), replace=True)) - np.median(RNG.choice(b, len(b), replace=True))
    pooled = np.concatenate([a, b])
    null = np.empty(n_perm)
    for i in range(n_perm):
        perm = RNG.permutation(pooled)
        null[i] = np.median(perm[:len(a)]) - np.median(perm[len(a):])
    return {
        "contrast": label,
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_difference_a_minus_b": observed,
        "bootstrap_ci95": [float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))],
        "permutation_p": float((np.sum(np.abs(null) >= abs(observed)) + 1) / (n_perm + 1)),
        "mann_whitney_p": float(mannwhitneyu(a, b, alternative="two-sided").pvalue),
    }


frame = pd.read_csv(INPUT, sep="\t").set_index("Unnamed: 0")
frame["treatment"] = frame["treatment"].astype(str)
frame["clone"] = frame["clone"].astype(str)
features = ["OXPHOS", "ROS", "MITO_COMPOSITE"]
adjusted = {}
for feature in features:
    residual, fit = residualize(frame, feature, "G2M")
    adjusted[feature] = {
        "residualization": fit,
        "dmso_resistant_minus_nonresistant": robust_contrast(
            residual[frame["clone"].eq("resistant") & frame["treatment"].eq("control")],
            residual[frame["clone"].eq("nonresistant") & frame["treatment"].eq("control")],
            "G2M-adjusted DMSO resistant minus nonresistant",
        ),
    }

# Difference-in-differences: (resistant mitotane - resistant DMSO) minus
# (nonresistant mitotane - nonresistant DMSO), evaluated on sample-level scores.
factorial = {}
for feature in features:
    residual, _ = residualize(frame, feature, "G2M")
    cells = {
        (clone, treatment): residual[frame["clone"].eq(clone) & frame["treatment"].eq(treatment)].dropna().to_numpy(float)
        for clone in ["resistant", "nonresistant"]
        for treatment in ["mitotane", "control"]
    }
    observed = (np.median(cells[("resistant", "mitotane")]) - np.median(cells[("resistant", "control")])) - (
        np.median(cells[("nonresistant", "mitotane")]) - np.median(cells[("nonresistant", "control")])
    )
    factorial[feature] = {
        "g2m_adjusted_difference_in_differences": float(observed),
        "cell_medians": {f"{k[0]}_{k[1]}": float(np.median(v)) for k, v in cells.items()},
        "interpretation": "descriptive factorial interaction; n=6 per cell, no clinical response endpoint",
    }

result = {
    "cohort": "GSE140818",
    "seed": SEED,
    "g2m_adjusted_contrasts": adjusted,
    "g2m_adjusted_factorial_interactions": factorial,
    "claim_boundary": "cell-model adjusted association and interaction context; not patient IC50/AUC or clinical response",
}
# Benjamini-Hochberg adjustment across the three mitochondrial features within
# each adjusted contrast family.
for family_key in ["dmso_resistant_minus_nonresistant"]:
    pvals = [adjusted[f][family_key]["permutation_p"] for f in features]
    order = np.argsort(pvals)
    q = np.empty(len(pvals), dtype=float)
    running = 1.0
    for rank in range(len(pvals) - 1, -1, -1):
        idx = order[rank]
        running = min(running, pvals[idx] * len(pvals) / (rank + 1))
        q[idx] = min(1.0, running)
    for feature, value in zip(features, q):
        adjusted[feature][family_key]["permutation_fdr_bh"] = float(value)
OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
