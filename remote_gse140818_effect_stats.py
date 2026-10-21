#!/usr/bin/env python3
"""Remote-only effect sizes and uncertainty for GSE140818 contrasts."""

import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


ROOT = pathlib.Path(sys.argv[1])
INPUT = ROOT / "results/public_mitotane_geo/GSE140818_HAC15_mito_scores.tsv"
OUT = ROOT / "results/public_mitotane_geo/GSE140818_HAC15_effect_stats.json"
SEED = 20261021
RNG = np.random.default_rng(SEED)


def cliff_delta(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    diff = a[:, None] - b[None, :]
    return float((np.sum(diff > 0) - np.sum(diff < 0)) / diff.size)


def contrast(frame, a_mask, b_mask, feature, label, n_boot=5000, n_perm=10000):
    a = frame.loc[a_mask, feature].dropna().to_numpy(float)
    b = frame.loc[b_mask, feature].dropna().to_numpy(float)
    if len(a) < 2 or len(b) < 2:
        return {"contrast": label, "feature": feature, "status": "insufficient"}
    observed = float(np.median(a) - np.median(b))
    delta = cliff_delta(a, b)
    boot = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        boot[i] = np.median(RNG.choice(a, len(a), replace=True)) - np.median(RNG.choice(b, len(b), replace=True))
    pooled = np.concatenate([a, b])
    n_a = len(a)
    null = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        perm = RNG.permutation(pooled)
        null[i] = np.median(perm[:n_a]) - np.median(perm[n_a:])
    p_perm = float((np.sum(np.abs(null) >= abs(observed)) + 1) / (n_perm + 1))
    mw = mannwhitneyu(a, b, alternative="two-sided")
    return {
        "contrast": label,
        "feature": feature,
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(np.median(a)),
        "median_b": float(np.median(b)),
        "median_difference_a_minus_b": observed,
        "median_difference_bootstrap_ci95": [float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))],
        "cliffs_delta_a_minus_b": delta,
        "mann_whitney_u": float(mw.statistic),
        "mann_whitney_p": float(mw.pvalue),
        "permutation_median_p": p_perm,
        "n_bootstrap": n_boot,
        "n_permutations": n_perm,
    }


frame = pd.read_csv(INPUT, sep="\t").set_index("Unnamed: 0")
frame["treatment"] = frame["treatment"].astype(str)
frame["clone"] = frame["clone"].astype(str)
features = ["OXPHOS", "ROS", "MITO_COMPOSITE", "G2M"]
contrasts = []
for feature in features:
    contrasts.append(contrast(frame, frame["clone"].eq("resistant") & frame["treatment"].eq("control"), frame["clone"].eq("nonresistant") & frame["treatment"].eq("control"), feature, "DMSO resistant minus nonresistant"))
    contrasts.append(contrast(frame, frame["clone"].eq("nonresistant") & frame["treatment"].eq("mitotane"), frame["clone"].eq("nonresistant") & frame["treatment"].eq("control"), feature, "nonresistant acute mitotane minus DMSO"))
    contrasts.append(contrast(frame, frame["clone"].eq("resistant") & frame["treatment"].eq("mitotane"), frame["clone"].eq("resistant") & frame["treatment"].eq("control"), feature, "resistant acute mitotane minus DMSO"))

result = {
    "cohort": "GSE140818",
    "seed": SEED,
    "design": "2 x 2 clone-by-treatment design; six arrays per cell of the factorial design",
    "contrasts": contrasts,
    "claim_boundary": "cell-model effect sizes and uncertainty; not patient IC50/AUC or clinical response",
}
OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
