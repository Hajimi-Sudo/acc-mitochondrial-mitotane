#!/usr/bin/env python3
"""Remote-only factorial analysis of the GSE140818 HAC15 experiment."""
import json
import pathlib
import sys

import numpy as np
import pandas as pd


ROOT = pathlib.Path(sys.argv[1])
INPUT = ROOT / "results/public_mitotane_geo/GSE140818_HAC15_mito_scores.tsv"
OUT = ROOT / "results/public_mitotane_geo/GSE140818_HAC15_factorial_stats.json"
SEED = 20261023
RNG = np.random.default_rng(SEED)


def residualize(frame, feature, control):
    use = frame[[feature, control]].dropna()
    design = np.column_stack([np.ones(len(use)), use[control].to_numpy(float)])
    coef = np.linalg.lstsq(design, use[feature].to_numpy(float), rcond=None)[0]
    out = pd.Series(index=frame.index, dtype=float)
    out.loc[use.index] = use[feature].to_numpy(float) - design @ coef
    return out, {"n": int(len(use)), "intercept": float(coef[0]), "control_slope": float(coef[1])}


def fit_ols(frame, value_col):
    use = frame[[value_col, "resistance", "treatment"]].dropna().copy()
    # resistance=1 and treatment=1 are the resistant and mitotane groups.
    x = np.column_stack([
        np.ones(len(use)),
        use["resistance"].to_numpy(float),
        use["treatment"].to_numpy(float),
        use["resistance"].to_numpy(float) * use["treatment"].to_numpy(float),
    ])
    y = use[value_col].to_numpy(float)
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    fitted = x @ beta
    resid = y - fitted
    df = len(y) - x.shape[1]
    mse = float(np.sum(resid ** 2) / df) if df > 0 else np.nan
    cov = mse * np.linalg.pinv(x.T @ x)
    se = np.sqrt(np.maximum(np.diag(cov), 0))
    t = np.divide(beta, se, out=np.full_like(beta, np.nan), where=se > 0)
    return {
        "n": int(len(y)),
        "coefficients": {
            "intercept_nonresistant_control": float(beta[0]),
            "resistant_minus_nonresistant_control": float(beta[1]),
            "mitotane_minus_control_nonresistant": float(beta[2]),
            "resistance_by_treatment_interaction": float(beta[3]),
        },
        "standard_errors": {str(i): float(v) for i, v in enumerate(se)},
        "t_statistics": {str(i): float(v) for i, v in enumerate(t)},
        "residual_mse": mse,
    }


def permutation_interaction(frame, value_col, n=20000):
    use = frame[[value_col, "resistance", "treatment"]].dropna().copy()
    values = use[value_col].to_numpy(float)
    resistance = use["resistance"].to_numpy(int)
    treatment = use["treatment"].to_numpy(int)
    cells = {(r, t): values[(resistance == r) & (treatment == t)] for r in [0, 1] for t in [0, 1]}
    observed = cells[(1, 1)].mean() - cells[(1, 0)].mean() - cells[(0, 1)].mean() + cells[(0, 0)].mean()
    # Treatment is permuted within each clone-background stratum. Since the
    # design is balanced, the OLS interaction equals the four-cell mean contrast.
    null = np.empty(n, dtype=float)
    for i in range(n):
        means = {}
        for r in [0, 1]:
            vals = values[resistance == r]
            perm_treatment = RNG.permutation(treatment[resistance == r])
            means[(r, 0)] = vals[perm_treatment == 0].mean()
            means[(r, 1)] = vals[perm_treatment == 1].mean()
        null[i] = means[(1, 1)] - means[(1, 0)] - means[(0, 1)] + means[(0, 0)]
    p = (np.sum(np.abs(null) >= abs(observed)) + 1) / (n + 1)
    return {"observed_interaction": float(observed), "permutation_p": float(p), "n_permutations": int(n)}


def stratified_bootstrap(frame, value_col, n=5000):
    cells = {(r, t): frame[(frame["resistance"] == r) & (frame["treatment"] == t)][value_col].dropna().to_numpy(float) for r in [0, 1] for t in [0, 1]}
    boot_means = {(r, t): np.array([RNG.choice(cells[(r, t)], len(cells[(r, t)]), replace=True).mean() for _ in range(n)]) for r in [0, 1] for t in [0, 1]}
    estimates = np.column_stack([
        boot_means[(0, 0)],
        boot_means[(1, 0)] - boot_means[(0, 0)],
        boot_means[(0, 1)] - boot_means[(0, 0)],
        boot_means[(1, 1)] - boot_means[(1, 0)] - boot_means[(0, 1)] + boot_means[(0, 0)],
    ])
    keys = ["intercept_nonresistant_control", "resistant_minus_nonresistant_control", "mitotane_minus_control_nonresistant", "resistance_by_treatment_interaction"]
    return {key: {"median": float(np.median(estimates[:, i])), "ci95": [float(np.quantile(estimates[:, i], 0.025)), float(np.quantile(estimates[:, i], 0.975))]} for i, key in enumerate(keys)}


frame = pd.read_csv(INPUT, sep="\t").set_index("Unnamed: 0")
frame["resistance"] = frame["clone"].eq("resistant").astype(int)
frame["treatment"] = frame["treatment"].eq("mitotane").astype(int)
features = ["OXPHOS", "ROS", "MITO_COMPOSITE", "G2M"]
results = {}
for feature in features:
    raw = frame.rename(columns={feature: "value"}).copy()
    raw_result = fit_ols(raw, "value")
    raw_result["interaction_test"] = permutation_interaction(raw, "value")
    raw_result["bootstrap"] = stratified_bootstrap(raw, "value")
    if feature == "G2M":
        adj_result = {"status": "not_applicable", "reason": "G2M is the adjustment variable"}
    else:
        residual, residualization = residualize(frame, feature, "G2M")
        adj = frame.assign(value=residual)
        adj_result = fit_ols(adj, "value")
        adj_result["interaction_test"] = permutation_interaction(adj, "value")
        adj_result["bootstrap"] = stratified_bootstrap(adj, "value")
        adj_result["residualization"] = residualization
    results[feature] = {"raw": raw_result, "g2m_adjusted": adj_result}

payload = {
    "schema": "acc_mitochondrial_mitotane.gse140818_factorial.v1",
    "cohort": "GSE140818",
    "design": "balanced 2 x 2 clone-background by acute-treatment design; six arrays per cell",
    "features": results,
    "interpretation_boundary": "sample-level HAC15 factorial context; interaction is not a patient-level response biomarker",
    "seed": SEED,
}
OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
print(json.dumps(payload, indent=2))
