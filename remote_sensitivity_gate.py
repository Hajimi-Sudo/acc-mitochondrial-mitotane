#!/usr/bin/env python
"""Remote-only sensitivity and evidence-gate checks for the ACC study."""
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler


ROOT = pathlib.Path(sys.argv[1])
FORMAL = ROOT / "results/formal"
OUT = ROOT / "results/gate"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261020
rng = np.random.default_rng(SEED)


def permutation_spearman(x, y, n=10000):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    observed = float(spearmanr(x, y).statistic)
    null = np.empty(n, dtype=float)
    for i in range(n):
        null[i] = spearmanr(x, rng.permutation(y)).statistic
    p = (np.sum(np.abs(null) >= abs(observed)) + 1) / (n + 1)
    return {"n": int(len(x)), "rho": observed, "permutation_p": float(p), "n_permutations": n}


def state_event_permutation(state, event, n=10000):
    state = np.asarray(state, dtype=bool)
    event = np.asarray(event, dtype=float)
    keep = np.isfinite(event)
    state, event = state[keep], event[keep]
    observed = float(event[state].mean() - event[~state].mean())
    null = np.empty(n, dtype=float)
    for i in range(n):
        perm = rng.permutation(state)
        null[i] = event[perm].mean() - event[~perm].mean()
    p = (np.sum(np.abs(null) >= abs(observed)) + 1) / (n + 1)
    return {"n": int(len(event)), "event_rate_difference_high_minus_low": observed, "permutation_p": float(p), "n_permutations": n}


def bootstrap_logistic(frame, outcome, features, n=1000):
    use = frame[[outcome] + features].dropna().copy()
    y = use[outcome].astype(int).to_numpy()
    x = use[features].to_numpy(float)
    if len(use) < 20 or len(np.unique(y)) != 2:
        return {"n": int(len(use)), "status": "insufficient"}
    estimates = []
    for _ in range(n):
        idx = rng.integers(0, len(use), len(use))
        if len(np.unique(y[idx])) != 2:
            continue
        xs = StandardScaler().fit_transform(x[idx])
        model = LogisticRegression(penalty=None, solver="lbfgs", max_iter=2000).fit(xs, y[idx])
        estimates.append(model.coef_[0])
    estimates = np.asarray(estimates)
    return {
        "n": int(len(use)),
        "features": features,
        "n_bootstrap": int(len(estimates)),
        "coefficient_median": dict(zip(features, np.median(estimates, axis=0).tolist())),
        "coefficient_ci95": {
            f: [float(np.quantile(estimates[:, i], 0.025)), float(np.quantile(estimates[:, i], 0.975))]
            for i, f in enumerate(features)
        },
    }


def gse_gate():
    data = pd.read_csv(FORMAL / "gse76019_formal_scores.tsv", sep="\t").set_index("Unnamed: 0")
    high = data["MITO_STATE"].eq("MITO_HIGH")
    event = pd.to_numeric(data["efs.event"], errors="coerce")
    efs = pd.to_numeric(data["efs.time"], errors="coerce")
    result = {
        "continuous_efs": permutation_spearman(data["MITO_COMPOSITE"], efs),
        "continuous_event": permutation_spearman(data["MITO_COMPOSITE"], event),
        "state_event": state_event_permutation(high, event),
        "state_efs_mannwhitney": {},
    }
    keep = high.notna() & efs.notna()
    a, b = efs[keep & high], efs[keep & ~high]
    u = mannwhitneyu(a, b, alternative="two-sided")
    result["state_efs_mannwhitney"] = {"n_high": int(len(a)), "n_low": int(len(b)), "U": float(u.statistic), "p": float(u.pvalue)}
    return result


def tcga_gate():
    scores = pd.read_csv(FORMAL / "tcga_acc_formal_scores.tsv.gz", sep="\t").set_index("Unnamed: 0")
    records = json.loads((ROOT / "data/raw/clinical/cbioportal_acc_tcga_patient_clinical_data.json").read_text())
    attrs = {}
    for record in records:
        attrs.setdefault(record["patientId"], {})[record["clinicalAttributeId"]] = record.get("value")
    manifest = pd.read_csv(ROOT / "data/processed/TCGA_ACC_file_sample_manifest.tsv", sep="\t").set_index("sample_id")
    clinical = pd.DataFrame([attrs.get(row["case_id"], {}) for _, row in manifest.iterrows()], index=manifest.index)
    data = scores.join(clinical, how="inner")
    raw = data.get("HISTORY_ADRENAL_HORMONE_EXCESS", pd.Series(index=data.index))
    functional = raw.map(lambda x: np.nan if pd.isna(x) else (0 if str(x).strip() in {"", "None"} else 1))
    keep = functional.notna()
    tab = pd.crosstab(data.loc[keep, "MITO_STATE"], functional[keep])
    result = {"n_functional_complete": int(keep.sum()), "state_functional_table": tab.to_dict()}
    if tab.shape == (2, 2):
        table = tab.reindex(index=["MITO_LOW", "MITO_HIGH"], columns=[0, 1], fill_value=0).to_numpy()
        odds, p = fisher_exact(table)
        result["state_functional_fisher"] = {"table": table.tolist(), "odds_ratio": float(odds), "p": float(p)}
    frame = pd.DataFrame({"functional": functional, "MITO_COMPOSITE": data["MITO_COMPOSITE"], "G2M_CONTROL": data["G2M_CONTROL"]})
    result["g2m_adjusted_functional_bootstrap"] = bootstrap_logistic(frame, "functional", ["MITO_COMPOSITE", "G2M_CONTROL"])
    result["claim_boundary"] = "TCGA hormone field is observational and incomplete; this is not a clinical response endpoint."
    return result


result = {
    "seed": SEED,
    "gse76019": gse_gate(),
    "tcga": tcga_gate(),
    "evidence_gate": {
        "direct_mitotane_response_verified": False,
        "gsa_human_accession": "HRA006596",
        "gsa_human_access": "controlled_access",
        "promotion_blocker": "Sample-level IC50/AUC labels and matched RNA-seq are not publicly retrievable without DAC authorization.",
        "current_claim_ceiling": "Associative mitochondrial-state analysis and hypothesis-generating pharmacogenomic rationale only.",
    },
}
(OUT / "sensitivity_gate_summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
