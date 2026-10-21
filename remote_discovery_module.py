#!/usr/bin/env python
"""Remote-only discovery module for mitochondrial-state biology in TCGA-ACC."""
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


ROOT = pathlib.Path(sys.argv[1])
FORMAL = ROOT / "results/formal"
OUT = ROOT / "results/discovery"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261020
rng = np.random.default_rng(SEED)


def bh(p_values):
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    q = np.minimum.accumulate((ranked[::-1] * len(p) / np.arange(len(p), 0, -1))[::-1])
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def gene_list(*names):
    return sorted(set(names))


def gene_level_comparison(expr, state, genes):
    rows = []
    high = state.eq("MITO_HIGH")
    for gene in genes:
        if gene not in expr.index:
            continue
        a = expr.loc[gene, high].to_numpy(float)
        b = expr.loc[gene, ~high].to_numpy(float)
        if len(a) < 3 or len(b) < 3:
            continue
        test = mannwhitneyu(a, b, alternative="two-sided")
        # Rank-biserial effect: positive means higher in MITO_HIGH.
        delta = float(2 * test.statistic / (len(a) * len(b)) - 1)
        rows.append({"gene": gene, "n_high": len(a), "n_low": len(b), "high_median": np.median(a), "low_median": np.median(b), "rank_biserial": delta, "p": float(test.pvalue)})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr_bh"] = bh(out["p"].to_numpy())
        out = out.sort_values(["fdr_bh", "p", "gene"])
    return out


def bootstrap_state_stability(scores, n=1000):
    use = scores[["MITO_RESPIRATION", "MITO_OXSTRESS"]].dropna()
    x = StandardScaler().fit_transform(use)
    ref = KMeans(n_clusters=2, n_init=100, random_state=SEED).fit(x)
    ref_labels = ref.labels_
    ref_high = int(pd.Series(use["MITO_RESPIRATION"].to_numpy()).groupby(ref_labels).mean().idxmax())
    ref_high_labels = ref_labels == ref_high
    exact = 0
    overlaps = []
    for _ in range(n):
        idx = rng.integers(0, len(use), len(use))
        model = KMeans(n_clusters=2, n_init=20, random_state=int(rng.integers(0, 2**31 - 1))).fit(x[idx])
        high = int(pd.Series(use["MITO_RESPIRATION"].to_numpy()[idx]).groupby(model.labels_).mean().idxmax())
        pred = model.labels_ == high
        # Compare bootstrap assignments to reference assignments on sampled rows.
        ref_sample = ref_high_labels[idx]
        overlap = max(np.mean(pred == ref_sample), np.mean((~pred) == ref_sample))
        overlaps.append(overlap)
        exact += int(overlap >= 0.8)
    return {"n_bootstrap": n, "mean_assignment_overlap": float(np.mean(overlaps)), "ci95_assignment_overlap": [float(np.quantile(overlaps, 0.025)), float(np.quantile(overlaps, 0.975))], "fraction_overlap_ge_0.8": float(exact / n)}


def load_clinical():
    records = json.loads((ROOT / "data/raw/clinical/cbioportal_acc_tcga_patient_clinical_data.json").read_text())
    attrs = {}
    for record in records:
        attrs.setdefault(record["patientId"], {})[record["clinicalAttributeId"]] = record.get("value")
    manifest = pd.read_csv(ROOT / "data/processed/TCGA_ACC_file_sample_manifest.tsv", sep="\t").set_index("sample_id")
    return pd.DataFrame([attrs.get(row["case_id"], {}) for _, row in manifest.iterrows()], index=manifest.index)


counts = pd.read_csv(ROOT / "data/processed/TCGA_ACC_protein_coding_unstranded_counts.tsv.gz", sep="\t", compression="gzip").set_index("gene_name")
counts = counts.drop(columns=["gene_id"]).loc[lambda x: ~x.index.duplicated()]
expr = np.log2(counts.astype(float) + 1.0)
scores = pd.read_csv(FORMAL / "tcga_acc_formal_scores.tsv.gz", sep="\t").set_index("Unnamed: 0")
state = scores["MITO_STATE"]

mitochondrial = gene_list("NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6", "NDUFA7", "NDUFA8", "NDUFA9", "NDUFA10", "NDUFA11", "NDUFA12", "NDUFA13", "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10", "NDUFB11", "NDUFB12", "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8", "NDUFSH", "SDHA", "SDHB", "SDHC", "SDHD", "UQCRC1", "UQCRC2", "UQCRFS1", "UQCRB", "UQCRQ", "COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1", "COX7A2", "COX7C", "COX8A", "ATP5F1A", "ATP5F1B", "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5PD", "ATP5PO")
steroid = gene_list("STAR", "CYP11A1", "CYP11B1", "CYP11B2", "CYP17A1", "CYP21A2", "CYP19A1", "HSD3B1", "HSD3B2", "HSD17B1", "HSD17B3", "AKR1C3", "SCARB1", "MC2R", "NR5A1")
lipid = gene_list("CYP27A1", "ABCA1", "ABCG1", "SOAT1", "LIPE", "PNPLA2", "PLIN2", "PLIN3", "NR1H3", "SREBF1", "SREBF2", "LDLR", "LIPA")
cell_cycle = gene_list("MKI67", "PCNA", "TOP2A", "CCNB1", "CCNB2", "CDC20", "CDC25C", "CDK1", "AURKA", "AURKB", "BUB1", "BUB1B", "TPX2")

comparisons = {
    "mitochondrial_core": gene_level_comparison(expr, state, mitochondrial),
    "steroidogenesis": gene_level_comparison(expr, state, steroid),
    "lipid_handling": gene_level_comparison(expr, state, lipid),
    "cell_cycle_control": gene_level_comparison(expr, state, cell_cycle),
}
for name, table in comparisons.items():
    table.to_csv(OUT / f"tcga_{name}_mito_high_low.tsv", sep="\t", index=False)

clinical = load_clinical()
missingness = {column: {"n_present": int(clinical[column].notna().sum()), "n_missing": int(clinical[column].isna().sum()), "n_unique": int(clinical[column].dropna().astype(str).nunique())} for column in clinical.columns}

summary = {
    "cohort": "TCGA-ACC",
    "n_samples": int(len(expr.columns)),
    "state_counts": state.value_counts().to_dict(),
    "bootstrap_state_stability": bootstrap_state_stability(scores),
    "top_discovery_hits": {name: table.head(10).to_dict(orient="records") for name, table in comparisons.items()},
    "clinical_field_missingness": missingness,
    "interpretation": "Discovery associations only; no direct mitotane-response endpoint is present in TCGA-ACC.",
}
(OUT / "tcga_discovery_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
print(json.dumps(summary, indent=2, default=str))
