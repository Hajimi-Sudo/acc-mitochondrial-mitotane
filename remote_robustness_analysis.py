#!/usr/bin/env python3
"""Remote-only robustness analysis for the TCGA-ACC mitochondrial state.

This analysis tests whether the observed state and endocrine association depend
on the exact two-score construction. It writes only aggregate summaries.
"""
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


ROOT = pathlib.Path(sys.argv[1])
OUT = ROOT / "results" / "robustness"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261022


def read_gene_set(path):
    genes = []
    for line in pathlib.Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("VERSION"):
            continue
        gene = line.split("\t")[-1].strip().upper()
        if gene:
            genes.append(gene)
    return sorted(set(genes))


def score(expr, genes):
    present = sorted(set(genes) & set(expr.index))
    block = expr.loc[present].astype(float)
    z = block.sub(block.mean(axis=1), axis=0).div(
        block.std(axis=1, ddof=1).replace(0, np.nan), axis=0
    )
    return z.mean(axis=0), len(present)


def state_from_frame(frame, features, seed=SEED):
    use = frame[features].dropna()
    x = StandardScaler().fit_transform(use.to_numpy(float))
    labels = KMeans(n_clusters=2, n_init=100, random_state=seed).fit_predict(x)
    means = use.assign(label=labels).groupby("label")[features[0]].mean()
    high = int(means.idxmax())
    out = pd.Series(index=frame.index, dtype=object)
    out.loc[use.index] = np.where(labels == high, "MITO_HIGH", "MITO_LOW")
    return out


def overlap(a, b):
    use = a.notna() & b.notna()
    return float(np.mean(a[use].to_numpy() == b[use].to_numpy()))


def rank_residual(x, y, control):
    frame = pd.DataFrame({"x": x, "y": y, "control": control}).dropna()
    if len(frame) < 10:
        return {"n": int(len(frame)), "rho": None, "p": None}
    rx = frame["x"].rank().to_numpy(float)
    ry = frame["y"].rank().to_numpy(float)
    rc = frame["control"].rank().to_numpy(float)
    design = np.column_stack([np.ones(len(frame)), rc])
    xr = rx - design @ np.linalg.lstsq(design, rx, rcond=None)[0]
    yr = ry - design @ np.linalg.lstsq(design, ry, rcond=None)[0]
    rho, p = spearmanr(xr, yr)
    return {"n": int(len(frame)), "rho": float(rho), "p": float(p)}


def category_summary(state, field):
    frame = pd.DataFrame({"state": state, "field": field}).dropna()
    frame["functional"] = ~frame["field"].astype(str).str.strip().str.upper().isin({"", "NONE"})
    table = pd.crosstab(frame["state"], frame["functional"])
    result = {
        "n_complete": int(len(frame)),
        "categories": frame["field"].astype(str).value_counts().sort_index().to_dict(),
        "functional_table": table.to_dict(),
    }
    if table.shape == (2, 2):
        ordered = table.reindex(index=["MITO_LOW", "MITO_HIGH"], columns=[False, True], fill_value=0).to_numpy()
        odds, p = fisher_exact(ordered)
        result["functional_fisher"] = {"table": ordered.tolist(), "odds_ratio": float(odds), "p": float(p)}
    return result


counts = pd.read_csv(
    ROOT / "data/processed/TCGA_ACC_protein_coding_unstranded_counts.tsv.gz",
    sep="\t", compression="gzip",
).set_index("gene_name")
counts = counts.drop(columns=["gene_id"]).loc[lambda x: ~x.index.duplicated()]
expr = np.log2(counts.astype(float) + 1.0)
sets = {
    "OXPHOS": read_gene_set(ROOT / "data/raw/gene_sets/hallmark_oxphos.dat"),
    "ROS": read_gene_set(ROOT / "data/raw/gene_sets/hallmark_ros.dat"),
    "G2M": read_gene_set(ROOT / "data/raw/gene_sets/hallmark_g2m.dat"),
}
scores = pd.DataFrame(index=expr.columns)
overlaps = {}
for name in sets:
    scores[name], overlaps[name] = score(expr, sets[name])
scores["MITO_COMPOSITE"] = scores[["OXPHOS", "ROS"]].mean(axis=1)
scores["G2M_CONTROL"] = scores["G2M"]
reference = pd.read_csv(ROOT / "results/formal/tcga_acc_formal_scores.tsv.gz", sep="\t").set_index("Unnamed: 0")["MITO_STATE"]

variants = {
    "two_score_kmeans": state_from_frame(scores, ["OXPHOS", "ROS"]),
    "oxphos_only_kmeans": state_from_frame(scores, ["OXPHOS"]),
    "ros_only_kmeans": state_from_frame(scores, ["ROS"]),
    "composite_median_split": pd.Series(
        np.where(scores["MITO_COMPOSITE"] >= scores["MITO_COMPOSITE"].median(), "MITO_HIGH", "MITO_LOW"),
        index=scores.index,
    ),
}

variant_rows = []
for name, state in variants.items():
    variant_rows.append({
        "variant": name,
        "n_high": int((state == "MITO_HIGH").sum()),
        "n_low": int((state == "MITO_LOW").sum()),
        "overlap_with_reference": overlap(state, reference),
    })

leave_one_rows = []
for set_name in ["OXPHOS", "ROS"]:
    genes = sorted(set(sets[set_name]) & set(expr.index))
    baseline = scores[set_name]
    for gene in genes:
        reduced = [g for g in genes if g != gene]
        if len(reduced) < 3:
            continue
        value, _ = score(expr, reduced)
        leave_one_rows.append({
            "gene_set": set_name,
            "omitted_gene": gene,
            "score_spearman_with_baseline": float(spearmanr(value, baseline).statistic),
        })
leave_one = pd.DataFrame(leave_one_rows)

records = json.loads((ROOT / "data/raw/clinical/cbioportal_acc_tcga_patient_clinical_data.json").read_text())
attrs = {}
for record in records:
    attrs.setdefault(record["patientId"], {})[record["clinicalAttributeId"]] = record.get("value")
manifest = pd.read_csv(ROOT / "data/processed/TCGA_ACC_file_sample_manifest.tsv", sep="\t").set_index("sample_id")
clinical = pd.DataFrame([attrs.get(row["case_id"], {}) for _, row in manifest.iterrows()], index=manifest.index)
field = clinical.get("HISTORY_ADRENAL_HORMONE_EXCESS", pd.Series(index=scores.index))

association_rows = []
for name, state in variants.items():
    association = category_summary(state, field)
    association_rows.append({"variant": name, **association})

result = {
    "schema": "acc_mitochondrial_mitotane.tcga_robustness.v1",
    "cohort": "TCGA-ACC",
    "n_samples": int(len(scores)),
    "gene_set_overlap": overlaps,
    "state_variants": variant_rows,
    "functional_hormone_by_variant": association_rows,
    "continuous_partial_spearman": {
        "OXPHOS_vs_functional_proxy_after_G2M": rank_residual(scores["OXPHOS"], field.ne("None").where(field.notna()), scores["G2M_CONTROL"]),
        "ROS_vs_functional_proxy_after_G2M": rank_residual(scores["ROS"], field.ne("None").where(field.notna()), scores["G2M_CONTROL"]),
    },
    "leave_one_gene_out": {
        "OXPHOS_min_rho": float(leave_one.loc[leave_one["gene_set"].eq("OXPHOS"), "score_spearman_with_baseline"].min()),
        "OXPHOS_median_rho": float(leave_one.loc[leave_one["gene_set"].eq("OXPHOS"), "score_spearman_with_baseline"].median()),
        "ROS_min_rho": float(leave_one.loc[leave_one["gene_set"].eq("ROS"), "score_spearman_with_baseline"].min()),
        "ROS_median_rho": float(leave_one.loc[leave_one["gene_set"].eq("ROS"), "score_spearman_with_baseline"].median()),
    },
    "claim_boundary": "robustness of an expression-state construction; no direct patient mitotane response endpoint",
}
(OUT / "tcga_robustness_summary.json").write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
pd.DataFrame(variant_rows).to_csv(OUT / "state_construction_sensitivity.tsv", sep="\t", index=False)
leave_one.to_csv(OUT / "leave_one_gene_out_scores.tsv", sep="\t", index=False)
print(json.dumps(result, indent=2, default=str))
