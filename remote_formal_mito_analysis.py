#!/usr/bin/env python
"""Remote-only formal mitochondrial/OXPHOS analysis for TCGA-ACC and GSE76019."""
import csv
import gzip
import json
import pathlib
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, kruskal, mannwhitneyu
from sklearn.linear_model import LogisticRegression
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


ROOT = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path(
    "/root/autodl-tmp/acc-mitochondrial-mitotane"
)
OUT = ROOT / "results" / "formal"
OUT.mkdir(parents=True, exist_ok=True)


def read_gene_set(path):
    genes = []
    for line in pathlib.Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("VERSION"):
            continue
        parts = line.split("\t")
        gene = parts[-1].strip().upper()
        if gene and re.fullmatch(r"[A-Z0-9-]+", gene):
            genes.append(gene)
    return sorted(set(genes))


SETS = {
    "OXPHOS": read_gene_set(ROOT / "data/raw/gene_sets/hallmark_oxphos.dat"),
    "ROS": read_gene_set(ROOT / "data/raw/gene_sets/hallmark_ros.dat"),
    "G2M": read_gene_set(ROOT / "data/raw/gene_sets/hallmark_g2m.dat"),
}


def rank_z(values):
    ranks = pd.Series(values).rank(method="average").to_numpy(dtype=float)
    return (ranks - np.nanmean(ranks)) / np.nanstd(ranks, ddof=1)


def score_expression(expr, sets):
    """Score genes with within-gene z normalization; rows are genes."""
    scores = pd.DataFrame(index=expr.columns, dtype=float)
    overlap = {}
    for name, requested in sets.items():
        genes = sorted(set(requested) & set(expr.index))
        overlap[name] = {"requested": len(requested), "present": len(genes), "genes": genes}
        if len(genes) < 3:
            scores[name] = np.nan
            continue
        block = expr.loc[genes].astype(float)
        normalized = block.sub(block.mean(axis=1), axis=0).div(
            block.std(axis=1, ddof=1).replace(0, np.nan), axis=0
        )
        scores[name] = normalized.mean(axis=0)
    scores["MITO_COMPOSITE"] = scores[["OXPHOS", "ROS"]].mean(axis=1)
    scores["MITO_RESPIRATION"] = scores["OXPHOS"]
    scores["MITO_OXSTRESS"] = scores["ROS"]
    scores["G2M_CONTROL"] = scores["G2M"]
    return scores, overlap


def residual_spearman(x, y, controls):
    frame = pd.DataFrame({"x": x, "y": y})
    if controls:
        frame = frame.join(pd.DataFrame(controls), how="inner")
    frame = frame.dropna()
    if len(frame) < 5:
        return {"n": int(len(frame)), "rho": None, "p": None}
    design = np.column_stack([np.ones(len(frame)), frame[list(controls)].to_numpy(float)])
    xr = frame["x"].to_numpy(float) - design @ np.linalg.lstsq(design, frame["x"], rcond=None)[0]
    yr = frame["y"].to_numpy(float) - design @ np.linalg.lstsq(design, frame["y"], rcond=None)[0]
    rho, p = spearmanr(xr, yr)
    return {"n": int(len(frame)), "rho": float(rho), "p": float(p)}


def state_labels(scores):
    use = scores[["MITO_RESPIRATION", "MITO_OXSTRESS"]].dropna()
    labels = KMeans(n_clusters=2, n_init=100, random_state=20261020).fit_predict(
        StandardScaler().fit_transform(use)
    )
    means = use.assign(label=labels).groupby("label")["MITO_RESPIRATION"].mean()
    high = int(means.idxmax())
    out = pd.Series(index=scores.index, dtype=object)
    out.loc[use.index] = np.where(labels == high, "MITO_HIGH", "MITO_LOW")
    return out, {"MITO_HIGH": int((out == "MITO_HIGH").sum()), "MITO_LOW": int((out == "MITO_LOW").sum())}


def adjusted_logistic(frame, outcome, features):
    """Small-sample standardized logistic model; missing outcomes are excluded."""
    use = frame[[outcome] + features].copy()
    use[outcome] = use[outcome].astype(str).str.upper().map({"YES": 1, "TRUE": 1, "NO": 0, "FALSE": 0})
    use = use.dropna()
    if len(use) < 10 or use[outcome].nunique() != 2:
        return {"n": int(len(use)), "status": "insufficient_binary_data"}
    x = use[features].to_numpy(float)
    x = StandardScaler().fit_transform(x)
    y = use[outcome].to_numpy(int)
    model = LogisticRegression(penalty=None, solver="lbfgs", max_iter=2000).fit(x, y)
    probs = model.predict_proba(x)[:, 1]
    # Training AUC is descriptive only; the primary quantity is the adjusted coefficient.
    ranks = pd.Series(probs).rank().to_numpy()
    auc = float((ranks[y == 1].mean() - (len(y[y == 1]) + 1) / 2) / len(y[y == 0]))
    return {
        "n": int(len(use)),
        "positive": int(y.sum()),
        "features": features,
        "standardized_coefficients": {f: float(c) for f, c in zip(features, model.coef_[0])},
        "odds_ratio_per_sd": {f: float(np.exp(c)) for f, c in zip(features, model.coef_[0])},
        "descriptive_training_auc": auc,
    }


def load_tcga():
    path = ROOT / "data/processed/TCGA_ACC_protein_coding_unstranded_counts.tsv.gz"
    counts = pd.read_csv(path, sep="\t", compression="gzip").set_index("gene_name")
    counts = counts.drop(columns=["gene_id"]).loc[lambda x: ~x.index.duplicated()]
    expr = np.log2(counts.astype(float) + 1.0)
    scores, overlap = score_expression(expr, SETS)
    labels, counts_state = state_labels(scores)
    scores["MITO_STATE"] = labels
    scores.to_csv(OUT / "tcga_acc_formal_scores.tsv.gz", sep="\t", compression="gzip")
    # Patient-level clinical data are JSON records keyed by clinicalAttributeId.
    records = json.loads((ROOT / "data/raw/clinical/cbioportal_acc_tcga_patient_clinical_data.json").read_text())
    attrs = {}
    for record in records:
        attrs.setdefault(record["patientId"], {})[record["clinicalAttributeId"]] = record.get("value")
    manifest = pd.read_csv(ROOT / "data/processed/TCGA_ACC_file_sample_manifest.tsv", sep="\t").set_index("sample_id")
    clinical_rows = []
    for sample_id, row in manifest.iterrows():
        d = {"sample_id": sample_id, "patient_id": row["case_id"]}
        d.update(attrs.get(row["case_id"], {}))
        clinical_rows.append(d)
    clinical_df = pd.DataFrame(clinical_rows).set_index("sample_id")
    joined = scores.join(clinical_df, how="left")
    tests = {}
    for field in ["HISTORY_ADRENAL_HORMONE_EXCESS", "PHARM_TX_MITOTANE_INDICATOR", "PHARM_TX_MITOTANE_THERAPUTIC_LEVELS"]:
        if field not in joined:
            continue
        x = joined[["MITO_COMPOSITE", field]].dropna()
        if len(x) >= 5:
            high = x["MITO_COMPOSITE"] >= x["MITO_COMPOSITE"].median()
            tab = pd.crosstab(high, x[field].astype(str))
            tests[field] = {"n": int(len(x)), "table": tab.to_dict()}
            if tab.shape == (2, 2):
                a, b = tab.iloc[0, 0], tab.iloc[0, 1]
                c, d = tab.iloc[1, 0], tab.iloc[1, 1]
                tests[field]["odds_ratio"] = float((a * d + 0.5) / (b * c + 0.5))
    tests["mito_vs_g2m"] = residual_spearman(joined["MITO_COMPOSITE"], joined["G2M_CONTROL"], {})
    tests["respiration_vs_g2m"] = residual_spearman(joined["MITO_RESPIRATION"], joined["G2M_CONTROL"], {})
    def encode_functional(value):
        if pd.isna(value):
            return np.nan
        return "NO" if str(value).strip().upper() in {"NONE", ""} else "YES"

    joined["FUNCTIONAL_HORMONE"] = joined.get("HISTORY_ADRENAL_HORMONE_EXCESS", pd.Series(index=joined.index)).map(encode_functional)
    adjusted = {
        "functional_hormone_mito_only": adjusted_logistic(joined, "FUNCTIONAL_HORMONE", ["MITO_COMPOSITE"]),
        "functional_hormone_mito_plus_g2m": adjusted_logistic(joined, "FUNCTIONAL_HORMONE", ["MITO_COMPOSITE", "G2M_CONTROL"]),
        "mitotane_exposure_mito_plus_g2m": adjusted_logistic(joined, "PHARM_TX_MITOTANE_INDICATOR", ["MITO_COMPOSITE", "G2M_CONTROL"]),
        "therapeutic_level_mito_plus_g2m": adjusted_logistic(joined, "PHARM_TX_MITOTANE_THERAPUTIC_LEVELS", ["MITO_COMPOSITE", "G2M_CONTROL"]),
    }
    result = {
        "cohort": "TCGA-ACC",
        "n_samples": int(len(scores)),
        "gene_set_overlap": {k: {"requested": v["requested"], "present": v["present"]} for k, v in overlap.items()},
        "state_counts": counts_state,
        "tests": tests,
        "g2m_adjusted_models": adjusted,
        "claim_boundary": "TCGA treatment fields are exposure/therapeutic-level proxies, not direct mitotane response labels",
    }
    return result


def parse_gse_matrix():
    path = ROOT / "data/raw/gse76019/GSE76019_series_matrix.txt.gz"
    with gzip.open(path, "rt", errors="replace") as handle:
        lines = handle.readlines()
    begin = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_begin"))
    end = next(i for i, line in enumerate(lines[begin + 1:], begin + 1) if line.startswith("!series_matrix_table_end"))
    header = lines[begin + 1].rstrip("\n").split("\t")
    rows = []
    for line in lines[begin + 2:end]:
        fields = line.rstrip("\n").split("\t")
        if len(fields) == len(header):
            rows.append(fields)
    data = pd.DataFrame(rows, columns=header).set_index(header[0])
    data.index = data.index.str.strip('"')
    data = data.apply(lambda col: pd.to_numeric(col.str.strip('"'), errors="coerce"))
    samples = [x.strip('"') for x in header[1:]]
    data.columns = samples
    metadata = {}
    for line in lines[:begin]:
        if line.startswith("!Sample_geo_accession") or line.startswith("!Sample_characteristics_ch1"):
            parts = line.rstrip("\n").split("\t")
            metadata.setdefault(parts[0], []).append([x.strip('"') for x in parts[1:]])
    return data, metadata


def load_gse():
    matrix, metadata = parse_gse_matrix()
    ann_path = ROOT / "data/raw/gse76019/GPL13158_full.txt"
    ann = pd.read_csv(ann_path, sep="\t", skiprows=15085, dtype=str, low_memory=False)
    ann = ann.rename(columns={ann.columns[0]: "probe"})
    ann = ann[["probe", "Gene Symbol"]].fillna("")
    ann["gene"] = ann["Gene Symbol"].str.split("///").str[0].str.strip().str.upper()
    ann = ann[ann["gene"].ne("")].drop_duplicates("probe")
    matrix.index = matrix.index.astype(str).str.strip('"')
    common = matrix.index.intersection(ann["probe"])
    mapped = ann.set_index("probe").loc[common]
    expr_probe = matrix.loc[common]
    expr = expr_probe.assign(gene=mapped["gene"].to_numpy()).groupby("gene").median()
    scores, overlap = score_expression(expr, SETS)
    labels, counts_state = state_labels(scores)
    scores["MITO_STATE"] = labels
    accession = metadata.get("!Sample_geo_accession", [[]])[0]
    chars = metadata.get("!Sample_characteristics_ch1", [])
    clinical = pd.DataFrame(index=accession)
    for row in chars:
        for i, value in enumerate(row):
            if ":" in value:
                key, val = value.split(":", 1)
                clinical.loc[accession[i], key.strip()] = val.strip()
    clinical["efs.time"] = pd.to_numeric(clinical.get("efs.time"), errors="coerce")
    clinical["efs.event"] = pd.to_numeric(clinical.get("efs.event"), errors="coerce")
    joined = scores.join(clinical, how="left")
    tests = {}
    for field in ["efs.time", "efs.event"]:
        if field in joined:
            x = pd.to_numeric(joined[field], errors="coerce")
            keep = x.notna() & joined["MITO_COMPOSITE"].notna()
            if keep.sum() >= 5:
                rho, p = spearmanr(joined.loc[keep, "MITO_COMPOSITE"], x[keep])
                tests[field] = {"n": int(keep.sum()), "rho": float(rho), "p": float(p)}
    if "Stage" in joined:
        groups = [joined.loc[joined["Stage"] == stage, "MITO_COMPOSITE"].dropna().to_numpy() for stage in sorted(joined["Stage"].dropna().unique())]
        groups = [x for x in groups if len(x) >= 2]
        if len(groups) >= 2:
            h, p = kruskal(*groups)
            tests["stage_kruskal"] = {"n_groups": len(groups), "H": float(h), "p": float(p)}
    tests["mito_vs_g2m"] = residual_spearman(joined["MITO_COMPOSITE"], joined["G2M_CONTROL"], {})
    joined.to_csv(OUT / "gse76019_formal_scores.tsv", sep="\t")
    result = {
        "cohort": "GSE76019",
        "n_samples": int(len(scores)),
        "n_expression_rows": int(matrix.shape[0]),
        "n_genes_after_probe_collapse": int(expr.shape[0]),
        "gene_set_overlap": {k: {"requested": v["requested"], "present": v["present"]} for k, v in overlap.items()},
        "state_counts": counts_state,
        "tests": tests,
        "cohort_role": "pediatric expression/outcome context only; no direct mitotane response labels",
    }
    return result


results = {"tcga": load_tcga(), "gse76019": load_gse()}
(OUT / "formal_mito_analysis_summary.json").write_text(json.dumps(results, indent=2, default=str), encoding="utf-8")
print(json.dumps(results, indent=2, default=str))
