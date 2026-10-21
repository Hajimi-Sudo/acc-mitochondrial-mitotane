#!/usr/bin/env python3
"""Remote-only mitochondrial scoring for the GSE73417 H295R xenograft study."""

import gzip
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


ROOT = pathlib.Path(sys.argv[1])
MATRIX = ROOT / "data/raw/public_mitotane_geo/GSE73417_series_matrix.txt.gz"
ANNOTATION = ROOT / "data/raw/public_mitotane_geo/annotation/GPL13497_full.txt"
OUT = ROOT / "results/public_mitotane_geo"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261021
RNG = np.random.default_rng(SEED)


def parse_matrix(path):
    metadata = {}
    rows = []
    in_table = False
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                continue
            if line == "!series_matrix_table_end":
                in_table = False
                continue
            if in_table:
                rows.append(line)
            elif line.startswith("!Sample_") or line.startswith("!Series_"):
                fields = line.split("\t")
                metadata[fields[0][1:]] = [x.strip().strip('"') for x in fields[1:]]
    frame = pd.read_csv(pd.io.common.StringIO("\n".join(rows)), sep="\t")
    header = [x.strip().strip('"') for x in rows[0].split("\t")]
    frame.iloc[:, 0] = frame.iloc[:, 0].astype(str).str.strip('"')
    frame = frame.set_index(frame.columns[0])
    frame.columns = header[1:]
    return frame.apply(pd.to_numeric, errors="coerce"), metadata


def load_annotation(path):
    started = False
    header = None
    gene_index = None
    rows = []
    with path.open("r", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line == "!platform_table_begin":
                started = True
                continue
            if line == "!platform_table_end":
                break
            if not started or not line:
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                gene_index = header.index("GENE_SYMBOL")
                continue
            if len(fields) == len(header):
                rows.append((fields[0], fields[gene_index]))
    # GPL13497 contains a filtered ~12.8k-probe platform table.
    if len(rows) < 10000:
        raise RuntimeError(f"GPL13497 annotation parse failed: {len(rows)} rows")
    annotation = pd.DataFrame(rows, columns=["probe", "gene"])
    annotation["gene"] = annotation["gene"].astype(str).str.split("///").str[0].str.upper().str.strip()
    return annotation[annotation["gene"].ne("") & annotation["gene"].ne("---")].drop_duplicates("probe")


def read_sets(root):
    sets = {}
    for name in ("hallmark_oxphos", "hallmark_ros", "hallmark_g2m"):
        lines = (root / f"data/raw/gene_sets/{name}.dat").read_text().splitlines()[2:]
        sets[name.removeprefix("hallmark_").upper()] = sorted({x.strip().upper() for x in lines if x.strip()})
    return sets


def score(expr, gene_sets):
    scores = pd.DataFrame(index=expr.columns)
    overlap = {}
    for name, genes in gene_sets.items():
        present = sorted(set(genes) & set(expr.index))
        overlap[name] = {"requested": len(genes), "present": len(present), "genes": present}
        block = expr.loc[present]
        z = block.sub(block.mean(axis=1), axis=0).div(block.std(axis=1).replace(0, np.nan), axis=0)
        scores[name] = z.mean(axis=0)
    scores["MITO_COMPOSITE"] = scores[["OXPHOS", "ROS"]].mean(axis=1)
    return scores, overlap


def contrast(frame, a_label, b_label, feature, n_perm=10000, n_boot=5000):
    a = frame.loc[frame["group"].eq(a_label), feature].dropna().to_numpy(float)
    b = frame.loc[frame["group"].eq(b_label), feature].dropna().to_numpy(float)
    observed = float(np.median(a) - np.median(b))
    pooled = np.concatenate([a, b])
    null = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        perm = RNG.permutation(pooled)
        null[i] = np.median(perm[: len(a)]) - np.median(perm[len(a) :])
    boot = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        boot[i] = np.median(RNG.choice(a, len(a), replace=True)) - np.median(RNG.choice(b, len(b), replace=True))
    test = mannwhitneyu(a, b, alternative="two-sided")
    return {
        "feature": feature,
        "group_a": a_label,
        "group_b": b_label,
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(np.median(a)),
        "median_b": float(np.median(b)),
        "median_difference_a_minus_b": observed,
        "bootstrap_ci95": [float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))],
        "permutation_p": float((np.sum(np.abs(null) >= abs(observed)) + 1) / (n_perm + 1)),
        "mann_whitney_p": float(test.pvalue),
        "n_permutations": n_perm,
        "n_bootstrap": n_boot,
    }


expr_probe, metadata = parse_matrix(MATRIX)
annotation = load_annotation(ANNOTATION)
common = expr_probe.index.intersection(annotation["probe"])
mapped = annotation.set_index("probe").loc[common]
expr_gene = expr_probe.loc[common].copy()
expr_gene["gene"] = mapped["gene"].to_numpy()
expr_gene = expr_gene.groupby("gene").median()
scores, overlap = score(expr_gene, read_sets(ROOT))

sample_ids = metadata["Sample_geo_accession"]
titles = metadata["Sample_title"]
descriptions = metadata["Sample_description"]
source = metadata["Sample_source_name_ch1"]
groups = []
for title, desc, source_name in zip(titles, descriptions, source):
    text = f"{title} {desc} {source_name}".lower()
    if "combined" in text:
        groups.append("combined")
    elif "mitotane" in text:
        groups.append("mitotane")
    elif "control" in text or "corn oil" in text:
        groups.append("control")
    elif "9-cisra" in text or "9cisra" in text:
        groups.append("9-cisRA")
    else:
        raise RuntimeError(f"Unresolved treatment group: {title} | {desc}")
scores.index = sample_ids
scores["group"] = groups
scores.to_csv(OUT / "GSE73417_H295R_xenograft_mito_scores.tsv", sep="\t")

features = ["OXPHOS", "ROS", "MITO_COMPOSITE", "G2M"]
result = {
    "cohort": "GSE73417",
    "title": metadata["Series_title"][0],
    "platform": metadata["Series_platform_id"][0],
    "sample_organism": sorted(set(metadata["Sample_organism_ch1"])),
    "model": "human H295R adrenocortical carcinoma xenografts in SCID mice",
    "n_probe_rows": int(expr_probe.shape[0]),
    "n_mapped_probes": int(len(common)),
    "n_genes": int(expr_gene.shape[0]),
    "group_counts": pd.Series(groups).value_counts().to_dict(),
    "gene_set_overlap": {k: {"requested": v["requested"], "present": v["present"]} for k, v in overlap.items()},
    "mitotane_vs_control": [contrast(scores, "mitotane", "control", feature) for feature in features],
    "combined_vs_control": [contrast(scores, "combined", "control", feature) for feature in features],
    "claim_boundary": "in-vivo xenograft treatment context; not a patient-level response cohort or clinical validation",
    "limitations": [
        "small group sizes (four samples per treatment group)",
        "xenograft model rather than human clinical specimens",
        "expression was measured after treatment and cannot establish pretreatment predictive sensitivity",
    ],
}
(OUT / "GSE73417_H295R_xenograft_mito_summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
