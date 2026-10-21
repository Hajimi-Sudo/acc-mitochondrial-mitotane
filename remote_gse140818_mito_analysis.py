#!/usr/bin/env python3
"""Remote-only mitochondrial analysis of the public GSE140818 HAC15 study."""

import gzip
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


ROOT = pathlib.Path(sys.argv[1])
MATRIX = ROOT / "data/raw/public_mitotane_geo/GSE140818_series_matrix.txt.gz"
ANNOTATION = ROOT / "data/raw/public_mitotane_geo/annotation/GPL15207_full.txt"
OUT = ROOT / "results/public_mitotane_geo"
OUT.mkdir(parents=True, exist_ok=True)


def parse_matrix(path):
    metadata = {}
    table = []
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
                table.append(line)
            elif line.startswith("!Sample_") or line.startswith("!Series_"):
                fields = line.split("\t")
                metadata.setdefault(fields[0][1:], []).append(
                    [x.strip().strip('"') for x in fields[1:]]
                )
    header = [x.strip().strip('"') for x in table[0].split("\t")]
    frame = pd.read_csv(pd.io.common.StringIO("\n".join(table)), sep="\t")
    frame.iloc[:, 0] = frame.iloc[:, 0].astype(str).str.strip('"')
    frame = frame.set_index(frame.columns[0])
    frame.columns = header[1:]
    return frame.apply(pd.to_numeric, errors="coerce"), metadata


def load_annotation():
    # GEO's current GPL table is the authoritative annotation for this array.
    # The table is large, so retain only the probe and Gene Symbol columns.
    begin = False
    header = None
    gene_symbol_index = -1
    rows = []
    with ANNOTATION.open("r", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line == "!platform_table_begin":
                begin = True
                continue
            if line == "!platform_table_end":
                break
            if not begin or not line:
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                gene_symbol_index = header.index("Gene Symbol") if "Gene Symbol" in header else -1
                continue
            if len(fields) != len(header):
                continue
            rows.append((fields[0], fields[gene_symbol_index]))
    if not rows or header is None or gene_symbol_index < 0 or len(rows) < 40000:
        raise RuntimeError(f"GPL15207 table parse failed: rows={len(rows)}, header={header[:15] if header else None}")
    frame = pd.DataFrame(rows, columns=["probe", "gene"])
    frame["probe"] = frame["probe"].astype(str).str.strip()
    frame["gene"] = frame["gene"].astype(str).str.split("///").str[0].str.upper().str.strip()
    return frame[frame["gene"].ne("") & frame["gene"].ne("---")].drop_duplicates("probe")


def score(expr, genesets):
    scores = pd.DataFrame(index=expr.columns)
    overlap = {}
    for name, genes in genesets.items():
        present = sorted(set(genes) & set(expr.index))
        overlap[name] = {"requested": len(genes), "present": len(present), "genes": present}
        block = expr.loc[present]
        normalized = block.sub(block.mean(axis=1), axis=0).div(block.std(axis=1).replace(0, np.nan), axis=0)
        scores[name] = normalized.mean(axis=0)
    scores["MITO_COMPOSITE"] = scores[["OXPHOS", "ROS"]].mean(axis=1)
    return scores, overlap


def compare(values, groups, feature):
    treated = values.loc[groups == "mitotane", feature].dropna()
    control = values.loc[groups == "control", feature].dropna()
    test = mannwhitneyu(treated, control, alternative="two-sided")
    return {
        "feature": feature,
        "treated_n": int(len(treated)),
        "control_n": int(len(control)),
        "treated_median": float(treated.median()),
        "control_median": float(control.median()),
        "delta_median": float(treated.median() - control.median()),
        "U": float(test.statistic),
        "p": float(test.pvalue),
    }


def compare_subset(values, groups, subset, feature):
    keep = subset.reindex(values.index).fillna(False)
    return compare(values.loc[keep], groups.loc[keep], feature)


expr_probe, metadata = parse_matrix(MATRIX)
annotation = load_annotation()
common = expr_probe.index.intersection(annotation["probe"])
mapped = annotation.set_index("probe").loc[common]
expr_gene = expr_probe.loc[common].copy()
expr_gene["gene"] = mapped["gene"].to_numpy()
expr_gene = expr_gene.groupby("gene").median()

sets = {
    "OXPHOS": pathlib.Path(ROOT / "data/raw/gene_sets/hallmark_oxphos.dat").read_text().splitlines()[2:],
    "ROS": pathlib.Path(ROOT / "data/raw/gene_sets/hallmark_ros.dat").read_text().splitlines()[2:],
    "G2M": pathlib.Path(ROOT / "data/raw/gene_sets/hallmark_g2m.dat").read_text().splitlines()[2:],
}
sets = {k: sorted(set(x.strip().upper() for x in v if x.strip())) for k, v in sets.items()}
scores, overlap = score(expr_gene, sets)

sample_ids = metadata["Sample_geo_accession"][0]
titles = metadata["Sample_title"][0]
descriptions = metadata["Sample_description"][0]
char_rows = metadata.get("Sample_characteristics_ch1", [])
therapy = next((row for row in char_rows if any(str(x).lower().startswith("treatment:") for x in row)), None)
if therapy is None or len(therapy) != len(sample_ids):
    raise RuntimeError("GSE140818 treatment metadata could not be resolved sample-by-sample")
groups = pd.Series(index=sample_ids, dtype=object)
clone = pd.Series(index=sample_ids, dtype=object)
for sid, title, desc, treatment in zip(sample_ids, titles, descriptions, therapy):
    text = f"{title} {desc}".lower()
    groups.loc[sid] = "mitotane" if "mitotane" in str(treatment).lower() else "control"
    if "mitotane resistant" in text:
        clone.loc[sid] = "resistant"
    elif "nonresistant" in text:
        clone.loc[sid] = "nonresistant"
    else:
        raise RuntimeError(f"unresolved clone label for {sid}: {desc}")
scores.index = sample_ids

results = {
    "cohort": "GSE140818",
    "platform": "GPL15207_PrimeView",
    "n_probe_rows": int(expr_probe.shape[0]),
    "n_mapped_probes": int(len(common)),
    "n_genes": int(expr_gene.shape[0]),
    "gene_set_overlap": {k: {"requested": v["requested"], "present": v["present"]} for k, v in overlap.items()},
    "design": {
        "acute_treatment": {"mitotane": int((groups == "mitotane").sum()), "control": int((groups == "control").sum())},
        "long_term_clone": {"resistant": int((clone == "resistant").sum()), "nonresistant": int((clone == "nonresistant").sum())},
    },
    "acute_treatment_tests_all_clones": [compare(scores, groups, feature) for feature in ["OXPHOS", "ROS", "MITO_COMPOSITE", "G2M"]],
    "acute_treatment_tests_nonresistant": [
        compare_subset(scores, groups, clone.eq("nonresistant"), feature)
        for feature in ["OXPHOS", "ROS", "MITO_COMPOSITE", "G2M"]
    ],
    "acute_treatment_tests_resistant": [
        compare_subset(scores, groups, clone.eq("resistant"), feature)
        for feature in ["OXPHOS", "ROS", "MITO_COMPOSITE", "G2M"]
    ],
    "dmso_baseline_resistant_vs_nonresistant": [
        compare_subset(
            scores,
            clone.rename("clone").replace({"resistant": "mitotane", "nonresistant": "control"}),
            groups.eq("control"),
            feature,
        )
        for feature in ["OXPHOS", "ROS", "MITO_COMPOSITE", "G2M"]
    ],
    "claim_boundary": "HAC15 cell-model treatment and resistance context, not patient IC50/AUC or clinical response",
}
scores.assign(treatment=groups, clone=clone).to_csv(OUT / "GSE140818_HAC15_mito_scores.tsv", sep="\t")
(OUT / "GSE140818_HAC15_mito_summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
print(json.dumps(results, indent=2))
