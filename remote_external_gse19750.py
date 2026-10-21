#!/usr/bin/env python3
"""Remote-only external bulk validation in GSE19750."""

import csv
import gzip
import json
import pathlib
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, kruskal, fisher_exact


ROOT = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path(
    "/root/autodl-tmp/acc-mitochondrial-mitotane"
)
RAW = ROOT / "data/raw/public_acc_bulk"
OUT = ROOT / "results/external_validation"
OUT.mkdir(parents=True, exist_ok=True)


def read_gene_set(path):
    genes = []
    for line in pathlib.Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("VERSION"):
            continue
        gene = line.split("\t")[-1].strip().upper()
        if re.fullmatch(r"[A-Z0-9-]+", gene):
            genes.append(gene)
    return sorted(set(genes))


SETS = {
    "OXPHOS": read_gene_set(ROOT / "data/raw/gene_sets/hallmark_oxphos.dat"),
    "ROS": read_gene_set(ROOT / "data/raw/gene_sets/hallmark_ros.dat"),
    "G2M": read_gene_set(ROOT / "data/raw/gene_sets/hallmark_g2m.dat"),
}


def bootstrap_median_difference(a, b, seed=20261021, n_bootstrap=10000):
    """Percentile CI for median(a) - median(b), preserving group sizes."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    rng = np.random.default_rng(seed)
    draws_a = rng.choice(a, size=(n_bootstrap, len(a)), replace=True)
    draws_b = rng.choice(b, size=(n_bootstrap, len(b)), replace=True)
    differences = np.median(draws_a, axis=1) - np.median(draws_b, axis=1)
    return [float(np.quantile(differences, 0.025)), float(np.quantile(differences, 0.975))]


def parse_metadata_and_matrix(path):
    metadata = {}
    table_lines = []
    in_table = False
    with gzip.open(path, "rt", errors="replace", newline="") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.lower() == "!series_matrix_table_begin":
                in_table = True
                continue
            if line.lower() == "!series_matrix_table_end":
                in_table = False
                continue
            if in_table:
                table_lines.append(line)
            elif line.startswith("!"):
                fields = next(csv.reader([line], delimiter="\t"))
                key = fields[0][1:]
                metadata.setdefault(key, []).append(fields[1:])

    if not table_lines:
        raise RuntimeError("GSE19750 filtered file contains no matrix table")
    table = pd.read_csv(
        pd.io.common.StringIO("\n".join(table_lines)),
        sep="\t",
        dtype=str,
        quotechar='"',
    )
    probe_col = table.columns[0]
    table[probe_col] = table[probe_col].astype(str).str.strip().str.strip('"')
    table = table.set_index(probe_col)
    table.columns = [str(x).strip().strip('"') for x in table.columns]
    table = table.apply(pd.to_numeric, errors="coerce")
    return metadata, table


def parse_platform(path):
    rows = []
    in_table = False
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line == "!platform_table_begin":
                in_table = True
                continue
            if line == "!platform_table_end":
                break
            if not in_table or not line:
                continue
            fields = line.split("\t")
            if fields[0] == "ID":
                continue
            if len(fields) < 3:
                continue
            gene = fields[2].split("///")[0].strip().upper()
            if gene and gene != "---":
                rows.append((fields[0].strip(), gene))
    if len(rows) < 30000:
        raise RuntimeError(f"GPL570 annotation parse failed: {len(rows)} rows")
    return pd.DataFrame(rows, columns=["probe", "gene"]).drop_duplicates("probe")


def score(expr):
    scores = pd.DataFrame(index=expr.columns)
    overlap = {}
    for name, requested in SETS.items():
        genes = sorted(set(requested) & set(expr.index))
        overlap[name] = {"requested": len(requested), "present": len(genes)}
        block = expr.loc[genes].astype(float)
        z = block.sub(block.mean(axis=1), axis=0).div(
            block.std(axis=1, ddof=1).replace(0, np.nan), axis=0
        )
        scores[name] = z.mean(axis=0)
    scores["MITO_COMPOSITE"] = scores[["OXPHOS", "ROS"]].mean(axis=1)
    return scores, overlap


def row_values(metadata, key, sample_ids, startswith=None):
    rows = metadata.get(key, [])
    if rows and not isinstance(rows[0], list):
        rows = [rows]
    if startswith is not None:
        rows = [
            row for row in rows
            if any(str(value).lower().startswith(startswith.lower()) for value in row)
        ]
    values = rows[0] if rows else []
    if len(values) != len(sample_ids):
        return pd.Series(index=sample_ids, dtype=object)
    return pd.Series(values, index=sample_ids)


def clean_value(value):
    value = str(value).strip().strip('"')
    return value.split(":", 1)[1].strip() if ":" in value else value


matrix_path = RAW / "GSE19750_filtered_data.txt.gz"
platform_path = RAW / "GPL570.annot.gz"
metadata, probe_expr = parse_metadata_and_matrix(matrix_path)
annotation = parse_platform(platform_path)
common = probe_expr.index.intersection(annotation["probe"])
mapped = annotation.set_index("probe").loc[common]
expr = probe_expr.loc[common].copy()
expr["gene"] = mapped["gene"].to_numpy()
expr = expr.groupby("gene").median()
scores, overlap = score(expr)

sample_ids = list(scores.index)
source = row_values(metadata, "Sample_source_name_ch1", sample_ids).str.lower()
sample_type = pd.Series(np.where(source.str.contains("normal"), "normal", "ACC"), index=sample_ids)
functional_raw = row_values(
    metadata, "Sample_characteristics_ch1", sample_ids, startswith="functional:"
)
functional = functional_raw.map(clean_value)
functional_status = functional.map(
    lambda x: np.nan if str(x).lower() in {"na", "unknown", "nan", ""}
    else ("nonfunctional" if str(x).lower() == "none" else "functional")
)

joined = scores.copy()
joined["sample_type"] = sample_type
joined["functional_raw"] = functional
joined["functional_status"] = functional_status
joined["sample_source"] = source
joined.to_csv(OUT / "GSE19750_external_scores.tsv", sep="\t")

acc = joined.loc[joined["sample_type"].eq("ACC")]
normal = joined.loc[joined["sample_type"].eq("normal")]
tests = {}
for feature in ["OXPHOS", "ROS", "MITO_COMPOSITE", "G2M"]:
    a = acc[feature].dropna()
    n = normal[feature].dropna()
    test = mannwhitneyu(a, n, alternative="two-sided")
    tests[f"ACC_vs_normal_{feature}"] = {
        "acc_n": int(len(a)),
        "normal_n": int(len(n)),
        "acc_median": float(a.median()),
        "normal_median": float(n.median()),
        "delta_acc_minus_normal": float(a.median() - n.median()),
        "U": float(test.statistic),
        "p": float(test.pvalue),
        "bootstrap_ci95_delta": bootstrap_median_difference(a, n),
    }

functional_frame = acc[["MITO_COMPOSITE", "G2M", "functional_status"]].dropna()
if functional_frame["functional_status"].nunique() == 2:
    low = functional_frame.loc[functional_frame["functional_status"].eq("nonfunctional"), "MITO_COMPOSITE"]
    high = functional_frame.loc[functional_frame["functional_status"].eq("functional"), "MITO_COMPOSITE"]
    test = mannwhitneyu(high, low, alternative="two-sided")
    tests["functional_vs_nonfunctional_MITO_COMPOSITE"] = {
        "functional_n": int(len(high)),
        "nonfunctional_n": int(len(low)),
        "functional_median": float(high.median()),
        "nonfunctional_median": float(low.median()),
        "delta_functional_minus_nonfunctional": float(high.median() - low.median()),
        "U": float(test.statistic),
        "p": float(test.pvalue),
        "bootstrap_ci95_delta": bootstrap_median_difference(high, low, seed=20261022),
    }

g2m_frame = acc[["MITO_COMPOSITE", "G2M"]].dropna()
if len(g2m_frame) >= 5:
    rho, p = spearmanr(g2m_frame["MITO_COMPOSITE"], g2m_frame["G2M"])
    tests["MITO_COMPOSITE_vs_G2M"] = {"n": int(len(g2m_frame)), "rho": float(rho), "p": float(p)}

state = pd.Series(np.where(
    acc["MITO_COMPOSITE"] >= acc["MITO_COMPOSITE"].median(), "MITO_HIGH", "MITO_LOW"
), index=acc.index)
state_functional = pd.crosstab(state, functional_status.loc[acc.index]).reindex(
    index=["MITO_LOW", "MITO_HIGH"], columns=["nonfunctional", "functional"], fill_value=0
)
state_test = None
if state_functional.shape == (2, 2):
    odds, p = fisher_exact(state_functional.to_numpy())
    state_test = {"table": state_functional.to_dict(), "odds_ratio": float(odds), "p": float(p)}

result = {
    "schema": "acc_mitochondrial_mitotane.gse19750_external.v1",
    "cohort": "GSE19750",
    "platform": "GPL570",
    "n_samples": int(len(scores)),
    "n_acc": int(len(acc)),
    "n_normal": int(len(normal)),
    "n_probe_rows": int(probe_expr.shape[0]),
    "n_mapped_probes": int(len(common)),
    "n_genes": int(expr.shape[0]),
    "gene_set_overlap": overlap,
    "functional_categories": functional.value_counts(dropna=False).to_dict(),
    "functional_status_counts_in_acc": functional_status.loc[acc.index].value_counts(dropna=False).to_dict(),
    "state_counts_acc_median_split": state.value_counts().to_dict(),
    "tests": tests,
    "state_by_functional": state_test,
    "cohort_role": "adult ACC bulk external expression validation; no direct mitotane response labels",
    "claim_boundary": "associative validation of mitochondrial expression state, not a mitotane sensitivity classifier",
}
(OUT / "GSE19750_external_summary.json").write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
print(json.dumps(result, indent=2, default=str))
