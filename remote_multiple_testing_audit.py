#!/usr/bin/env python3
"""Remote-only multiple-testing audit for the ACC mitochondrial project.

Reads only reviewed aggregate TSV outputs on the remote project and writes
aggregate q-value tables plus a machine-readable audit. No raw or patient-level
data are accessed.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path("/root/autodl-tmp/acc-mitochondrial-mitotane")
TABLES = ROOT / "results/publication_package/tables"
OUT = ROOT / "results/multiple_testing"
OUT.mkdir(parents=True, exist_ok=True)


def read_tsv(name: str) -> list[dict[str, str]]:
    with (TABLES / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def bh(p_values: list[float]) -> list[float]:
    order = sorted(range(len(p_values)), key=p_values.__getitem__)
    q = [1.0] * len(p_values)
    running = 1.0
    for rank, index in reversed(list(enumerate(order, start=1))):
        running = min(running, p_values[index] * len(p_values) / rank)
        q[index] = min(1.0, running)
    return q


def add_q(rows: list[dict[str, str]], p_key: str, q_key: str = "bh_q_value") -> list[dict[str, str]]:
    p_values = [float(row[p_key]) for row in rows]
    q_values = bh(p_values)
    for row, q_value in zip(rows, q_values):
        row[q_key] = f"{q_value:.12g}"
    return rows


def write_tsv(name: str, rows: list[dict[str, str]]) -> None:
    fields = list(rows[0])
    with (OUT / name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


families: dict[str, dict[str, object]] = {}


def process(name: str, output: str, family: str, p_key: str = "permutation_p") -> None:
    all_rows = read_tsv(name)
    rows = [row for row in all_rows if row.get(p_key, "").strip()]
    for row in all_rows:
        if not row.get(p_key, "").strip():
            row["bh_q_value"] = "NA"
    if rows:
        rows = add_q(rows, p_key)
    rows = rows + [row for row in all_rows if not row.get(p_key, "").strip()]
    write_tsv(output, rows)
    families[family] = {
        "source_table": name,
        "output_table": output,
        "n_tests": len([row for row in all_rows if row.get(p_key, "").strip()]),
        "n_descriptive_rows_excluded": len(all_rows) - len([row for row in all_rows if row.get(p_key, "").strip()]),
        "p_values": [float(row[p_key]) for row in all_rows if row.get(p_key, "").strip()],
        "q_values": [float(row["bh_q_value"]) for row in rows if row.get(p_key, "").strip()],
    }


process(
    "gse19750_adult_acc_external_validation.tsv",
    "gse19750_adult_acc_external_validation_bh.tsv",
    "GSE19750_ACC_vs_normal_four_prespecified_scores",
    "p_value",
)
process(
    "tcga_state_and_correlation.tsv",
    "tcga_state_and_correlation_bh.tsv",
    "TCGA_two_prespecified_score_correlations",
    "p_value",
)
resistance_rows = read_tsv("gse140818_contrasts.tsv")
for label, marker, output in [
    ("raw", "raw", "gse140818_contrasts_raw_bh.tsv"),
    ("G2M_adjusted", "G2M-adjusted", "gse140818_contrasts_g2m_adjusted_bh.tsv"),
]:
    rows = [row for row in resistance_rows if marker in row["contrast"]]
    rows = add_q(rows, "permutation_p")
    write_tsv(output, rows)
    families[f"GSE140818_DMSO_resistance_{label}_reported_scores"] = {
        "source_table": "gse140818_contrasts.tsv",
        "output_table": output,
        "n_tests": len(rows),
        "p_values": [float(row["permutation_p"]) for row in rows],
        "q_values": [float(row["bh_q_value"]) for row in rows],
    }

factorial_rows = read_tsv("gse140818_factorial_interactions.tsv")
for label, model, output in [
    ("raw", "raw", "gse140818_factorial_interactions_raw_bh.tsv"),
    ("G2M_adjusted", "G2M_adjusted", "gse140818_factorial_interactions_g2m_adjusted_bh.tsv"),
]:
    rows = [
        row for row in factorial_rows
        if row["model"] == model and row["feature"] != "G2M"
    ]
    rows = add_q(rows, "permutation_p")
    write_tsv(output, rows)
    families[f"GSE140818_{label}_factorial_three_mitochondrial_scores"] = {
        "source_table": "gse140818_factorial_interactions.tsv",
        "output_table": output,
        "n_tests": len(rows),
        "p_values": [float(row["permutation_p"]) for row in rows],
        "q_values": [float(row["bh_q_value"]) for row in rows],
    }

process(
    "gse73417_mitotane_vs_control.tsv",
    "gse73417_mitotane_vs_control_bh.tsv",
    "GSE73417_mitotane_vs_control_four_scores",
)
process(
    "gse186870_hac15_context.tsv",
    "gse186870_HAC15_treated_vs_control_bh.tsv",
    "GSE186870_HAC15_treated_vs_control_four_scores",
    "p_value",
)
process(
    "gse76019_context.tsv",
    "gse76019_context_bh.tsv",
    "GSE76019_two_exploratory_outcome_associations",
    "p_value",
)

audit = {
    "schema": "acc_mitochondrial_mitotane.multiple_testing_audit.v1",
    "execution": "remote_only",
    "method": "Benjamini-Hochberg within prespecified endpoint/design families",
    "alpha": 0.05,
    "primary_analysis": "TCGA state construction and robustness were prespecified descriptive/stability analyses; no primary inferential p-value was used for the state-construction endpoint.",
    "exploratory_policy": "Nominal p-values are retained, but inferential interpretation uses family-specific BH q-values and effect sizes with 95% bootstrap intervals.",
    "families": families,
    "limitations": [
        "The family definitions were fixed by cohort, endpoint, and design before interpreting the q-values.",
        "The audit does not create a patient-level mitotane-response endpoint.",
        "The TCGA endocrine coefficient table contains descriptive in-sample models without an inferential p-value in the reviewed aggregate table; those models remain exploratory and are not upgraded by multiplicity correction.",
    ],
}
(OUT / "multiple_testing_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
print(json.dumps(audit, indent=2))
