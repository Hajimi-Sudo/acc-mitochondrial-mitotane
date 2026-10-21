#!/usr/bin/env python3
"""Render reviewed aggregate tables to publication JPGs on the remote server.

This is manuscript packaging only. It does not read raw matrices or perform
new statistical analysis.
"""

import csv
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = pathlib.Path(sys.argv[1])
TABLES = ROOT / "results/publication_package/tables"
OUT = ROOT / "results/publication_package/figures_jpg"
OUT.mkdir(parents=True, exist_ok=True)

BLUE = "#2f6690"
ORANGE = "#d95f02"
TEAL = "#86b6c9"
GREEN = "#6a994e"
GRAY = "#8c98a4"


def read_tsv(name):
    with (TABLES / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def num(row, key):
    return float(row[key])


def save(fig, name):
    fig.savefig(OUT / f"{name}.jpg", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def forest(rows, labels, estimate, low, high, pkey, title, xlabel, name, xlim=None):
    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    y = np.arange(len(rows))
    est = np.array([num(r, estimate) for r in rows])
    lo = np.array([num(r, low) for r in rows])
    hi = np.array([num(r, high) for r in rows])
    ax.errorbar(est, y, xerr=[est - lo, hi - est], fmt="o", color=ORANGE,
                ecolor=BLUE, elinewidth=2.2, capsize=3, markersize=6)
    ax.axvline(0, color="#555555", linestyle="--", linewidth=1)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left", weight="bold")
    ax.grid(axis="x", color="#d8e0e5", linewidth=0.8)
    for idx, row in enumerate(rows):
        ax.text(1.01, idx, f"p={num(row, pkey):.3f}", transform=ax.get_yaxis_transform(), va="center", fontsize=8)
    if xlim:
        ax.set_xlim(*xlim)
    fig.tight_layout()
    save(fig, name)


cohorts = read_tsv("cohort_roles.tsv")
cohort_values = [12 if r["cohort"] == "GSE186870" else (0 if r["n"] == "Not analyzed" else float(str(r["n"]).split()[0])) for r in cohorts]
fig, ax = plt.subplots(figsize=(8.6, 4.5))
labels = [r["cohort"] for r in cohorts]
colors = [BLUE, TEAL, ORANGE, GRAY, GREEN, "#a0a7ad"]
ax.barh(labels[::-1], cohort_values[::-1], color=colors[::-1])
for y, value in enumerate(cohort_values[::-1]):
    ax.text(value + 1, y, str(int(value)), va="center", fontsize=9)
ax.set_xlabel("Samples or profiled specimens")
ax.set_title("Cohort roles and evidence boundaries", loc="left", weight="bold")
ax.grid(axis="x", color="#d8e0e5", linewidth=0.8)
fig.tight_layout()
save(fig, "Figure1_cohort_roles")

state = [("MITO_HIGH", 65), ("MITO_LOW", 14)]
fig, ax = plt.subplots(figsize=(6.2, 4.2))
ax.bar([x[0] for x in state], [x[1] for x in state], color=[BLUE, ORANGE], width=0.55)
ax.set_ylabel("Samples")
ax.set_title("TCGA-ACC mitochondrial state composition", loc="left", weight="bold")
ax.grid(axis="y", color="#d8e0e5", linewidth=0.8)
for x, (_, value) in enumerate(state):
    ax.text(x, value + 2, str(value), ha="center", weight="bold")
fig.tight_layout()
save(fig, "Figure2_TCGA_state_counts")

gse140 = read_tsv("gse140818_contrasts.tsv")
forest(gse140, [f'{r["feature"]}: {r["contrast"].replace(" DMSO resistant-minus-nonresistant", "")}' for r in gse140], "estimate", "ci_low", "ci_high", "permutation_p", "HAC15 resistance association and G2M sensitivity", "Median score difference", "Figure3_GSE140818_G2M_sensitivity", (-0.8, 0.2))

gse734 = read_tsv("gse73417_mitotane_vs_control.tsv")
forest(gse734, [r["feature"] for r in gse734], "estimate", "ci_low", "ci_high", "permutation_p", "H295R xenograft treatment context", "Median score difference", "Figure4_GSE73417_xenograft_context", (-1, 1))

gse186 = read_tsv("gse186870_hac15_context.tsv")
fig, ax = plt.subplots(figsize=(7.2, 4.2))
ax.bar([r["feature"] for r in gse186], [num(r, "delta_median") for r in gse186], color=[BLUE, TEAL, ORANGE, GRAY])
ax.axhline(0, color="#555555", linewidth=1)
ax.set_ylabel("Treated minus control median score")
ax.set_title("GSE186870 HAC15 treatment context", loc="left", weight="bold")
ax.grid(axis="y", color="#d8e0e5", linewidth=0.8)
fig.tight_layout()
save(fig, "Figure5_GSE186870_HAC15_context")

hormone = read_tsv("tcga_hormone_g2m_sensitivity.tsv")
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.bar([r["model"] for r in hormone], [num(r, "mito_standardized_coefficient") for r in hormone], color=[BLUE, ORANGE], width=0.55)
ax.axhline(0, color="#555555", linewidth=1)
ax.set_ylabel("Standardized MITO_COMPOSITE coefficient")
ax.set_title("Endocrine association attenuates after G2M adjustment", loc="left", weight="bold")
ax.grid(axis="y", color="#d8e0e5", linewidth=0.8)
fig.tight_layout()
save(fig, "Figure6_TCGA_hormone_G2M_sensitivity")

robust = read_tsv("tcga_state_construction_sensitivity.tsv")
fig, ax = plt.subplots(figsize=(7.6, 4.2))
ax.barh([r["variant"] for r in robust][::-1], [num(r, "overlap_with_reference") for r in robust][::-1], color=BLUE)
ax.set_xlim(0, 1.05)
ax.set_xlabel("Overlap with reference two-score labels")
ax.set_title("TCGA mitochondrial-state construction sensitivity", loc="left", weight="bold")
ax.grid(axis="x", color="#d8e0e5", linewidth=0.8)
fig.tight_layout()
save(fig, "Figure7_TCGA_state_robustness")

factor = read_tsv("gse140818_factorial_interactions.tsv")
forest(factor, [f'{r["feature"]} {r["model"]}' for r in factor], "estimate", "ci_low", "ci_high", "permutation_p", "GSE140818 mitochondrial-by-treatment interaction", "Interaction estimate", "Figure8_GSE140818_factorial_interaction", (-0.55, 0.55))

gse19750 = read_tsv("gse19750_adult_acc_external_validation.tsv")
forest(gse19750, [r["feature"] for r in gse19750], "delta_acc_minus_normal", "ci_low", "ci_high", "p_value", "Adult ACC external expression context", "Median difference: ACC minus normal adrenal", "Figure9_GSE19750_adult_external_context", (-0.9, 0.65))

print(f"Rendered {len(list(OUT.glob('*.jpg')))} aggregate JPG figures to {OUT}")
