#!/usr/bin/env python3
"""Build a small, claim-bound evidence manifest on the remote server."""

import hashlib
import json
import pathlib
import sys


ROOT = pathlib.Path(sys.argv[1])
OUT = ROOT / "results/evidence"
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load(rel):
    path = ROOT / rel
    return json.loads(path.read_text(encoding="utf-8"))


gse140 = load("results/public_mitotane_geo/GSE140818_HAC15_effect_stats.json")
gse140_adj = load("results/public_mitotane_geo/GSE140818_HAC15_adjusted_stats.json")
gse140_factorial = load("results/public_mitotane_geo/GSE140818_HAC15_factorial_stats.json")
tcga_robustness = load("results/robustness/tcga_robustness_summary.json")
gse734 = load("results/public_mitotane_geo/GSE73417_H295R_xenograft_mito_summary.json")
tcga = load("results/formal/formal_mito_analysis_summary.json")

manifest = {
    "schema": "acc_mitochondrial_mitotane.evidence_manifest.v1",
    "project_root": str(ROOT),
    "execution_boundary": "remote_only",
    "claim_status": {
        "patient_level_mitotane_sensitivity": "blocked_no_sample_matched_IC50_AUC_or_response_labels",
        "mitochondrial_state_reproducibility": "supported_as_expression_state",
        "endocrine_function_association": "exploratory_not_independent_after_G2M_sensitivity",
        "cell_model_context": "supported_bounded",
        "xenograft_context": "supported_bounded_null_mitochondrial_contrast",
    },
    "cohort_roles": {
        "TCGA-ACC": "human discovery and endocrine/clinical association",
        "GSE76019": "pediatric expression/outcome context only",
        "GSE140818": "HAC15 acute treatment and long-term resistance context",
        "GSE186870": "underpowered HAC15 treatment context",
        "GSE73417": "H295R xenograft post-treatment context",
        "HRA006596": "controlled-access candidate; not analyzed",
    },
    "key_results": {
        "TCGA_state_counts": tcga.get("tcga", {}).get("state_counts", {}),
        "GSE140818_effect_stats": gse140,
        "GSE140818_G2M_adjusted": gse140_adj,
        "GSE140818_factorial": gse140_factorial,
        "TCGA_robustness": tcga_robustness,
        "GSE73417_mitotane_vs_control": gse734.get("mitotane_vs_control", []),
        "GSE73417_group_counts": gse734.get("group_counts", {}),
    },
    "source_artifacts": {},
}

for rel in [
    "results/formal/formal_mito_analysis_summary.json",
    "results/public_mitotane_geo/GSE140818_HAC15_effect_stats.json",
    "results/public_mitotane_geo/GSE140818_HAC15_adjusted_stats.json",
    "results/public_mitotane_geo/GSE140818_HAC15_factorial_stats.json",
    "results/robustness/tcga_robustness_summary.json",
    "results/public_mitotane_geo/GSE73417_H295R_xenograft_mito_summary.json",
    "data/raw/public_mitotane_geo/annotation/GPL13497_full.txt",
]:
    path = ROOT / rel
    manifest["source_artifacts"][rel] = {
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }

(OUT / "evidence_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
print(json.dumps(manifest, indent=2))
