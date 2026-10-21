#!/usr/bin/env python3
"""Audit public mitotane-related GEO matrices on the remote server."""

import csv
import gzip
import hashlib
import json
import re
from pathlib import Path


ROOT = Path("/root/autodl-tmp/acc-mitochondrial-mitotane")
MATRIX_DIR = ROOT / "data/raw/public_mitotane_geo"
OUT = ROOT / "results/pilot/public_mitotane_data_availability.json"


def fields(line):
    return next(csv.reader([line], delimiter="\t"))


def clean(value):
    value = value.strip().strip('"')
    return re.sub(r"\\s+", " ", value)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit_matrix(path):
    metadata = {}
    samples = {}
    n_rows = 0
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
                if line.startswith('"ID_REF"'):
                    header = fields(line)
                    samples["sample_ids"] = header[1:]
                elif line and not line.startswith("!"):
                    n_rows += 1
                continue
            if line.startswith("!Series_"):
                key, value = line.split("\t", 1)
                metadata[key[1:]] = clean(value)
            elif line.startswith("!Sample_"):
                parts = fields(line)
                key = parts[0][1:]
                samples[key] = [clean(value) for value in parts[1:]]
    sample_ids = samples.get("sample_ids", [])
    treatment_rows = []
    for key, values in samples.items():
        if key.startswith("Sample_characteristics") or key in {"Sample_title", "Sample_description"}:
            text = " ".join(values).lower()
            if any(term in text for term in ("treatment", "mitotane", "dmso", "control", "9-cisra", "cisra")):
                treatment_rows.append((key, values))
    treatments = {}
    for key, values in treatment_rows:
        is_characteristic = key.startswith("Sample_characteristics")
        is_explicit_treatment = is_characteristic and any(
            value.lower().startswith("treatment:") for value in values
        )
        for sample_id, value in zip(sample_ids, values):
            if is_explicit_treatment:
                treatments[sample_id] = value
            elif sample_id not in treatments and key in {"Sample_title", "Sample_description"}:
                if any(term in value.lower() for term in ("mitotane", "control", "dmso", "9-cisra", "cisra")):
                    treatments[sample_id] = value
    treatment_values = [treatments.get(sample_id, "") for sample_id in sample_ids]
    return {
        "accession": metadata.get("Series_geo_accession"),
        "title": metadata.get("Series_title"),
        "summary": metadata.get("Series_summary"),
        "overall_design": metadata.get("Series_overall_design"),
        "platform": metadata.get("Series_platform_id"),
        "matrix_sha256": sha256(path),
        "matrix_bytes": path.stat().st_size,
        "n_probe_rows": n_rows,
        "n_samples": len(sample_ids),
        "sample_ids": sample_ids,
        "sample_titles": samples.get("Sample_title", []),
        "sample_source": samples.get("Sample_source_name_ch1", []),
        "sample_treatments": treatment_values,
        "sample_organisms": sorted(set(samples.get("Sample_organism_ch1", []))),
        "sample_descriptions": samples.get("Sample_description", []),
    }


def main():
    datasets = {}
    for path in sorted(MATRIX_DIR.glob("GSE*_series_matrix.txt.gz")):
        datasets[path.stem.split("_")[0]] = audit_matrix(path)
    report = {
        "execution_location": "remote_only",
        "remote_project_root": str(ROOT),
        "public_direct_response_cohort_verified": False,
        "datasets": datasets,
        "evidence_layers": {
            "GSE140818": "human HAC15 cell-line mitotane-resistant clones and acute 50 uM mitotane/DMSO expression profiling; not patient IC50 data",
            "GSE186870": "human HAC15 cell-line mitotane treatment/control expression profiling; treatment-response transcript context, not patient IC50 data",
            "GSE73417": "mouse H295R xenograft mitotane/9-cisRA treatment expression profiling; species-mismatched mechanistic context only",
        },
        "blocked_claim": "No public sample-matched patient transcriptome plus mitotane IC50/AUC/responder label has been verified. HRA006596 remains controlled-access.",
        "next_analysis": "Use TCGA-ACC for human endocrine-state discovery; use GSE140818 and GSE186870 as cell-model mitotane perturbation/resistance context; use GSE73417 only as mouse mechanistic context.",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
