#!/usr/bin/env python
import collections
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
out = root / "data/raw/mitotane_paper/gsa_hra006596"
individuals = json.loads((out / "individuals.json").read_text())
runs = json.loads((out / "runs.json").read_text())
individual_views = individuals.get("individualViews") or []
run_views = runs.get("runViews") or []
strategy_counts = collections.Counter(x.get("strategyName") for x in run_views)
samples = sorted({x.get("sampleAcc") for x in run_views if x.get("sampleAcc")})
experiments = sorted({x.get("expName") for x in run_views if x.get("expName")})
summary = {
    "accession": "HRA006596",
    "title": "Genetic features of adrenocortical carcinoma correlate with in vitro sensitivity to mitotane",
    "access": "controlled_access",
    "bioproject": "PRJCA022756",
    "request_path": "/gsa-human/browse/request/HRA006596",
    "public_individual_count": len(individual_views),
    "public_run_file_rows": len(run_views),
    "unique_samples_in_runs": len(samples),
    "unique_experiments": len(experiments),
    "strategy_counts": dict(strategy_counts),
    "sample_accessions": samples,
    "experiment_accessions": experiments,
    "has_public_processed_files": bool(runs.get("processedFiles")),
    "raw_data_download_permitted": False,
    "direct_mitotane_endpoint_in_public_metadata": False,
    "blocker": "DAC authorization is required before retrieving RNA-seq/WES files; no access request was submitted.",
}
target = root / "results/pilot/gsa_hra006596_access_audit.json"
target.write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
