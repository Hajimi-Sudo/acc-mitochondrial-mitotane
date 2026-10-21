#!/usr/bin/env python3
"""Remote-only QC for public mitotane-related GEO series matrices."""

import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/root/autodl-tmp/acc-mitochondrial-mitotane")
MATRIX_DIR = ROOT / "data/raw/public_mitotane_geo"
OUT = ROOT / "results/pilot/public_mitotane_geo_qc.json"


def parse_series(path):
    metadata = {}
    sample = {}
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
            elif line.startswith("!Series_") or line.startswith("!Sample_"):
                parts = line.split("\t")
                key = parts[0][1:]
                values = [x.strip().strip('"') for x in parts[1:]]
                if key.startswith("Series_"):
                    metadata[key] = " ".join(values)
                else:
                    sample[key] = values
    header = [x.strip().strip('"') for x in table[0].split("\t")]
    frame = pd.read_csv(pd.io.common.StringIO("\n".join(table)), sep="\t", dtype={0: str})
    frame = frame.replace({"null": np.nan, "": np.nan})
    frame.iloc[:, 0] = frame.iloc[:, 0].astype(str).str.strip('"')
    values = frame.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    return metadata, sample, values, header[1:]


def treatment_groups(sample, sample_ids):
    titles = dict(zip(sample_ids, sample.get("Sample_title", [""] * len(sample_ids))))
    descriptions = dict(zip(sample_ids, sample.get("Sample_description", [""] * len(sample_ids))))
    chars = " ".join(sample.get("Sample_characteristics_ch1", []))
    groups = {}
    for sid in sample_ids:
        text = f"{titles.get(sid, '')} {descriptions.get(sid, '')}".lower()
        if "mitotane" in text and "control" not in text and "dmso" not in text:
            groups[sid] = "mitotane"
        elif "dmso" in text or "control" in text or "corn oil" in text:
            groups[sid] = "control"
        elif "9-cisra" in text or "9cisra" in text:
            groups[sid] = "9-cisRA"
        else:
            groups[sid] = "unclassified"
    if "treatment:" in chars.lower():
        for sid, value in zip(sample_ids, sample.get("Sample_characteristics_ch1", [])):
            low = value.lower()
            if "mitotane" in low:
                groups[sid] = "mitotane"
            elif "dmso" in low or "control" in low or "corn oil" in low:
                groups[sid] = "control"
    return groups


def audit(path):
    metadata, sample, values, sample_ids = parse_series(path)
    groups = treatment_groups(sample, sample_ids)
    finite = values.notna().sum(axis=0)
    stats = {
        "n_probes": int(values.shape[0]),
        "n_samples": int(values.shape[1]),
        "missing_values": int(values.isna().sum().sum()),
        "sample_finite_fraction": {sid: float(finite.iloc[i] / len(values)) for i, sid in enumerate(sample_ids)},
        "sample_mean": {sid: float(values.iloc[:, i].mean()) for i, sid in enumerate(sample_ids)},
        "sample_median": {sid: float(values.iloc[:, i].median()) for i, sid in enumerate(sample_ids)},
        "sample_min": {sid: float(values.iloc[:, i].min()) for i, sid in enumerate(sample_ids)},
        "sample_max": {sid: float(values.iloc[:, i].max()) for i, sid in enumerate(sample_ids)},
        "group_counts": pd.Series(groups).value_counts().to_dict(),
    }
    return {
        "accession": metadata.get("Series_geo_accession"),
        "title": metadata.get("Series_title"),
        "platform": metadata.get("Series_platform_id"),
        "organisms": sorted(set(sample.get("Sample_organism_ch1", []))),
        "groups": groups,
        "qc": stats,
        "annotation_status": "probe_ids_unmapped_platform_annotation_not_publicly_retrieved",
    }


def main():
    report = {path.stem.split("_")[0]: audit(path) for path in sorted(MATRIX_DIR.glob("GSE*_series_matrix.txt.gz"))}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
