#!/usr/bin/env python3
"""Remote-only TCGA-ACC inventory and resumable downloader."""
import argparse
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import requests

ROOT = Path("/root/autodl-tmp/acc-mitochondrial-mitotane")
RAW = ROOT / "data" / "raw"
MANIFESTS = ROOT / "manifests"
LOGS = ROOT / "logs"
GDC = "https://api.gdc.cancer.gov"


def session():
    s = requests.Session()
    s.trust_env = False
    s.headers.update({"User-Agent": "ACC-mitochondrial-mitotane/0.1"})
    return s


def query(endpoint, payload):
    s = session()
    for attempt in range(1, 5):
        try:
            r = s.post(f"{GDC}/{endpoint}", json=payload, timeout=(10, 60))
            r.raise_for_status()
            return r.json()
        except Exception:
            if attempt == 4:
                raise
            time.sleep(attempt)


def file_filter():
    return {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.project.project_id", "value": "TCGA-ACC"}},
            {"op": "=", "content": {"field": "data_type", "value": "Gene Expression Quantification"}},
            {"op": "=", "content": {"field": "access", "value": "open"}},
            {"op": "=", "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"}},
        ],
    }


def inventory():
    MANIFESTS.mkdir(parents=True, exist_ok=True)
    fields = ",".join([
        "file_id", "file_name", "file_size", "md5sum", "data_format", "data_type",
        "analysis.workflow_type", "access", "cases.case_id", "cases.submitter_id",
        "cases.samples.submitter_id", "cases.samples.sample_type",
    ])
    file_payload = {"filters": file_filter(), "fields": fields, "format": "JSON", "size": "10000"}
    case_payload = {
        "filters": {"op": "=", "content": {"field": "project.project_id", "value": "TCGA-ACC"}},
        "fields": ",".join(["case_id", "submitter_id", "project.project_id", "demographic", "diagnoses", "samples.submitter_id", "samples.sample_type"]),
        "format": "JSON", "size": "10000",
    }
    files = query("files", file_payload)
    cases = query("cases", case_payload)
    (MANIFESTS / "tcga_acc_star_counts_inventory.json").write_text(json.dumps(files, indent=2), encoding="utf-8")
    (MANIFESTS / "tcga_acc_cases_inventory.json").write_text(json.dumps(cases, indent=2), encoding="utf-8")
    return files.get("data", {}).get("hits", [])


def download_one(rec):
    file_id = rec["file_id"]
    name = Path(rec.get("file_name") or f"{file_id}.dat").name
    dest = RAW / "tcga_star_counts" / f"{file_id}__{name}"
    dest.parent.mkdir(parents=True, exist_ok=True)
    expected = int(rec.get("file_size") or 0)
    last = ""
    for attempt in range(1, 5):
        try:
            existing = dest.stat().st_size if dest.exists() else 0
            if expected and existing > expected:
                dest.unlink()
                existing = 0
            headers = {"Range": f"bytes={existing}-"} if existing and expected and existing < expected else {}
            with session().get(f"{GDC}/data/{file_id}", headers=headers, stream=True, timeout=(15, 300)) as r:
                if r.status_code == 416 and expected and existing == expected:
                    pass
                else:
                    r.raise_for_status()
                    mode = "ab" if existing and r.status_code == 206 else "wb"
                    with dest.open(mode) as out:
                        for chunk in r.iter_content(chunk_size=4 * 1024 * 1024):
                            if chunk:
                                out.write(chunk)
            size = dest.stat().st_size
            if expected and size != expected:
                raise IOError(f"size mismatch {size} != {expected}")
            sha = hashlib.sha256()
            with dest.open("rb") as fh:
                for chunk in iter(lambda: fh.read(4 * 1024 * 1024), b""):
                    sha.update(chunk)
            return {"file_id": file_id, "file_name": name, "path": str(dest), "bytes": size, "sha256": sha.hexdigest(), "status": "ok", "attempts": attempt}
        except Exception as exc:
            last = repr(exc)
            time.sleep(min(2 ** attempt, 16))
    return {"file_id": file_id, "file_name": name, "path": str(dest), "status": "failed", "error": last}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--inventory-only", action="store_true")
    args = ap.parse_args()
    ROOT.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    recs = inventory()
    print(f"open TCGA-ACC STAR count files: {len(recs)}", flush=True)
    if args.inventory_only:
        return
    results = []
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = [pool.submit(download_one, rec) for rec in recs]
        for i, future in enumerate(as_completed(futures), 1):
            result = future.result()
            results.append(result)
            print(f"[{i}/{len(futures)}] {result['status']} {result['file_name']} {result.get('bytes', 0)}", flush=True)
    manifest = {"started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(), "workers": args.workers, "n_requested": len(results), "n_ok": sum(x["status"] == "ok" for x in results), "n_failed": sum(x["status"] != "ok" for x in results), "files": sorted(results, key=lambda x: x["file_id"])}
    (MANIFESTS / "tcga_star_counts_download.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    if manifest["n_failed"]:
        sys.exit(2)


if __name__ == "__main__":
    main()
