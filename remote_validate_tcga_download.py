import hashlib
import json
import pathlib
import sys

root = pathlib.Path("/root/autodl-tmp/acc-mitochondrial-mitotane")
inventory = json.loads((root / "manifests/tcga_acc_star_counts_inventory.json").read_text())["data"]["hits"]
data_dir = root / "data/raw/tcga_star_counts"
rows = []
for record in inventory:
    file_id = record["file_id"]
    name = pathlib.Path(record.get("file_name") or f"{file_id}.dat").name
    candidates = [data_dir / name, data_dir / f"{file_id}__{name}"]
    path = next((p for p in candidates if p.exists()), None)
    row = {"file_id": file_id, "file_name": name, "status": "missing"}
    if path is not None:
        digest = hashlib.md5()
        size = 0
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
                size += len(chunk)
                digest.update(chunk)
        row.update({"path": str(path), "bytes": size, "expected_bytes": int(record.get("file_size") or 0), "md5": digest.hexdigest(), "expected_md5": record.get("md5sum")})
        row["status"] = "ok" if size == int(record.get("file_size") or 0) and digest.hexdigest() == record.get("md5sum") else "mismatch"
    rows.append(row)
summary = {"n_inventory": len(rows), "n_ok": sum(x["status"] == "ok" for x in rows), "n_missing": sum(x["status"] == "missing" for x in rows), "n_mismatch": sum(x["status"] == "mismatch" for x in rows), "files": rows}
(root / "manifests/tcga_star_counts_validation.json").write_text(json.dumps(summary, indent=2))
print(json.dumps({k: summary[k] for k in ["n_inventory", "n_ok", "n_missing", "n_mismatch"]}))
if summary["n_missing"] or summary["n_mismatch"]:
    sys.exit(2)
