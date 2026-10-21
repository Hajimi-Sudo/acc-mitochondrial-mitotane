import json
import pathlib
import sys

inventory = json.load(open(sys.argv[1], encoding="utf-8"))["data"]["hits"]
output = pathlib.Path(sys.argv[2])
lines = []
for record in inventory:
    file_id = record["file_id"]
    name = pathlib.Path(record.get("file_name") or f"{file_id}.dat").name
    lines.append(f"https://api.gdc.cancer.gov/data/{file_id}")
output.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"aria2_items={len(lines)}")
