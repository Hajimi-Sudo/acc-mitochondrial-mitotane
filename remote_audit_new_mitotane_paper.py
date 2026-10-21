#!/usr/bin/env python3
"""Audit the newly indexed mitotane-response paper for public matched data."""

import json
import pathlib
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET


ROOT = pathlib.Path("/root/autodl-tmp/acc-mitochondrial-mitotane")
OUT = ROOT / "results/literature_verification/drug_sensitivity_search"
OUT.mkdir(parents=True, exist_ok=True)


def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": "ACC-mitochondrial-mitotane/1.0"})
    with urllib.request.urlopen(req, timeout=60) as handle:
        return handle.read()


pmid = "42712189"
europe_url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?{urllib.parse.urlencode({'query': 'EXT_ID:' + pmid, 'format': 'json', 'resultType': 'core'})}"
core = json.loads(fetch(europe_url))
full_url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/MED/{pmid}/fullTextXML"
full_text = ""
try:
    full_text = fetch(full_url).decode("utf-8", errors="replace")
except Exception as exc:
    full_text = f"FULLTEXT_ERROR: {exc}"

if full_text.startswith("<"):
    try:
        root = ET.fromstring(full_text)
        plain = " ".join("".join(root.itertext()).split())
    except Exception:
        plain = " ".join(full_text.split())
else:
    plain = full_text
keywords = [
    "GEO", "ArrayExpress", "PRIDE", "ProteomeXchange", "Zenodo", "GitHub",
    "data availability", "supplementary", "RNA-seq", "transcriptomic", "proteomic",
    "IC50", "AUC", "responder", "non-responder", "cell line", "sample",
]
snippets = {}
for keyword in keywords:
    low = plain.lower()
    idx = low.find(keyword.lower())
    if idx >= 0:
        snippets[keyword] = plain[max(0, idx - 500): idx + 1800]

payload = {
    "pmid": pmid,
    "europepmc_url": europe_url,
    "fulltext_url": full_url,
    "record": core.get("resultList", {}).get("result", []),
    "fulltext_bytes": len(full_text.encode("utf-8")),
    "keyword_snippets": snippets,
}
(OUT / "PMID_42712189_audit.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps(payload, indent=2, ensure_ascii=False)[:30000])
