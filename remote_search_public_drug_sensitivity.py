#!/usr/bin/env python3
"""Search public NCBI/Europe PMC records for sample-matched mitotane data."""

import json
import pathlib
import urllib.parse
import urllib.request


ROOT = pathlib.Path("/root/autodl-tmp/acc-mitochondrial-mitotane")
OUT = ROOT / "results/literature_verification/drug_sensitivity_search"
OUT.mkdir(parents=True, exist_ok=True)

queries = {
    "geo_mitotane": ("geo", "mitotane[All Fields]"),
    "geo_adrenocortical_mitotane": ("geo", "(mitotane[All Fields]) AND (adrenocortical[All Fields] OR adrenal[All Fields])"),
    "pubmed_mitotane_acc": ("pubmed", "mitotane AND adrenocortical carcinoma"),
    "pubmed_mitotane_ic50": ("pubmed", "mitotane AND (IC50 OR AUC OR sensitivity OR resistant)"),
}
results = {}
for name, (db, term) in queries.items():
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode({"db": db, "term": term, "retmode": "json", "retmax": 100})
    with urllib.request.urlopen(url, timeout=60) as handle:
        results[name] = json.load(handle)

for name, term in [
    ("europepmc_mitotane", "mitotane"),
    ("europepmc_mitotane_acc", "mitotane AND (adrenocortical OR adrenal carcinoma)"),
    ("europepmc_mitotane_ic50", "mitotane AND (IC50 OR AUC OR sensitivity)"),
]:
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": term, "format": "json", "pageSize": 100})
    with urllib.request.urlopen(url, timeout=60) as handle:
        results[name] = json.load(handle)

(OUT / "structured_public_drug_sensitivity_search.json").write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
for name, data in results.items():
    if "esearchresult" in data:
        print(name, data["esearchresult"].get("count"), data["esearchresult"].get("idlist", [])[:20])
    else:
        rows = data.get("resultList", {}).get("result", [])
        print(name, data.get("hitCount"), [(row.get("pmid"), row.get("title")) for row in rows[:10]])
