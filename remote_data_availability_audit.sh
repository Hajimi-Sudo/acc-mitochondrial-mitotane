#!/usr/bin/env bash
set -eu

PROJECT_ROOT=${PROJECT_ROOT:-/root/autodl-tmp/acc-mitochondrial-mitotane}
PYTHON_BIN=${PYTHON_BIN:-python3}
OUT="$PROJECT_ROOT/results/pilot"
mkdir -p "$OUT"
LOG="$PROJECT_ROOT/logs/data_availability_audit.log"
echo '=== remote data availability audit ===' | tee "$LOG"
date -u | tee -a "$LOG"

NCBI_CURL="env -u HTTP_PROXY -u HTTPS_PROXY -u ALL_PROXY -u http_proxy -u https_proxy -u all_proxy -u DOWNLOAD_PROXY curl --noproxy '*' -L --fail --connect-timeout 8 --max-time 30 -sS"

echo '--- GSE76019 series/sample annotations ---' | tee "$OUT/GSE76019_annotation_audit.txt"
eval "$NCBI_CURL 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76019&targ=self&form=text&view=quick'" | awk '/^\^SERIES|^!Series_title|^!Series_overall_design|^!Series_status|^!Series_pubmed_id|^!Series_sample_id/ {print}' | tee -a "$OUT/GSE76019_annotation_audit.txt"
PARTS="$OUT/GSE76019_sample_parts"
mkdir -p "$PARTS"
"$PYTHON_BIN" - "$PARTS" <<'PY'
import concurrent.futures, pathlib, sys, urllib.request
out = pathlib.Path(sys.argv[1])
prefix = ('^SAMPLE', '!Sample_title', '!Sample_source_name_ch1',
          '!Sample_characteristics_ch1', '!Sample_description')
def fetch(n):
    gsm = f'GSM{n}'
    url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=quick'
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'ACC-mito-pilot/1.0'})
        with urllib.request.urlopen(req, timeout=20) as h:
            text = h.read().decode('utf-8', 'replace')
        lines = [x for x in text.splitlines() if x.startswith(prefix)]
        (out / f'{gsm}.txt').write_text('\n'.join(lines) + '\n')
        return gsm, len(lines), ''
    except Exception as exc:
        return gsm, 0, str(exc)
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
    for gsm, count, error in pool.map(fetch, range(1972883, 1972917)):
        print(f'{gsm}\t{count}\t{error}')
PY
cat "$PARTS"/GSM*.txt 2>/dev/null | tee -a "$OUT/GSE76019_annotation_audit.txt" || true

echo '--- GSE76019 treatment/response keyword scan ---' | tee "$OUT/GSE76019_keyword_scan.txt"
grep -i -E 'mitotane|treatment|response|responder|non-responder|hormone|secretion|cortisol|androgen|steroid' "$OUT/GSE76019_annotation_audit.txt" | tee -a "$OUT/GSE76019_keyword_scan.txt" || true

echo '--- candidate paper supplementary file list ---' | tee "$OUT/mitotane_paper_data_availability.txt"
eval "$NCBI_CURL 'https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC11109426'" | tee -a "$OUT/mitotane_paper_data_availability.txt" || true

echo '--- TCGA-ACC project metadata only ---' | tee "$OUT/TCGA_ACC_project_metadata.json"
env -u HTTP_PROXY -u HTTPS_PROXY -u ALL_PROXY -u http_proxy -u https_proxy -u all_proxy -u DOWNLOAD_PROXY curl --noproxy '*' -L --fail --connect-timeout 8 --max-time 20 -sS 'https://api.gdc.cancer.gov/projects/TCGA-ACC?pretty=true' | tee -a "$OUT/TCGA_ACC_project_metadata.json" || echo 'GDC metadata request failed; retry through proxy or alternate endpoint.' | tee -a "$OUT/TCGA_ACC_project_metadata.json"

echo 'audit_status=metadata_and_endpoint_availability_checked' | tee -a "$LOG"
