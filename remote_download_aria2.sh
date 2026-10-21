#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/autodl-tmp/acc-mitochondrial-mitotane
INV=$ROOT/manifests/tcga_acc_star_counts_inventory.json
OUT=$ROOT/data/raw/tcga_star_counts
LIST=$ROOT/manifests/tcga_star_counts_aria2.list
LOG=$ROOT/logs/tcga_star_counts_aria2.log
mkdir -p "$OUT" "$ROOT/manifests"

/root/miniconda3/bin/python "$ROOT/code/make_aria2_list.py" "$INV" "$LIST"

env -u HTTP_PROXY -u HTTPS_PROXY -u ALL_PROXY -u http_proxy -u https_proxy -u all_proxy aria2c --no-proxy='*' --dir="$OUT" --input-file="$LIST" --continue=true --allow-overwrite=false --auto-file-renaming=false --file-allocation=none --max-concurrent-downloads=8 --max-connection-per-server=4 --split=4 --min-split-size=1M --connect-timeout=10 --timeout=90 --max-tries=4 --retry-wait=2 --summary-interval=10 --console-log-level=notice 2>&1 | tee "$LOG"

echo 'aria2_status=download_command_complete' | tee -a "$LOG"
