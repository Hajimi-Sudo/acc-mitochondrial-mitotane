#!/usr/bin/env bash
set -eu

PROJECT_ROOT=${PROJECT_ROOT:-/root/autodl-tmp/acc-mitochondrial-mitotane}
EXPECTED=/root/autodl-tmp/acc-mitochondrial-mitotane
if [ "$PROJECT_ROOT" != "$EXPECTED" ]; then
  echo 'ERROR: unexpected project root' >&2
  exit 2
fi

mkdir -p "$PROJECT_ROOT"/00_contract "$PROJECT_ROOT"/data/raw "$PROJECT_ROOT"/data/processed "$PROJECT_ROOT"/results/pilot "$PROJECT_ROOT"/logs "$PROJECT_ROOT"/code "$PROJECT_ROOT"/env "$PROJECT_ROOT"/manifests
LOG="$PROJECT_ROOT/logs/pilot_audit.log"
{
  echo '=== isolated remote pilot audit ==='
  date -u
  echo "host=$(hostname)"
  echo "user=$(id -un)"
  echo "project_root=$PROJECT_ROOT"
  echo '--- storage ---'
  df -h /root/autodl-tmp
  du -sh "$PROJECT_ROOT" /root/autodl-tmp/acc-steroid-metabolism 2>/dev/null || true
  echo '--- environment ---'
  for tool in python python3 R Rscript curl wget git; do
    printf '%s=' "$tool"
    command -v "$tool" || true
  done
  echo "conda=$(command -v /root/miniconda3/bin/conda || true)"
  /root/miniconda3/bin/conda env list 2>/dev/null || true
  echo '--- isolation ---'
  test -d /root/autodl-tmp/acc-steroid-metabolism
  test "$PROJECT_ROOT" != /root/autodl-tmp/acc-steroid-metabolism
  find "$PROJECT_ROOT" -maxdepth 2 -type d -print | sort
  echo '--- GSE76019 metadata ---'
  curl -L --fail --max-time 45 -sS 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76019&targ=self&form=text&view=quick' | sed -n '1,45p' | tee "$PROJECT_ROOT/results/pilot/GSE76019_series_metadata.txt"
  echo '--- mitotane paper metadata ---'
  curl -L --fail --max-time 45 -sS 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=38779454&retmode=xml' | tr '\n' ' ' | sed 's/<[^>]*>/ /g' | tr -s ' ' | cut -c1-12000 | tee "$PROJECT_ROOT/results/pilot/mitotane_paper_metadata.txt"
  echo 'pilot_status=metadata_audit_complete'
} 2>&1 | tee "$LOG"
