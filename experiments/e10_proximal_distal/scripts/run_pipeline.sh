#!/bin/bash
# 실험 3k 본 실행 (PROTOCOL_exp3k_2026-10-15): k02 main (14회) → k03 main
set -e
cd "$(dirname "$0")"
export PYTHONIOENCODING=utf-8
python -W ignore k02_sldsc.py main > ../results/k02_main.log 2>&1
python -W ignore k03_tests.py main > ../results/k03_main.log 2>&1
echo "완료 $(date '+%m-%d %H:%M')" > ${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3k/PIPELINE_STATUS.txt
