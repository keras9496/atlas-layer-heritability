#!/bin/bash
# 실험 3j 본 실행 (PROTOCOL_exp3j_2026-10-14): j02 main (383회, 순차) → j03 main
set -e
cd "$(dirname "$0")"
export PYTHONIOENCODING=utf-8
S=${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3j/PIPELINE_STATUS.txt
echo "j02 main 시작 $(date '+%m-%d %H:%M')" > $S
python -W ignore j02_sldsc.py main > ../results/j02_main.log 2>&1
echo "j03 main 시작 $(date '+%m-%d %H:%M')" > $S
python -W ignore j03_tests.py main > ../results/j03_main.log 2>&1
echo "완료 $(date '+%m-%d %H:%M')" > $S
