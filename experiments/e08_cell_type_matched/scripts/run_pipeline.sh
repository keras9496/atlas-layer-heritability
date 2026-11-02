#!/bin/bash
# 실험 3i 본 실행: i01 main → i02 dry → i02 main (14회) → i03 dry → i03 main. 순차. 상태는 PIPELINE_STATUS.txt
D=${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3i; S=$D/PIPELINE_STATUS.txt; mkdir -p $D
cd "$(dirname "$0")"
export PYTHONIOENCODING=utf-8
step() { echo "$(date '+%F %T') 시작 $*" >> $S; python -W ignore "$@" >> $D/run.log 2>&1 || { echo "$(date '+%F %T') 실패 $*" >> $S; exit 1; }; echo "$(date '+%F %T') 완료 $*" >> $S; }
step i01_matched_annot.py
step i02_sldsc.py dry
step i03_tests.py dry
step i02_sldsc.py main
step i03_tests.py main
echo "$(date '+%F %T') 전체 완료" >> $S
