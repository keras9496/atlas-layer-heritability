#!/bin/bash
# 실험 3h 본 실행: h01 main (33회, 순차) → h02 main. 상태는 PIPELINE_STATUS.txt
S=${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3h/PIPELINE_STATUS.txt
cd "$(dirname "$0")"
export PYTHONIOENCODING=utf-8
echo "$(date '+%F %T') 시작 h01 main" >> $S
python h01_sldsc.py main >> ${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3h/h01_main.log 2>&1 || { echo "$(date '+%F %T') h01 실패" >> $S; exit 1; }
echo "$(date '+%F %T') h01 완료 → h02 main" >> $S
python h02_tests.py main >> ${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3h/h02_main.log 2>&1 || { echo "$(date '+%F %T') h02 실패" >> $S; exit 1; }
echo "$(date '+%F %T') 전체 완료" >> $S
