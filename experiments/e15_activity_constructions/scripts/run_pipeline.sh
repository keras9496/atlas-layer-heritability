#!/usr/bin/env bash
# 실험 3p 순차 실행: p01 check가 통과한 뒤에 돌린다. 단계마다 실패하면 멈춘다.
set -e
cd "$(dirname "$0")"
D=${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3p
export PYTHONIOENCODING=utf-8
echo "p01 main 시작 $(date)"      >  $D/pipeline.log
python -u p01_annot_ldscores.py main >> $D/pipeline.log 2>&1
echo "p02 dry 시작 $(date)"       >> $D/pipeline.log
python -u p02_sldsc.py dry          >> $D/pipeline.log 2>&1
python -u p03_tests.py dry          >> $D/pipeline.log 2>&1
echo "p02 main 시작 $(date)"      >> $D/pipeline.log
python -u p02_sldsc.py main         >> $D/pipeline.log 2>&1
echo "p03 main 시작 $(date)"      >> $D/pipeline.log
python -u p03_tests.py main         >> $D/pipeline.log 2>&1
python -u p04_fig_export.py         >> $D/pipeline.log 2>&1
echo "끝 $(date)"                  >> $D/pipeline.log
