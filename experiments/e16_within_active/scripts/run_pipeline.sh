#!/usr/bin/env bash
# 실험 3q 순차 실행: q01 check가 통과한 뒤에 돌린다. 단계마다 실패하면 멈춘다.
set -e
cd "$(dirname "$0")"
D=${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3q
export PYTHONIOENCODING=utf-8
echo "q01 main 시작 $(date)"      >  $D/pipeline.log
python -u q01_annot_ldscores.py main >> $D/pipeline.log 2>&1
echo "q02 dry 시작 $(date)"       >> $D/pipeline.log
python -u q02_sldsc.py dry          >> $D/pipeline.log 2>&1
python -u q03_tests.py dry          >> $D/pipeline.log 2>&1
echo "q02 main 시작 $(date)"      >> $D/pipeline.log
python -u q02_sldsc.py main         >> $D/pipeline.log 2>&1
echo "q03 main 시작 $(date)"      >> $D/pipeline.log
python -u q03_tests.py main         >> $D/pipeline.log 2>&1
echo "끝 $(date)"                  >> $D/pipeline.log
