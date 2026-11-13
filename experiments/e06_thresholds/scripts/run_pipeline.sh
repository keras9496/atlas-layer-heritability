#!/usr/bin/env bash
# 실험 3g 순차 실행기: u01 → u02 시험(Z 섞은 가짜 형질, 결과 아님) → u03 시험 → u02 본 실행 → u03 본 검정 (3f run_pipeline.sh 와 같은 구조).
#  - 상태: ${AGA_DATA:-D:/AD_GWAS_Atlas_data}/exp3g/PIPELINE_STATUS.txt, 단계별 로그 exp3g/<단계>.log
#  - 실패하면 그 단계에서 멈춘다. 다시 실행하면 끝난 S-LDSC(.results)와 LD 점수는 건너뛴다.
# 사용: bash scripts/run_pipeline.sh   (e06_thresholds 폴더에서)
cd "$(dirname "$0")/.." || exit 1
G=${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3g; mkdir -p $G; ST=$G/PIPELINE_STATUS.txt
export PYTHONIOENCODING=utf-8
st() { echo "$(date '+%F %T') $*" | tee -a "$ST"; }
step() { st "$1 시작"; if python -W ignore scripts/$2 $3 > $G/$1.log 2>&1; then st "$1 완료"; else st "$1 실패 (로그 $G/$1.log)"; exit 1; fi; }
st "실행기 시작"
step u01_threshold_annot u01_threshold_annot.py
step u02_dry u02_sldsc.py dry
step u03_dry u03_tests.py dry
step u02_main u02_sldsc.py main
step u03_main u03_tests.py main
st "완료 → results/u03_main_tests.txt, figures/exp3g_threshold.png"
