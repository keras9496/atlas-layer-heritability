#!/usr/bin/env bash
# 실험 3f 1단계 순차 실행기: v01 → v02 시험(Z 섞은 가짜 형질, 결과 아님) → v03 시험 → v02 본 실행 → v03 본 검정 (3e run_pipeline.sh 와 같은 구조).
#  - 상태: ${AGA_DATA:-D:/AD_GWAS_Atlas_data}/exp3f/PIPELINE_STATUS.txt, 단계별 로그 exp3f/<단계>.log
#  - 실패하면 그 단계에서 멈춘다. 다시 실행하면 끝난 S-LDSC(.results)와 LD 점수는 건너뛴다.
#  - 3단계(검증)는 1단계 분류로 예측을 프로토콜에 적은 뒤 따로 실행한다.
# 사용: bash scripts/run_pipeline.sh   (e05_layer_map 폴더에서)
cd "$(dirname "$0")/.." || exit 1
F=${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3f; mkdir -p $F; ST=$F/PIPELINE_STATUS.txt
export PYTHONIOENCODING=utf-8
st() { echo "$(date '+%F %T') $*" | tee -a "$ST"; }
step() { st "$1 시작"; if python -W ignore scripts/$2 $3 > $F/$1.log 2>&1; then st "$1 완료"; else st "$1 실패 (로그 $F/$1.log)"; exit 1; fi; }
st "실행기 시작 (1단계)"
step v01_resample_annot v01_resample_annot.py
step v02_dry v02_sldsc.py dry
step v03_dry v03_tests.py dry
step v02_main v02_sldsc.py main
step v03_main v03_tests.py main
st "1단계 완료 → results/v03_main_tests.txt, figures/exp3f_layer_map.png"
