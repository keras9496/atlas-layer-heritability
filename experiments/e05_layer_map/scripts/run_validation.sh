#!/usr/bin/env bash
# 실험 3f 3단계(검증) 순차 실행기: v02 val (99회) → v03 val. 예측은 프로토콜 "3단계 예측" 표에 실행 전 고정했다.
# 사용: bash scripts/run_validation.sh   (e05_layer_map 폴더에서)
cd "$(dirname "$0")/.." || exit 1
F=${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3f; ST=$F/PIPELINE_STATUS.txt
export PYTHONIOENCODING=utf-8
st() { echo "$(date '+%F %T') $*" | tee -a "$ST"; }
step() { st "$1 시작"; if python -W ignore scripts/$2 $3 > $F/$1.log 2>&1; then st "$1 완료"; else st "$1 실패 (로그 $F/$1.log)"; exit 1; fi; }
st "실행기 시작 (3단계 검증)"
step v02_val v02_sldsc.py val
step v03_val v03_tests.py val
st "3단계 완료 → results/v03_val_tests.txt"
