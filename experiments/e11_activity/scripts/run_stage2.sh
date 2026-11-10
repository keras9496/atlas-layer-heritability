#!/usr/bin/env bash
# 실험 3l 2단계 실행기 (프로토콜 2026-10-18 고정 후): l01 조회 완료를 기다린 뒤 l02 main → l03 dry → l04 dry → l03 main → l04 main 을 순차 실행.
#  - l01 완료 판단: PIPELINE_STATUS.txt 에 "l01 완료" 줄 (run_pipeline.sh 가 쓴다). 10분마다 확인, 최대 4일.
#  - 단계가 실패하면 멈춘다(로그 exp3l/<단계>.log). 고친 뒤 다시 실행하면 끝난 실행(npz, LD 점수)은 건너뛴다.
#  - 상태: ${AGA_DATA:-D:/AD_GWAS_Atlas_data}/exp3l/PIPELINE_STATUS.txt
# 사용: bash scripts/run_stage2.sh   (e11_activity 폴더에서)
cd "$(dirname "$0")/.." || exit 1
Y=${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3l; ST=$Y/PIPELINE_STATUS.txt
export PYTHONIOENCODING=utf-8
st() { echo "$(date '+%F %T') $*" | tee -a "$ST"; }
st "2단계 실행기 시작 — 조회 끝날 때까지 기다림"
done_l01() { grep -qE "^[0-9-]+ [0-9:]+ l01 완료 \(parts" "$ST" && [ "$(ls $Y/atlas_active/part_*.parquet 2>/dev/null | wc -l)" -eq 583 ]; }   # run_pipeline.sh 의 줄 형식 + 조각 583개 (2026-10-18 수정: 상태 문구 오인)
w=0; while ! done_l01; do w=$((w + 1)); [ $w -gt 576 ] && { st "l01 4일 대기 초과 → 중단"; exit 1; }; sleep 600; done
for step in "l02_annot_ldscores main" "l03_sldsc dry" "l04_tests dry" "l03_sldsc main" "l04_tests main"; do
  set -- $step; name="$1_$2"; st "$name 시작"
  if python -W ignore scripts/$1.py $2 > $Y/$name.log 2>&1; then st "$name 완료"; else st "$name 실패 (로그 $Y/$name.log)"; exit 1; fi
done
st "3l 전체 완료 → results/l04_main_tests.txt"
