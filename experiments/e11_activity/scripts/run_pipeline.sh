#!/usr/bin/env bash
# 실험 3l 실행기 (1단계: l01 Atlas 활성도 조회만). 분석 단계(l02–)는 따로 실행한다.
#  - l01 이 완료 전에 죽으면 다시 시작한다(재시작 가능, 최대 30회). l01 은 끝에 os._exit(0) 로 스스로 종료한다.
#  - 상태: ${AGA_DATA:-D:/AD_GWAS_Atlas_data}/exp3l/PIPELINE_STATUS.txt, 로그 exp3l/l01.log
# 사용: bash scripts/run_pipeline.sh   (e11_activity 폴더에서)
cd "$(dirname "$0")/.." || exit 1
Y=${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3l; ST=$Y/PIPELINE_STATUS.txt; NPART=583
export PYTHONIOENCODING=utf-8
st() { echo "$(date '+%F %T') $*" | tee -a "$ST"; }
nparts() { ls $Y/atlas_active/part_*.parquet 2>/dev/null | wc -l; }
st "실행기 시작 (parts $(nparts)/$NPART)"
tries=0
while [ "$(nparts)" -lt $NPART ]; do
  tries=$((tries + 1)); [ $tries -gt 30 ] && { st "l01 재시작 30회 초과 → 중단"; exit 1; }
  st "l01 실행 (parts $(nparts), 시도 $tries)"
  python -W ignore scripts/l01_atlas_active.py >> $Y/l01.log 2>&1
done
st "l01 완료 (parts $(nparts)) — 조회 끝. 분석은 프로토콜 고정 후"
