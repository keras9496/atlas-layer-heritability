#!/usr/bin/env bash
# 실험 3 순차 실행기: x02(조회) 완료 대기 → x03 → x04 → x05.
#  - 실행 중인 x02 가 있으면 끝날 때까지 기다린다. x02 가 완료 전에 죽으면 다시 시작한다(재시작 가능, 최대 20회).
#  - 각 단계 상태를 ${AGA_DATA:-D:/AD_GWAS_Atlas_data}/exp3/PIPELINE_STATUS.txt 에 남긴다. 실패하면 그 단계에서 멈춘다.
# 사용: bash scripts/run_pipeline.sh   (e01_targets_sldsc 폴더에서)
cd "$(dirname "$0")/.." || exit 1
X=${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3; ST=$X/PIPELINE_STATUS.txt; NPART=583          # 5,828,323 / 10,000 → 583 parts
export PYTHONIOENCODING=utf-8
st() { echo "$(date '+%F %T') $*" | tee -a "$ST"; }
x02_running() { [ "$(powershell -NoProfile -Command "@(Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | Where-Object { \$_.CommandLine -like '*x02_atlas_query*' }).Count" | tr -d '\r')" -gt 0 ]; }
nparts() { ls $X/atlas/part_*.parquet 2>/dev/null | wc -l; }
st "실행기 시작 (parts $(nparts)/$NPART)"
tries=0
while [ "$(nparts)" -lt $NPART ]; do
  if x02_running; then sleep 600; continue; fi
  tries=$((tries + 1)); [ $tries -gt 20 ] && { st "x02 재시작 20회 초과 → 중단"; exit 1; }
  st "x02 실행 중 아님 (parts $(nparts)) → 재시작 $tries"
  python scripts/x02_atlas_query.py >> $X/x02.log 2>&1
done
for i in $(seq 10); do x02_running || break; sleep 60; done       # s03b 처럼 끝난 뒤 종료가 멈추는 경우: 10분 뒤 강제 종료
x02_running && { st "x02 가 완료 후 종료되지 않음 → 프로세스 종료"; powershell -NoProfile -Command "Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | Where-Object { \$_.CommandLine -like '*x02_atlas_query*' } | ForEach-Object { Stop-Process -Id \$_.ProcessId -Force }"; }
st "x02 완료 (parts $(nparts))"
for s in x03_annot_ldscores x04_sldsc x05_summary; do
  st "$s 시작"
  if python -W ignore scripts/$s.py > $X/$s.log 2>&1; then st "$s 완료"; else st "$s 실패 (로그 $X/$s.log)"; exit 1; fi
done
st "전체 완료 → results/RESULTS_EXP3_*.md, figures/exp3_sldsc.png"
