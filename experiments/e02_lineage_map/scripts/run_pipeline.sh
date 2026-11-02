#!/usr/bin/env bash
# 실험 3c 순차 실행기: z02(조회) 완료 대기 → z03 → z04 → z05. 대화 compact·세션 종료와 무관하게 진행 (실험 3 run_pipeline.sh 와 같은 구조).
#  - z02 가 완료 전에 죽으면 다시 시작한다(재시작 가능, 최대 30회). z02 는 끝에 os._exit(0) 로 스스로 종료한다.
#  - z03 → z04 → z05 는 순차. z04/z05 가 아직 없으면 그 단계에서 기다린다(파일이 생길 때까지 10분마다 확인, 최대 3일).
#  - 상태: ${AGA_DATA:-D:/AD_GWAS_Atlas_data}/exp3c/PIPELINE_STATUS.txt, 단계별 로그 exp3c/<단계>.log
# 사용: bash scripts/run_pipeline.sh   (e02_lineage_map 폴더에서)
cd "$(dirname "$0")/.." || exit 1
Y=${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3c; ST=$Y/PIPELINE_STATUS.txt; NPART=583
export PYTHONIOENCODING=utf-8
st() { echo "$(date '+%F %T') $*" | tee -a "$ST"; }
nparts() { ls $Y/atlas/part_*.parquet 2>/dev/null | wc -l; }
st "실행기 시작 (parts $(nparts)/$NPART)"
tries=0
while [ "$(nparts)" -lt $NPART ]; do
  tries=$((tries + 1)); [ $tries -gt 30 ] && { st "z02 재시작 30회 초과 → 중단"; exit 1; }
  st "z02 실행 (parts $(nparts), 시도 $tries)"
  python -W ignore scripts/z02_atlas_query.py >> $Y/z02.log 2>&1
done
st "z02 완료 (parts $(nparts))"
for s in z03_annot_ldscores z04_sldsc z05_summary; do
  w=0; while [ ! -f scripts/$s.py ]; do [ $w -eq 0 ] && st "$s.py 없음 → 대기"; w=$((w + 1)); [ $w -gt 432 ] && { st "$s.py 3일 대기 초과 → 중단"; exit 1; }; sleep 600; done
  st "$s 시작"
  if python -W ignore scripts/$s.py > $Y/$s.log 2>&1; then st "$s 완료"; else st "$s 실패 (로그 $Y/$s.log)"; exit 1; fi
done
st "전체 완료 → results/RESULTS_EXP3c_*.md, figures/"
