#!/usr/bin/env bash
# 실험 3d 순차 실행기: y03 → y04 A → y05 (중간 요약) → y04 B → y05 (최종). 세션 종료와 무관하게 진행 (실험 3c run_pipeline.sh 와 같은 구조).
#  - 상태: ${AGA_DATA:-D:/AD_GWAS_Atlas_data}/exp3d/PIPELINE_STATUS.txt, 단계별 로그 exp3d/<단계>.log
#  - 실패하면 그 단계에서 멈춘다. 다시 실행하면 끝난 S-LDSC(.results)와 LD 점수는 건너뛴다.
# 사용: bash scripts/run_pipeline.sh   (e03_single_layer_map 폴더에서)
cd "$(dirname "$0")/.." || exit 1
Y=${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3d; mkdir -p $Y; ST=$Y/PIPELINE_STATUS.txt
export PYTHONIOENCODING=utf-8
st() { echo "$(date '+%F %T') $*" | tee -a "$ST"; }
step() { st "$1 시작"; if python -W ignore scripts/$2 $3 > $Y/$1.log 2>&1; then st "$1 완료"; else st "$1 실패 (로그 $Y/$1.log)"; exit 1; fi; }
st "실행기 시작"
step y03_annot_ldscores y03_annot_ldscores.py
step y04_A y04_sldsc.py A
step y05_A y05_summary.py
step y04_B y04_sldsc.py B
step y05_B y05_summary.py
st "전체 완료 → results/RESULTS_EXP3d_*.md, figures/"
