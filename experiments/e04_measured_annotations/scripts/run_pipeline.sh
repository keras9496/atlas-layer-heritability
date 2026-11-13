#!/usr/bin/env bash
# 실험 3e 순차 실행기: w02 → w03 시험(Z 섞은 가짜 형질, 결과 아님) → w03 본 실행 → w04 (실험 3d run_pipeline.sh 와 같은 구조).
#  - 상태: ${AGA_DATA:-D:/AD_GWAS_Atlas_data}/exp3e/PIPELINE_STATUS.txt, 단계별 로그 exp3e/<단계>.log
#  - 실패하면 그 단계에서 멈춘다. 다시 실행하면 끝난 S-LDSC(.results)와 LD 점수는 건너뛴다.
#  - w01(자료 받기)은 따로 먼저 실행한다.
# 사용: bash scripts/run_pipeline.sh   (e04_measured_annotations 폴더에서)
cd "$(dirname "$0")/.." || exit 1
W=${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3e; mkdir -p $W; ST=$W/PIPELINE_STATUS.txt
export PYTHONIOENCODING=utf-8
st() { echo "$(date '+%F %T') $*" | tee -a "$ST"; }
step() { st "$1 시작"; if env $3 python -W ignore scripts/$2 > $W/$1.log 2>&1; then st "$1 완료"; else st "$1 실패 (로그 $W/$1.log)"; exit 1; fi; }
st "실행기 시작"
step w02_annot_ldscores w02_annot_ldscores.py
step w03_dryrun w03_sldsc.py W03_DRYRUN=1
step w03_sldsc w03_sldsc.py W03_DRYRUN=0
step w04_summary w04_summary.py
st "전체 완료 → results/w04_summary.txt, figures/exp3e_conditional.png"
