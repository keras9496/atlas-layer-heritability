#!/usr/bin/env bash
# 실험 3r: 2026 AD 합의 메타분석 proxy 제외 요약통계(GWAS Catalog GCST90704647, GRCh38) 내려받기와 md5 확인
set -e
D=${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3r/raw
B=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90704001-GCST90705000/GCST90704647
cd $D
curl -s -o GCST90704647.tsv.gz-meta.yaml $B/GCST90704647.tsv.gz-meta.yaml
curl -s -o md5sum.txt $B/md5sum.txt
curl -s -C - -o GCST90704647.tsv.gz $B/GCST90704647.tsv.gz
md5sum GCST90704647.tsv.gz
cat md5sum.txt
