#!/usr/bin/env bash
# 실험 3n n00: 최신 GWAS 요약통계 내려받기 (2026-10-21). 공개 파일만, 이어받기 가능. 분석 전 형식 확인용.
R=${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3n/raw; L=$R/download.log; cd $R || exit 1
get() { echo "$(date '+%F %T') 시작 $1" >> $L; curl -L --fail --retry 5 --retry-delay 10 -C - -o "$1" "$2" >> $L 2>&1 && echo "$(date '+%F %T') 완료 $1 $(stat -c %s "$1") bytes" >> $L || echo "$(date '+%F %T') 실패 $1" >> $L; }
G=https://giant-consortium.web.broadinstitute.org/images
get README_Yengo_2022.txt "$G/e/ee/Readme_Yengo_et_al_2022_Nature.txt"
get README_Yengo_2018.txt "$G/0/01/README_summary_statistics_Yengo_et_al_2018.txt"
get README_bip2024.txt https://ndownloader.figshare.com/files/49760763
get README_PGC3_SCZ_v3.2.pdf https://ndownloader.figshare.com/files/34865091
get README_GLGC2021.txt https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/README.txt
get SCZ_PGC3_wave3_european.vcf.tsv.gz https://ndownloader.figshare.com/files/34517828
get BIP_bip2024_eur_no23andMe.gz https://ndownloader.figshare.com/files/49760772
get BMI_Yengo2018_UPDATED.txt.gz "$G/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz"
get HEIGHT_Yengo2022_EUR.gz "$G/f/f7/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz"
get LDL_GLGC2021_EUR.gz https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz
echo "$(date '+%F %T') 전체 끝" >> $L
