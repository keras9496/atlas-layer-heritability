#!/usr/bin/env bash
# 실험 3s: ADSP FunGen-xQTL Atlas (NIAGADS NG00184.v1, 공개, GRCh38) 단일 맥락 정밀지도화 결과 내려받기
# 공개 포털 경로(이메일·로그인 없음): https://st1.niagads.org/portal/v1/download-public/NG00184/fileset/<file_ID>
#   57059 = eQTL.single_context_finemapping_all.tar (bulk·정제 미세아교·단핵구 eQTL, 약 365 MB)
#   57109 = snuc-eQTL.single_context_finemapping_all.tar (ROSMAP 단일핵 6 세포형, 약 314 MB)
#   57112/57113/57114 = 메타데이터(json, txt), readme
set -e
D=${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}/exp3s/raw
U=https://st1.niagads.org/portal/v1/download-public/NG00184/fileset
cd $D
for id in 57113 57114 57058 57108 57059 57109; do
  curl -s -L -OJ "$U/$id" || echo "실패 $id"
done
ls -la
md5sum ADSP_FunGen_xQTL.v1.eQTL.single_context_finemapping_all.tar ADSP_FunGen_xQTL.v1.snuc-eQTL.single_context_finemapping_all.tar
