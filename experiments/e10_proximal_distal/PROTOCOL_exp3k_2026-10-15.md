# 실험 3k: 근위·원위 분리 — CAGE형 프로필이 프로모터 근접만으로 설명되는가 (2026-10-15, 결과 보기 전 고정)

## 배경
- 공동저자 피드백(2026-10-14) 7번 셋째 항목: CAGE 우위의 성격을 확인하라(연속형 점수 또는 프로모터 근처·원위 분리). 사용자 결정(2026-10-14): 근위·원위 분리를 한다.
- 원고 한계 넷째: "CAGE predictions may partly reflect proximity to promoters … we did not separate proximal from distal CAGE effects." 서열 모델은 원위 인핸서보다 프로모터의 발현 결정 요인을 더 잘 잡는다(Karollus 2023).
- 질문: TSS에서 먼 SNP만 보아도 형질 사이 층위 차이가 남고 CAGE형 형질에서 CAGE가 1위인가? 남으면 CAGE형은 프로모터 근접만의 결과가 아니다(인핸서 전사 시작 등). 남지 않으면 CAGE형은 "프로모터 근처의 예측 효과"로 좁혀 해석해야 한다.

## 사후성 명시
- 이미 본 것: 3c–3j의 모든 결과. 기준선 모형에는 이미 프로모터·TSS 주석(Promoter_UCSC, TSS_Hoffman 등)이 있다.
- 아직 보지 않은 것: 근위·원위로 나눈 주석, 그 크기, 그 S-LDSC 결과.

## 주석 (k01)
- TSS: GENCODE v46 종합 주석(GRCh38)에서 단백질 암호화·lncRNA 유전자의 모든 전사체 시작 위치(+ 가닥은 start, − 가닥은 end), 1–22번 염색체.
  - 모든 전사체를 쓰므로 대체 프로모터도 근위에 들어간다. 그래서 원위는 엄격한 정의다.
- 근위 = 가장 가까운 TSS까지 1,000 bp 이하, 원위 = 1,000 bp 초과. 대상 SNP 5,828,323개(GRCh38 위치)에서 계산한다.
- 주 지도와 같은 전체 트랙 상위 5 % 주석(3d `L5_{g}_{m}`)을 그대로 나눈다: `P_{g}_{m}` = L5 ∧ 근위, `D_{g}_{m}` = L5 ∧ 원위. 다시 순위를 매기지 않는다(주 결과를 분해하는 것이 목적이다).
  - 계통 g = 각 형질의 자기 계통 5개(골수계, 뇌 조직, 림프계, 간, 결합조직), 층위 m = CAGE, DNase, H3K27ac → 30개 주석.
- 확인(어긋나면 멈춘다): (1) 모든 SNP에서 P + D = L5; (2) 같은 방식으로 다시 계산한 L5의 sd가 3d `annot_sd_common.csv`와 같다(차이 < 10⁻⁹).
- 기술: 대상 SNP 중 근위 비율, 각 L5 주석의 근위 비율.
- LD 점수: 3d·3f·3i와 같은 호출(1 cM, 1000 Genomes EUR GRCh38, HapMap3 회귀 SNP).

## 모형 (k02)
- **PD_{g}** = baselineLD v2.2 + `P_{g}_{CAGE, DNASE, H3K27ac}` + `D_{g}_{CAGE, DNASE, H3K27ac}` (6개 동시). 주 지도의 3층위 동시 모형에서 각 주석을 근위·원위로 나눈 것이다.
- 형질 14개 × 자기 계통 = 14회, 단일 실행(전체 트랙; 주 지도에서 전체 트랙과 재추출의 1위 층위가 11/11 같았다).
- jackknife: 3j의 정렬 블록(형질 사이 공분산 반영). 호출과 설정은 3j와 같다. 시험: 가짜 형질 1개(Z 섞음, 결과 아님).

## 가설과 판정 (α 사전 고정)
- 원위 대비: c1_D = CAGE_D − DNase_D, c2_D = CAGE_D − H3K27ac_D (τ\*). 근위도 같은 식(c1_P, c2_P).

| | 내용 | 판정 |
|---|---|---|
| **K-H1 (주)** | 원위에서 (a) 1위 층위: 주 지도 CAGE형 7개(AD, 단핵구 수, 교육 연수, 조현병, 양극성장애, BMI, LDL) 중 CAGE_D 1위 수, 키 DNase_D 1위 여부; (b) 지도 형질 11개의 (c1_D, c2_D) 다변량 Q (공동 jackknife, χ²) | (a) 7개 중 5개 이상 그리고 키 DNase **와** (b) P < 0.05 → **"원위에서 유지"**; 둘 중 하나 → "부분"; 둘 다 아님 → "원위에서 유지 안 됨" |
| K-H2 | 원위의 사전 대비 3개(AD 대 키 c1_D, 뇌 형질 4개 대 키 c1_D, IBD·RA 대 AD·단핵구 수 c2_D; 양측 α 0.05/3)와 검증 대비(호중구 수·총콜레스테롤 대 골밀도 c1_D; 단측) | 기술 |
| K-H3 | 근위의 같은 검정(1위 층위, Q, 대비) | 기술 |
| K-H4 | 형질마다 원위 1위가 2위보다 큰지(단측 P), F 근사 Q, 주석 크기 | 기술 |

- 기준 근거: 원위 주석은 원래 주석보다 작아 SE가 커진다. 그래서 1위 층위는 7개 중 2개까지 바뀌는 것을 허용하고, 개별 대비의 유의성은 판정에 넣지 않는다.
- τ\*는 주석 크기에 따라 달라지므로 근위와 원위의 τ\*를 서로 비교하지 않는다. 비교는 같은 구간(근위끼리, 원위끼리)의 층위 사이에서만 한다. 같은 구간 안에서도 층위마다 크기가 조금 다를 수 있으므로 크기를 함께 보고한다.

## 읽는 규칙
- **"원위에서 유지"**
  - 쓸 문장: "When each layer annotation was split at 1 kb from the nearest annotated transcription start site, the layer profiles of distal SNPs still differed between traits (*Q* = …, *P* = …), and CAGE remained the most informative distal layer for [n] of the 7 CAGE-type traits, with DNase for height. The CAGE type therefore does not reflect promoter proximity alone."
  - 한계 넷째를 이 결과로 바꾼다.
- **"부분" / "원위에서 유지 안 됨"**
  - CAGE 우위가 어느 구간에 있는지 보고한다. 초록·논의의 해석을 "predicted effects on transcription initiation near promoters"로 좁히고, 한계 넷째에 결과를 쓴다.
- K-H2–K-H4는 판정을 돕는 기술이다. 원위 주석이 작아 대비가 유의하지 않아도 그것만으로 결론을 바꾸지 않는다. 그 점을 함께 쓴다.

## 한계 (미리 적음)
- 1 kb 기준은 하나만 쓴다. 근위에 CAGE 점수 창(501 bp)이 TSS와 겹치는 SNP가 들어가도록 고른 값이다.
- 주석 TSS에는 없는 전사 시작 위치(인핸서 RNA 등)는 원위에 들어간다. 원위 CAGE 신호가 인핸서 RNA인지, 주석되지 않은 프로모터인지는 이 실험으로 가릴 수 없다.
- 근위·원위 주석은 원래 주석보다 작고, 그 크기가 층위마다 다르다.

## 산출물
- 스크립트: `scripts/k01_split_annot.py`, `k02_sldsc.py`, `k03_tests.py`, `run_pipeline.sh`
- 자료: `D:\AD_GWAS_Atlas_data\exp3k\`, TSS 원자료 `D:\AD_GWAS_Atlas_data\ref\gencode.v46.annotation.gtf.gz` (2026-10-15 EBI FTP에서 받음)
- 결과: `results/RESULTS_EXP3k_*.md`

## 실행 기록 (결과 보기 전)
- 2026-10-15 GENCODE v46 종합 주석 GTF(51,006,804 바이트)를 받았다.
- 2026-10-15 k01 시험(chr22) 통과 후 본 실행: 확인 통과(P + D = L5 모든 SNP, L5 sd 재현 최대 차이 8.3 × 10⁻¹⁷). TSS 고유 위치 190,108, 대상 SNP 근위 313,875 (5.4 %). L5 주석의 근위 비율: CAGE 16–19 %, DNase 5–9 %, H3K27ac 9–11 %. 원위 주석 크기: CAGE 4.0–4.2 %, DNase 4.6–4.8 %, H3K27ac 4.4–4.6 % (근위 0.2–1.0 %).
- 2026-10-15 dry(가짜 형질 1개, 결과 아님): k02·k03 형식 확인.
- 2026-10-15 10:39 사용자 확인으로 고정 → 본 실행(run_pipeline.sh: k02 main 14회 → k03 main).

## 결과 후 기록
- 2026-10-15 완료. **K-H1 판정: 원위에서 유지** — 원위 CAGE형 7/7 CAGE 1위, 키 DNase; 원위 Q 63.7, P 1.9 × 10⁻⁶ (F 1.0 × 10⁻⁴). 원위 AD 대 키 P 0.003, 뇌 대 키 7 × 10⁻⁵, 검증 대비 2 × 10⁻⁷. 근위(0.2–1.0 %)는 Q P 0.54로 형질이 갈리지 않음. 결과 문서: `results/RESULTS_EXP3k_2026-10-15.md`
