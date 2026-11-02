# 실험 3r: 더 큰 AD GWAS로 본 AD의 지도 위치와 활성도를 넘는 변이 효과 (2026-10-25 작성, 2026-10-25 결과 보기 전 고정)

> **고정(2026-10-25, 사용자 확인)**: "3r 3s 이렇게 분석 갑시다!". 이 시점까지 새 자료에서 본 것은 파일 머리글·행 수·md5(3r)와 파일 형식·md5(3s)뿐이다.

## 배경
- 2026-10-25 받은 AD 확장 질문 3번: "단핵구·대식세포 등에서 예측한 기준 allele의 활성도를 충분히 조건화한 뒤에도 CAGE와 DNase의 변이 효과 점수가 추가적인 유전력 정보를 제공하는지 … 최신 AD GWAS를 이용하고 … 활성도가 비슷한 영역 안에서도 예측 효과가 큰 변이에 위험 신호가 집중되는지 확인".
- 사용자 결정(2026-10-25): "B (권장): 글 + S-LDSC 확장 둘". 3r은 그중 질문 3(가), 3s는 질문 2(A)다. 검토서: `Atlas_AD_Genomewide 실험 (1-3)/paper/manuscript_v1/검토_AD확장질문3_2026-10-25.md`.
- 지금까지 Bellenguez 2022 stage I(APOE 제외)로 얻은 AD 결과(이미 본 것): 골수계 CAGE 효과 τ\* E +2.91 ± 0.92, EA +0.42 ± 0.94, EAc +2.33 ± 0.89, S −0.21 ± 0.65, SA +0.36, EAM +0.52; 활성도 상위 10 % 안 Δ(HL10) CAGE −0.19 ± 1.66, DNase +2.21(단측 P 0.028, 기술); HL20 CAGE −0.50 ± 1.17, DNase +1.35. h² 0.026, h² z 10.7. 걸림돌은 AD의 정보량이다.

## 사후성 명시
- 이미 본 것: 3c–3q 전체 결과(위 AD 수치 포함). 특히 HL10 DNase Δ의 단서(P 0.028)를 본 뒤에 이 프로토콜을 쓴다.
- 아직 보지 않은 것: 새 AD GWAS의 모든 결과. munge(평균 χ² 포함)와 S-LDSC는 이 프로토콜을 고정한 뒤에만 돌린다. 고정 전에는 머리글과 행 수만 확인한다.

## 자료 (r00, 공개)
| 항목 | 내용 |
|---|---|
| 연구 | 2026 Nature Genetics AD 합의 메타분석(DOI 10.1038/s41588-026-02583-1; ADGC, Bonn, CHARGE, EADB, EADI, FinnGen, GERAD, GR@ACE/DEGESCO, PGC-ALZ 등 유럽계 52개 연구) |
| 파일 | GWAS Catalog **GCST90704647** (proxy 제외판), `GCST90704647.tsv.gz`, GRCh38, GWAS-SSF, md5 55e24a975113bb767fe11605eec87af2 |
| 표현형 | "Alzheimer's disease or related dementias", 환자 75,638명 / 대조 1,043,805명, proxy 제외 |
| 열(머리글 확인) | chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency, p_value, variant_id, rs_id, N_cases, N_controls, Neff_cases, Neff_controls, Neff_total, HetISq, HetPVal |
| 선택 이유 | 세 공개 파일(주 분석·proxy 제외·biobank 제외) 중 저자들이 LDSC에 쓴 판. proxy 희석이 없다 |
| 원래 자료와 겹침 | Bellenguez 2022(EADB 등)를 포함한다 → **독립 재현이 아니라 "더 큰, 겹치는 GWAS로 본 강건성"**(3n과 같은 위치)으로만 보고한다 |

## 요약통계 처리 (r01; 원래 AD 처리 `code/00_prepare_sumstats.py`와 같은 규칙)
- SNV만(effect·other allele 한 글자), **APOE 영역 제외(chr19:43.9–45.9 Mb, GRCh38; 주 분석과 같음)**, rs_id로 HapMap3 no-MHC SNP만, rsID 중복은 첫 행.
- Z = beta / standard_error, P = p_value, FRQ = effect_allele_frequency.
- N = Neff_total(SNP별 유효 표본 수; 3n의 환자–대조 형질과 같은 규칙). 원래 AD는 N_cases + N_controls(원 표본 수)를 썼다. τ\*는 N의 상수배에 영향받지 않는다.
- `munge_sumstats.py --signed-sumstats Z,0`(3n n01과 같다). 확인: HapMap3 no-MHC 목록의 80 % 이상이 남아야 한다.
- 형질 이름: `N_AD`.

## 모형 (r02; 새 주석 없음, 모두 기존 주석; 계통 = 골수계)
설정은 3j와 같다: baselineLD v2.2, 정렬 블록 200개, 같은 회귀 가중치와 빈도 파일, 같은 실행기.

| 모형 | 주석 | 출처 | 실행 수 |
|---|---|---|---|
| T3_01–T3_20, T3all | 층위 지도(재추출 20회 + 전체 트랙) | 3f, 3d | 21 |
| M1_<계통> | 12개 계통 단일 주석 | 3c | 12 |
| EA, EAM | 효과 + 이진 활성도(+ 실측) | 3l | 2 |
| EAc, S, SA | 연속 활성도 / 활성도 층화 효과 / 층화 + 이진 | 3p | 3 |
| HL10, HL20 | 활성도 상위 10 %·20 % 안 활성도 맞춤 효과 상·하 반분 | 3q | 2 |

- 새 실행: 40회. 확인 단계(check): 원래 AD GWAS로 T3all·EA·EAc·HL10을 이 실행기로 다시 돌려 기존 npz와 계수·delete 값 상대 차이 < 10⁻⁹이어야 한다.

## 가설과 판정 (α 사전 고정)
| | 내용 | 판정 |
|---|---|---|
| **R-H1 (주)** | 활성도 상위 10 % 안에서 활성도를 맞춘 효과 상·하 반분의 차이 Δ = τ\*(H) − τ\*(Lo), 골수계 **CAGE**와 **DNase** 두 검정(단측, 각각 α 0.025) | 층위마다: P < 0.025 → **"검출"**; 아니고 95 % 구간 상한(Δ + 1.96 SE) < +2.28(3q 키 DNase Δ) → **"키 크기의 효과는 배제"**; 그 밖 → **"미해결(검정력 부족)"** |
| R-H2 | 자기 층위(CAGE) 효과 τ\*: EA, EAc, S, SA, EAM (단측 P, 유지 비율 = τ\*/E) | 기술. 3p 읽는 규칙을 그대로 쓴다: 공변량(EA·EAc)과 층화(S)가 갈리면 "구성에 의존" |
| R-H3 | HL20(상위 20 %)의 Δ, CAGE·DNase | 기술(민감도, R-H1과 같은 세 범주로 적음) |
| R-H4 | 지도 위치: 12개 계통 τ\*와 골수계 순위; 층위 τ\*와 c1(CAGE − DNase), c2(CAGE − H3K27ac)(재추출 평균); AD 대 키 c1(원래 키, 3n 새 키 각각; 공동 jackknife); 새 AD로 바꾼 지도의 Q(나머지 10형질은 주 지도, 그리고 3n 갱신 지도 두 판) | 기술 |
| R-H0 | 정보량: h², SE, h² z, 평균 χ², 절편, 비율, 회귀 SNP 수 | 기술. Bellenguez(h² z 10.7)와 나란히 |

- 기준 근거
  - R-H1의 두 층위: 질문이 CAGE와 DNase를 함께 묻는다. DNase는 Bellenguez에서 단서(P 0.028)를 본 뒤라 같은 참가자를 포함하는 자료에서 다시 보는 것이다. 그래서 결과를 "확인"이 아니라 "더 큰 GWAS에서도 보이는가"로만 적는다.
  - "키 크기 배제" 범주: 검출되지 않을 때도 정보가 남게 하려는 것이다. 키 DNase Δ(+2.28)는 3q에서 가장 강건했던 활성 요소 안 효과다. 상한이 그보다 작으면 "AD에서 활성 요소 안 효과 정보는 있더라도 키보다 작다"고 쓸 수 있다.
  - 원래 GWAS와의 차이 검정은 하지 않는다(같은 참가자 포함).

## 읽는 규칙
- R-H1 결과에 따라 원고의 AD 문장을 정한다(수정안 v2e).
  - CAGE "검출" → "더 큰 AD GWAS에서는 활성 요소 안 CAGE 효과 정보가 보였다(겹치는 자료)".
  - CAGE "키 크기 배제" → "활성 요소 안 CAGE 효과 정보는 검출되지 않았고, 있더라도 키의 접근성 효과보다 작다".
  - CAGE "미해결" → "검정력 부족으로 배제할 수 없다"(지금 원고 수정 방향 그대로).
  - DNase도 같은 세 범주로 따로 적는다.
- R-H0에서 h² z가 10.7보다 크지 않으면 "더 큰 GWAS가 AD 정보량을 늘리지 않았다"고 함께 적는다.
- 지도 본체(주 지도)는 원래 GWAS 그대로 둔다. 3r은 강건성·사용 예시로만 보고한다.

## 한계 (미리 적음)
- 새 GWAS는 Bellenguez를 포함하므로 독립 검증이 아니다. FinnGen도 들어 있다.
- 표현형이 "AD 또는 관련 치매"다. 원래 AD(임상 + proxy)와 정의가 다르다.
- 예측 트랙에 미세아교가 없다. 골수계(단핵구·대식세포)가 대리다.
- S-LDSC는 게놈 전역 평균이라 효과가 몰린 좌위를 가리키지 못한다.

## 산출물
- 스크립트: `scripts/r00_download.sh`(완료), `r01_munge.py`, `r02_sldsc.py`(3n n02·3l l03·3p p02·3q q02 실행기 통합), `r03_tests.py`(3p p_common 재사용)
- 자료: `D:\AD_GWAS_Atlas_data\exp3r\`
- 결과: `results/RESULTS_EXP3r_*.md`

## 실행 기록 (결과 보기 전)
- 2026-10-25 r00 내려받기 완료(공개 FTP, 598 MB). md5 55e24a975113bb767fe11605eec87af2가 공개 목록과 일치했다. 행 수는 머리글 포함 14,014,755다. 머리글과 행 수만 확인했고 값은 보지 않았다.
- 2026-10-25 고정 뒤 r01: APOE 영역 10,611행 제외, HapMap3 no-MHC SNV 1,203,752개(99 %, 기준 80 % 통과), munge 통과 1,178,160 SNP. Neff_total 중앙값 197,002(범위 61,196–197,002). P가 0인 SNP 없음.
- r02 check 통과: 원래 AD T3all·EA·EAc·HL10의 계수·delete 값 차이 0. r02 main 40회 완료(약 9분).

## 결과 후 기록
- 2026-10-25 r03 검정. 결과 문서 `results/RESULTS_EXP3r_2026-10-25.md`.
- R-H0: h² z 10.7 → 13.1. **R-H1 CAGE Δ +0.35 ± 1.44, DNase Δ +0.76 ± 1.01 → 둘 다 "미해결"**(Bellenguez DNase 단서 +2.21은 줄어듦). R-H3(상위 20 %): CAGE +0.00 ± 0.97, DNase −0.02 ± 0.87 → 둘 다 "키 크기의 효과는 배제".
- R-H2: E +2.91, EA +0.66, EAc +2.16(P 0.005), S −0.41, SA +0.05, EAM +0.78 → Bellenguez와 같은 "구성에 의존". EA 골수계 CAGE 활성도 +1.63 ± 0.60.
- R-H4: 골수계 1위 유지, c1 +3.76, AD 대 키 +5.19(P 4 × 10⁻⁴; 3n 새 키와는 +4.36, P 0.002), AD를 바꾼 지도 Q 73.6.
