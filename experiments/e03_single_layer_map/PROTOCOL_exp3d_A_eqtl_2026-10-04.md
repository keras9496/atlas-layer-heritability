# 실험 3d-A: Atlas CAGE 대 eQTL·기존 전사 주석 비교 (2026-10-04, 사후 기술 분석, 비교 방법 먼저 고정)

## 질문
"Atlas 예측 CAGE(국소 전사 시작) 효과가 큰 SNP에 유전율이 몰린다"는 결과는 기존 GWAS 해석 도구보다 무엇을 더 알려 주는가? 첫 단계로, 기존의 대표 도구인 **fine-mapped eQTL**과 **실측 CAGE 인핸서**에 비교한다.

## 사후성 명시
- 3c·3d 결과를 본 뒤 정한 **사후 기술 분석**이다. 판정(α)은 두지 않는다.
- 프로토콜을 쓰기 전, 결과 파일 형식을 확인하다가 `AD__D1__L5_myeloid_CAGE.results` 한 파일의 baseline 계수 줄을 봤다(GTEx_eQTL_MaxCPP 계수 z 약 2.9, Enhancer_Andersson z 약 0.6). 다른 형질의 값은 보지 않았다.

## 비교 대상 (baselineLD v2.2에 이미 들어 있는 주석)
| 주석 | 뜻 |
|---|---|
| GTEx_eQTL_MaxCPP | GTEx 44개 조직 eQTL fine-mapping의 최대 인과 확률 (Hormozdiari 2018) |
| BLUEPRINT_H3K27acQTL / H3K4me1QTL / DNA_methylation_MaxCPP | 혈액 면역세포 분자 QTL의 최대 인과 확률 |
| Enhancer_Andersson | FANTOM5 CAGE로 정의한 인핸서(조직 합산) |
| TSS_Hoffman, Promoter_UCSC | 일반 TSS·프로모터 |

## Atlas 쪽 (형질마다 미리 고정, 3c에서 유의했던 계통)
| 형질 | 계통 |
|---|---|
| AD, 단핵구 수 | 골수계 |
| 교육 연수, 조현병, 양극성장애, 신경증, BMI | 뇌 조직 |
| IBD, RA | 림프계 |
| LDL | 간 |
| 키 | 결합조직 |

- 각 계통의 **CAGE**와 **RNA**를 쓴다.
  - D1 단일 모형(baselineLD + L5_{g,m})과, 트랙 맞춤 L5k를 같은 방식으로 넣은 단일 모형(새 실행)을 쓴다.

## 실행
- **새 실행**
  - 기준 모형(baselineLD만): 형질 11개, 11회
  - 트랙 맞춤 단일 모형(L5k_{g,CAGE}, L5k_{g,RNA}): 형질 11 × 2 = 22회
  - 순차, 총 33회
- 기존 실행: D1의 L5_{g,CAGE}, L5_{g,RNA} (22개 결과 파일)
- τ\* 계산: τ · sd(주석) · M_5_50 / h²_obs
  - sd는 MAF ≥ 5 % 참조 SNP(1000G EUR hg38 frq)에서 baselineLD annot 파일로 계산한다.
  - Atlas 주석의 sd는 3d `annot_sd_common.csv`를 쓴다.

## 보고 (기술)
1. 형질마다 같은 모형 안의 τ\*를 나란히 보인다.
   - Atlas CAGE
   - Atlas RNA
   - GTEx_eQTL_MaxCPP
   - BLUEPRINT QTL 3개
   - Enhancer_Andersson
2. **eQTL 감쇠**: GTEx_eQTL_MaxCPP τ\*가 기준 모형에서 Atlas CAGE를 넣은 모형으로 얼마나 바뀌는가. 크게 줄면 두 주석이 같은 유전율을 나눠 갖는다는 뜻이다.
3. 읽는 규칙 (미리 고정)
   - Atlas CAGE τ\*가 eQTL τ\*보다 크고 트랙 맞춤에서도 유지되면: "유전율은 fine-mapped eQTL보다 예측 국소 전사 시작 효과에 더 많이 정렬된다"
   - 비슷하거나 작으면: "Atlas CAGE는 eQTL과 비슷한 수준의 정보를 준다"
   - 형질 11개 중 몇 개가 첫째 경우인지 센다. 판정이 아닌 기술이다.
   - RNA는 같은 비교에서 대조로 쓴다(예측 유전자 수준 발현).

## 한계
- τ\*끼리의 비교는 표준화 척도라 가능하지만, SE가 겹치면 "더 크다"고 쓰지 않는다. 차이의 SE는 계산하지 않는다(서로 다른 실행이라 공분산을 모른다). 같은 모형 안의 쌍만 part_delete로 차이 SE를 계산한다.
- GTEx eQTL MaxCPP는 여러 조직을 합친 값이다. 조직 특이 eQTL과의 비교는 아니다.

## 출력
- `scripts/y06_eqtl_compare.py`, `results/y06_eqtl_compare.csv`, `results/y06_eqtl_compare.txt`

## 결과 후 기록 (2026-10-04)
- 실행: 새 S-LDSC 33회(순차) 완료
  - 골수계 RNA는 트랙이 2개라 L5k와 L5가 같다. 같은 결과가 나오는 것이 정상이다.
- 읽는 규칙 적용
  - Atlas CAGE τ\* > eQTL τ\*: 점추정 11/11, 같은 모형 차이 p < 0.05는 D1 7/11, D1k 9/11
  - eQTL 감쇠는 6–25 %(RA −45 %, LDL −60 %)
- 결과 문서: `results/RESULTS_EXP3d_A_eqtl_2026-10-04.md`
