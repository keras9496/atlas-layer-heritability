# 실험 3: AlphaGenome Atlas 변이 효과 주석 × Bellenguez 계층적 LD 점수 회귀 (S-LDSC) (2026-09-25, 결과 보기 전 고정)

## 질문
AD 유전 위험(SNP 유전율)은 **Atlas가 특정 세포형·특정 조절 모달리티(CAGE, DNase, 히스톤, RNA)에서 큰 효과를 예측한 흔한 변이**에 농축되어 있는가?
- 개인 자료는 쓰지 않는다. Bellenguez 요약통계(약 49만 명) 전체의 검정력으로 게놈 전역에서 본다.
- 선행 방법: Dey et al. 2020(DeepSEA·Basenji 변이 효과 주석 S-LDSC). 우리는 이를 AlphaGenome Atlas, AD, 세포형 × 모달리티로 확장한다.
- 실험 1·2와의 관계
  - 실험 1: 유전자좌 안에서만 봤다.
  - 실험 2: 개인 표현형으로 검정했는데, 검정력이 부족했다.
  - 실험 3: 같은 질문을 게놈 전역·요약통계 수준에서 검정한다.

## 자료
- GWAS: GCST90027158 stage I (hg38)
  - HapMap3 SNP으로 `munge_sumstats`, N = n_cases + n_controls
  - APOE(chr19:43.9–45.9 Mb)는 미리 제외한다. MHC는 가중치 파일에서 제외된다.
- 참조: Zenodo 10515792 S-LDSC 참조 묶음 (GRCh38)
  - 1000G Phase 3 EUR plink, baselineLD v2.2 (hg38), HapMap3 가중치(MHC 제외), frq
  - 도구: LDSC Python 3 포트
- **Atlas 주석 대상**: 참조 plink의 SNP 중 **MAF ≥ 5 %** 전부(약 600만 예상).
  - REF/ALT는 1000G 고심도 hg38 pvar로 맞춘다.
  - 스코어러: DNASE, CHIP_HISTONE (H3K27ac, H3K4me1, H3K4me3), CAGE, RNA_SEQ, AVI_SCORE
  - 맥락 5개: myeloid, neuron, astrocyte (1단계 정규식) / liver, nonliver (S3 정규식, 전체 풀)
  - 값: 트랙 max |ALT − REF|
  - MAF < 5 % SNP은 모든 Atlas 주석 = 0이다. 주석 정의 자체가 "흔한 SNP 중 Atlas 효과 상위"다. 유전율도 표준대로 MAF ≥ 5 % 기준(`--frqfile`)으로 추정한다.

## 주석 (결과 보기 전 정의)
- 백분위: 주석 대상 흔한 SNP 전체에서 (맥락, 모달리티)별로 매긴다.
- **세포형 주석 5개 (주)**: A_c = 1[S_c 백분위 ≥ 0.95]
  - S_c = DNase·H3K27ac·H3K4me1·CAGE 백분위의 평균(실험 1·2와 같은 정의)
  - 이진, 상위 5 %
- **모달리티 주석 30개**: A_{c,m} = 1[m 백분위(맥락 c) ≥ 0.95]
- 대조 주석: AVI 상위 5 %(세포형 비특이 Atlas 종합 점수)

## 모형과 통계량
- 각 주석은 **baselineLD v2.2(97개 주석)에 조건부**로 하나씩 넣는다.
  - 결과 지표: 농축도(유전율 비율 / SNP 비율)와 **표준화 효과 τ\***(Gazal 2017)
  - τ\*는 baseline을 넘는 고유 기여이며, 주 지표다.
- p값: block jackknife(200 블록) z 검정, τ > 0 단측.

## 주 가설 (α 사전 고정)
| | 내용 | 검정 | α |
|---|---|---|---|
| **E3-H1** | 세포형 주석 5개 각각의 τ\* > 0 (baselineLD 조건부) | 단측 | 0.05/5 = 0.01. **골수계가 주 관심.** |
| **E3-H2 (특이성)** | baselineLD + 세포형 주석 5개를 **동시에** 넣은 모형에서 τ\*_골수계 − τ\*_뉴런 > 0 | 제거 블록 값으로 차이의 jackknife SE | 0.05 |
| **E3-H3 (모달리티)** | 모달리티 주석 30개 각각 τ\* > 0 (baselineLD 조건부) | 단측, Bonferroni | 0.05/30 |

- 세포형마다 어떤 모달리티(CAGE vs 히스톤 vs DNase vs RNA)가 기여하는지 기술한다.

## 부 분석 (기술)
- AVI 상위 5 %를 추가로 조건부로 넣기: 세포형 정보가 Atlas 종합 점수를 넘는가(실험 1 S1과 같은 질문)
- 임계값 민감도: 상위 1 %, 상위 10 %
- 연속 주석: 백분위 값 그대로
- 골수계 vs 간 차이(실험 1 S3)
- APOE를 포함한 분석
- 실험 1 결과와 방향이 일치하는지 비교

## 해석 주의
- Atlas에는 미세아교세포 트랙이 없다(골수계 = CD14 단핵구 대리).
- 변이 효과 주석은 LD로 퍼진다. S-LDSC의 baselineLD 조건부가 LD·MAF·기존 기능 주석을 보정하지만, 완전하지 않을 수 있다.
- Bellenguez에는 대리 사례(가족력)가 포함되어 있어 관측 척도 h²가 작다. 농축도와 τ\*는 척도와 무관하게 해석한다.

## 실행 (순차)
| 순서 | 내용 | 예상 시간 |
|---|---|---|
| x01 | 참조 묶음 준비, 흔한 SNP의 REF/ALT 결정 | |
| x02 | Atlas 조회 (재시작 가능, 할당량 대응) | 약 2.6일 |
| x03 | 주석 파일 + LD 점수 (`ldsc --l2`) | |
| x04 | munge + S-LDSC | |
| x05 | 요약·그림 | |

- 코드·결과: `실험3_SLDSC/`. 대용량 자료: `D:\AD_GWAS_Atlas_data\exp3\`.

## 실행 기록 (결과 보기 전)
- x01: 참조의 MAF ≥ 5 % SNP 5,956,110개 중 REF/ALT를 정한 5,828,323개(97.9 %)를 Atlas로 조회한다. 나머지는 주석 0.
- 도구
  - LD 점수는 abdenlab/ldsc-python3에 패치를 적용해 계산한다(pandas `.ix` 제거, Windows gzip, M 파일 쓰기 모드).
    - 검증: chr22 "base" 주석의 L2가 공식 baselineLD의 baseL2와 완전히 일치했다(상관 1.0, 최대 차이 0). M_5_50도 일치했다.
  - 회귀는 brainvisa/ldscore 포크로 한다. abdenlab 포트의 회귀 부분은 염색체별 파일을 읽을 때 염색체 수 인자가 빠지는 버그가 여러 곳에 있어서 쓰지 않는다.
    - 검증: 무작위 5 % 가짜 주석의 농축도 1.80 ± 2.50 (p 0.75)로 귀무와 일치했다. h²_obs 0.026, 절편 1.049.
    - part_delete로 계산한 jackknife SE가 보고된 SE와 일치했다.
- munge: HapMap3(MHC 제외) 1,177,203 SNP, APOE 제외, 평균 χ² 1.269.
- 자동 순차 실행 (2026-09-25 22:16, 결과 보기 전): `scripts/run_pipeline.sh`를 Windows 분리 프로세스로 띄웠다(세션·compact와 무관).
  - x02 완료(583 parts)를 기다린다. x02가 도중에 죽으면 다시 시작한다. 끝났는데 멈춰 있으면 10분 뒤 종료한다. 그다음 x03 → x04 → x05를 순서대로 실행한다.
  - 상태: `D:\AD_GWAS_Atlas_data\exp3\PIPELINE_STATUS.txt`. 단계별 로그: `exp3\x0N_*.log`.
  - 사전 점검: x03을 시험 폴더(`exp3\test\annot`)에서 chr2로 돌렸다. 주석당 약 5 %, M·M_5_50·ct5set 형식 정상. x05(요약·그림)는 이때 작성했다.
