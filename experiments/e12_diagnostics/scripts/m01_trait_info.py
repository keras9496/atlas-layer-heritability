"""실험 3m m01: 형질별 정보량 지표를 기존 실행 기록에서 뽑는다 (교수님 피드백 2026-10-17, 보충 표 S1 확장; 새 실행 없음).
 원자료: 3j 정렬 블록 재실행의 주 모형 로그 exp3j/runs/<형질>__T3all.log (baselineLD v2.2 + 자기 계통 3층위, 전체 트랙).
 뽑는 값: 관측 척도 h²(SE), Z = h²/SE, 평균 χ², λGC, 절편(SE), 비율(SE), 회귀 SNP 수(χ² 상한 제거 뒤).
 h²·절편은 모형에 따라 조금 달라지지만(주석 추가) 차이는 작다. 표에는 이 모형의 값을 쓴다고 밝힌다.
출력: results/m01_trait_info.csv
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re
import pandas as pd
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; R = Path(f"{DATA}/exp3j/runs"); (EXP/"results").mkdir(exist_ok=True)
TRAITS = ["AD", "MONO", "EA", "SCZ", "BIP", "NEURO", "BMI", "IBD", "RA", "LDL", "HEIGHT", "BMD", "NEUT", "TC"]
num = r"([-\d.eE+]+)"
rows = []
for t in TRAITS:
    s = (R/f"{t}__T3all.log").read_text(encoding="utf-8", errors="replace")
    g = lambda p: re.search(p, s).groups()
    h2, h2se = map(float, g(rf"Total Observed scale h2: {num} \({num}\)"))
    icpt, icse = map(float, g(rf"Intercept: {num} \({num}\)"))
    m = re.search(rf"Ratio: {num} \({num}\)", s); ratio, rse = (map(float, m.groups()) if m else (float("nan"), float("nan")))
    if not m and "Ratio < 0" in s: ratio = rse = float("nan")
    nreg = int(re.findall(r"\((\d+) SNPs remain\)", s)[-1])
    rows.append(dict(trait=t, h2=h2, h2_se=h2se, h2_z=h2 / h2se, mean_chi2=float(g(rf"Mean Chi\^2: {num}")[0]), lambda_gc=float(g(rf"Lambda GC: {num}")[0]),
                     intercept=icpt, intercept_se=icse, ratio=ratio, ratio_se=rse, n_reg_snp=nreg))
D = pd.DataFrame(rows); D.to_csv(EXP/"results"/"m01_trait_info.csv", index=False)
print(D.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
