"""실험 3w w09b: w09의 모형에서 참 이질성 전체(τ² + σ²)의 크기와 검정 (사후 분석; w09의 τ와 σ는 GWAS 대부분이 다섯 묶음 가운데 하나에만 농축되어 서로 잘 갈리지 않으므로 합을 따로 적는다).
 검정: τ² = σ² = 0의 우도비, 경계 보정(χ²0 ¼, χ²1 ½, χ²2 ¼의 혼합). 표준오차는 블록 jackknife. 추정값이 경계(분산 0)에 있으면 우도비는 0이고 P = 1이다(w09의 lrt).
실행: python -I -X utf8 w09b_total.py
출력: results/w09b_summary.txt
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy import stats
from w09_weighted import Model, jse, lrt, half, PG, LIN, FIVE, NON, NAME, RES
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
def mix(lr): return float(0.5 * stats.chi2.sf(lr, 1) + 0.25 * stats.chi2.sf(lr, 2)) if lr > 0 else 1.0
P("== 조직 묶음 주석의 평균 둘레의 참 이질성(추정 오차를 뺀 것): sqrt(τ² + σ²)")
for lab, groups, note in (("all", LIN, "전체 트랙, 12묶음"), ("all", FIVE, "전체 트랙, 다섯 묶음"), ("res", FIVE, "세 단위, 다섯 묶음"), ("mat", FIVE, "세포형 일치, 다섯 묶음"),
                          ("all", NON, "전체 트랙, stromal을 뺀 네 묶음"), ("res", NON, "세 단위, stromal을 뺀 네 묶음"), ("mat", NON, "세포형 일치, stromal을 뺀 네 묶음")):
    M = Model(lab, groups); f = M.fit(); f0 = M.fit(fix_t2=0.0, fix_s2=0.0); lr = lrt(f["ll"], f0["ll"]); tot = f["t2"] + f["s2"]
    dl = [M.fit(y=M.DV[b], start=[max(f["t2"], 1.0), max(f["s2"], 1.0)]) for b in range(M.DV.shape[0])]; se = jse([x["t2"] + x["s2"] for x in dl])
    rms = float(np.sqrt(np.mean([np.diag(S).mean() for S in M.S]))); sd = float(np.std(M.y - M.X @ np.linalg.lstsq(M.X, M.y, rcond=None)[0], ddof=len(M.groups)))
    P(f"  [{note}] 쌍 {len(M.pairs)}개, GWAS {len(M.ids)}개: 참 표준편차 {np.sqrt(tot):.1f}° (τ² + σ² = {tot:.0f} ± {se:.0f}); 이질성 없음의 우도비 = {lr:.1f}, P = {mix(lr):.1e}; 견줌: 묶음 평균을 뺀 잔차의 표준편차 {sd:.1f}°")
P("== 주석 하나 안의 참 이질성(GWAS마다 쌍이 하나이므로 σ만 추정; 이질성 없음의 우도비는 χ²0과 χ²1의 반반 혼합)")
for g in FIVE:
    out = []
    for lab in ("all", "res", "mat"):
        M = Model(lab, [g]); f = M.fit(fix_t2=0.0); f0 = M.fit(fix_t2=0.0, fix_s2=0.0); pv = half(lrt(f["ll"], f0["ll"]))
        out.append(f"{NAME[lab]} {np.sqrt(f['s2']):.1f}° (P = {pv:.1e}; 평균 {f['beta'][0]:.1f}°; 쌍 {len(M.pairs)})")
    P(f"  {g}: " + " | ".join(out))
(RES/"w09b_summary.txt").write_text(chr(10).join(log), encoding="utf-8")
print(f"가장 큰 투영 기울기 {PG[0]:.1e}", file=sys.stderr)
