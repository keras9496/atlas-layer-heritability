"""실험 3w w10: 가중 모형(w09)으로 낸 다섯 주석의 평균과 stromal − 나머지 네 주석의 차이 — 다섯 가지 주석 구성(전체 트랙, 세 단위, 세포형 일치, 원위, 근위)에서 같은 방식으로 (원고의 대표값; 사후 분석, 새 적합 없음).
 모형과 추정은 w09와 같다(주석 평균 + GWAS 효과 + 쌍 고유 효과 + 추정 오차). 표준오차는 블록 jackknife(delete 값으로 다시 추정).
 12주석(전체 트랙) 모형의 평균과 중앙값, 쌍 수도 한 표로 낸다(그림 4와 보충 표용).
실행: python -I -X utf8 w10_constructions.py
출력: results/w10_means.csv (주석 × 구성: 모형 평균, s.e., 중앙값, 쌍 수), w10_contrasts.csv, w10_summary.txt
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy import stats
from w09_weighted import Model, jse, jcov, PG, LIN, FIVE, NON, NAME, KO, RES, D, valid
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
def fitjk(M):
    f = M.fit(); dl = [M.fit(y=M.DV[b], start=[max(f["t2"], 1.0), max(f["s2"], 1.0)]) for b in range(M.DV.shape[0])]; return f, np.array([x["beta"] for x in dl])
rows = []; crow = []
def one(lab, groups, setname):
    M = Model(lab, groups); f, B = fitjk(M); cov = np.linalg.inv(f["A"]); gi = {g: k for k, (kind, g) in enumerate(M.cols)}
    for g in M.groups:
        ys = np.array([D[(i, gg, lab)]["e"] for i, gg in M.pairs if gg == g]); rows.append(dict(set=setname, construction=lab, group=g, n_pairs=len(ys), model_mean=f["beta"][gi[g]], se_jackknife=jse(B[:, gi[g]]), se_model=float(np.sqrt(cov[gi[g], gi[g]])), median=float(np.median(ys))))
    tot = float(np.sqrt(f["t2"] + f["s2"])); P(f"  [{setname}; {NAME[lab]}] 쌍 {len(M.pairs)}개, GWAS {len(M.ids)}개 | " + "; ".join(f"{KO[g]} {f['beta'][gi[g]]:.1f}° ± {jse(B[:, gi[g]]):.1f}°" for g in M.groups) + f" | 참 표준편차 {tot:.1f}°")
    if "stromal" in M.groups and len(M.groups) == 5:
        c = np.zeros(len(M.cols)); c[gi["stromal"]] = 1
        for g in NON: c[gi[g]] = -0.25
        d = float(c @ f["beta"]); se = jse(B @ c); sm = float(np.sqrt(c @ cov @ c)); P(f"     stromal − 네 주석 평균: {d:+.1f}° ± {se:.1f}° (모형 s.e. {sm:.1f}°; P = {2 * stats.norm.sf(abs(d / se)):.1e})")
        crow.append(dict(construction=lab, contrast="stromal − mean of four", estimate=d, se_jackknife=se, se_model=sm, p=2 * stats.norm.sf(abs(d / se))))
        for g in NON:
            c2 = np.zeros(len(M.cols)); c2[gi["stromal"]] = 1; c2[gi[g]] = -1; d2 = float(c2 @ f["beta"]); crow.append(dict(construction=lab, contrast=f"stromal − {g}", estimate=d2, se_jackknife=jse(B @ c2), se_model=float(np.sqrt(c2 @ cov @ c2)), p=2 * stats.norm.sf(abs(d2 / jse(B @ c2)))))
        C = np.zeros((3, len(M.cols)));
        for r, g in enumerate(NON[1:]): C[r, gi[NON[0]]] = -1; C[r, gi[g]] = 1
        dd = C @ f["beta"]; Wd = float(dd @ np.linalg.solve(C @ cov @ C.T, dd)); Wj = float(dd @ np.linalg.solve(jcov(B @ C.T), dd))
        P(f"     네 주석의 평균이 같다: Wald = {Wd:.1f} (모형 공분산; P = {stats.chi2.sf(Wd, 3):.1e}), jackknife 공분산으로는 {Wj:.1f} (P = {stats.chi2.sf(Wj, 3):.1e})")
        crow.append(dict(construction=lab, contrast="equality of four non-stromal (Wald, 3 d.f.)", estimate=Wd, p=float(stats.chi2.sf(Wd, 3)), wald_jackknife=Wj, p_jackknife=float(stats.chi2.sf(Wj, 3))))
    return M, f, B
P("== 다섯 주석, 다섯 구성: 가중 모형의 평균 (± jackknife s.e.)")
K = {}
for lab in ("all", "res", "mat", "dist", "prox"): K[lab] = one(lab, FIVE, "five")
# 원위 − 근위: 두 모형의 delete 추정값을 같은 블록끼리 견준다
def con(lab):
    M, f, B = K[lab]; gi = {g: k for k, (kind, g) in enumerate(M.cols)}; c = np.zeros(len(M.cols)); c[gi["stromal"]] = 1
    for g in NON: c[gi[g]] = -0.25
    return float(c @ f["beta"]), B @ c
d1, v1 = con("dist"); d2, v2 = con("prox"); se = jse(v1 - v2); P(f"  원위의 (stromal − 네 주석) − 근위의 (stromal − 네 주석): {d1 - d2:+.1f}° ± {se:.1f}° (P = {2 * stats.norm.sf(abs((d1 - d2) / se)):.1e})")
crow.append(dict(construction="dist − prox", contrast="stromal − mean of four", estimate=d1 - d2, se_jackknife=se, p=2 * stats.norm.sf(abs((d1 - d2) / se))))
for g in FIVE:
    a = K["dist"]; b = K["prox"]; ga = [gg for _, gg in a[0].cols].index(g); gb = [gg for _, gg in b[0].cols].index(g); d = a[1]["beta"][ga] - b[1]["beta"][gb]; s = jse(a[2][:, ga] - b[2][:, gb]); P(f"  {KO[g]}: 원위 − 근위 = {d:+.1f}° ± {s:.1f}°")
    crow.append(dict(construction="dist − prox", contrast=g, estimate=d, se_jackknife=s, p=2 * stats.norm.sf(abs(d / s))))
P("  구성 사이의 차이: (stromal − 네 주석)의 차이, 같은 블록의 delete 추정값으로")
for x, y in (("mat", "all"), ("res", "all"), ("mat", "res")):
    dx, vx = con(x); dy, vy = con(y); se = jse(vx - vy); P(f"  {NAME[x]} − {NAME[y]}: {dx - dy:+.1f}° ± {se:.1f}° (P = {2 * stats.norm.sf(abs((dx - dy) / se)):.1e})")
    crow.append(dict(construction=f"{x} − {y}", contrast="stromal − mean of four", estimate=dx - dy, se_jackknife=se, p=2 * stats.norm.sf(abs((dx - dy) / se))))
P("\n== 12주석, 전체 트랙"); one("all", LIN, "twelve")
pd.DataFrame(rows).to_csv(RES/"w10_means.csv", index=False); pd.DataFrame(crow).to_csv(RES/"w10_contrasts.csv", index=False); (RES/"w10_summary.txt").write_text("\n".join(log), encoding="utf-8")
print(f"가장 큰 투영 기울기 {PG[0]:.1e}", file=sys.stderr)
