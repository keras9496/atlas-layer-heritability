"""실험 3w w09: 추정 오차를 넣은 가중 모형 — 조직 묶음 주석의 평균을 뺀 뒤 GWAS 사이의 참 차이 (PROTOCOL_exp3w_추가 B; 사후 분석, 새 적합 없음).
 모형: y_ig = μ_g + u_i + v_ig + e_ig;  Var(u_i) = τ²(GWAS), Var(v_ig) = σ²(쌍 고유), Cov(e_i) = S_i(같은 GWAS의 쌍들 사이의 jackknife 공분산; 알려진 값으로 둔다). GWAS 사이의 오차는 독립으로 둔다.
 추정: REML(τ² ≥ 0, σ² ≥ 0). τ² = 0, σ² = 0의 우도비 검정(경계: χ²0과 χ²1의 반반 혼합). 조직 묶음 평균이 같다는 Wald 검정. τ의 95% 구간은 프로파일 우도. 표준오차는 블록 jackknife(delete 값으로 y를 바꿔 다시 추정).
 수렴: 우도의 기울기를 닫힌 꼴로 넣어 투영 기울기가 0이 될 때까지 맞춘다(수치 기울기로는 우도가 평평한 곳에서 τ², σ²가 넷째 유효숫자부터 흔들려, S-LDSC를 다시 돌렸을 때 생기는
 마지막 자리 차이가 결과에 드러난다). 추정값이 경계(분산 0)에 있으면 우도비는 0이고 P = 1이다; LR0보다 작은 우도비는 두 로그우도의 반올림 오차이므로 0으로 둔다.
실행: python -I -X utf8 w09_weighted.py
출력: results/w09_summary.txt, w09_model.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy import optimize, stats
from wlib_panel import *
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LR0 = 1e-6; PG = [0.0]                                                          # PG: 맞춘 모형들의 투영 기울기 가운데 가장 큰 값(수렴 확인용, 표준 오류로만 적는다)
def lrt(ll, ll0): lr = 2 * (ll - ll0); return float(lr) if lr > LR0 else 0.0
def design(pairs, groups, with_cat):
    cols = [("g", g) for g in groups]; cats = sorted({cat(i) for i, _ in pairs})
    if with_cat: cols += [("c", c) for c in cats if c != "other"]
    Xm = np.zeros((len(pairs), len(cols)))
    for r, (i, g) in enumerate(pairs):
        Xm[r, cols.index(("g", g))] = 1
        if with_cat and cat(i) != "other": Xm[r, cols.index(("c", cat(i)))] = 1
    return Xm, cols
class Model:
    def __init__(self, lab, groups, with_cat=False, minpairs=1):
        v = valid(lab, groups); n = pd.Series([i for i, _ in v]).value_counts(); v = [(i, g) for i, g in v if n[i] >= minpairs]
        self.pairs = sorted(v); self.groups = [g for g in groups if any(gg == g for _, gg in self.pairs)]; self.ids = sorted({i for i, _ in self.pairs})
        self.X, self.cols = design(self.pairs, self.groups, with_cat); self.y = np.array([D[(i, g, lab)]["e"] for i, g in self.pairs]); self.DV = np.column_stack([D[(i, g, lab)]["dv"] for i, g in self.pairs])
        self.idx = [np.array([r for r, (j, _) in enumerate(self.pairs) if j == i]) for i in self.ids]; self.S = [np.atleast_2d(jcov(self.DV[:, ix])) for ix in self.idx]
    def parts(self, t2, s2, y, grad=False):
        A = 0; b = 0; ld = 0; Vi = []
        for ix, S in zip(self.idx, self.S):
            V = S + t2 + s2 * np.eye(len(ix)); L = np.linalg.cholesky(V); ld += 2 * np.log(np.diag(L)).sum(); W_ = np.linalg.inv(V); Vi.append(W_); Xi = self.X[ix]; A = A + Xi.T @ W_ @ Xi; b = b + Xi.T @ W_ @ y[ix]
        beta = np.linalg.solve(A, b); R = [y[ix] - self.X[ix] @ beta for ix in self.idx]; q = sum(float(r @ W_ @ r) for r, W_ in zip(R, Vi)); ll = -0.5 * (ld + np.linalg.slogdet(A)[1] + q)
        if not grad: return ll, beta, A
        Ai = np.linalg.inv(A); g = np.zeros(2)                                  # REML 로그우도의 기울기: ∂V/∂τ² = 11', ∂V/∂σ² = I
        for ix, W_, r in zip(self.idx, Vi, R):
            XW = self.X[ix].T @ W_; w1 = W_.sum(1); x1 = XW.sum(1); Wr = W_ @ r
            g[0] += w1.sum() - x1 @ Ai @ x1 - Wr.sum() ** 2; g[1] += np.trace(W_) - np.trace(Ai @ (XW @ XW.T)) - Wr @ Wr
        return ll, beta, A, -0.5 * g
    def fit(self, y=None, fix_t2=None, fix_s2=None, start=None):
        y = self.y if y is None else y; free = [k for k, v in enumerate((fix_t2, fix_s2)) if v is None]; nfree = len(free)
        def full(p):
            v = [fix_t2, fix_s2]
            for k, x in zip(free, p): v[k] = float(x)
            return v
        def nll(p): ll, _, _, g = self.parts(*full(p), y, grad=True); return -ll, -g[free]
        if nfree == 0: ll, beta, A = self.parts(fix_t2, fix_s2, y); return dict(ll=ll, t2=fix_t2, s2=fix_s2, beta=beta, A=A)
        starts = [start] if start is not None else [[1.0] * nfree, [100.0] * nfree, [400.0] * nfree] + ([[400.0, 1.0], [1.0, 400.0]] if nfree == 2 else [])
        best = None
        for s0 in starts:
            r = optimize.minimize(nll, s0, jac=True, method="L-BFGS-B", bounds=[(0, 1e5)] * nfree, options=dict(ftol=1e-15, gtol=1e-11, maxiter=5000, maxls=60))
            if best is None or r.fun < best.fun: best = r
        x = np.where(np.asarray(best.x) < 1e-9, 0.0, best.x); t2, s2 = full(x); ll, beta, A, g = self.parts(t2, s2, y, grad=True); gn = -g[free]          # 경계에 닿은 분산(1e-9 아래)은 0으로 둔다
        pg = float(np.max(np.abs(np.where(x <= 0, np.minimum(gn, 0), gn)))); PG[0] = max(PG[0], pg); assert pg < 1e-6, ("수렴하지 않음", t2, s2, pg)   # 경계에서는 안쪽으로 향하는 성분만 센다
        return dict(ll=ll, t2=t2, s2=s2, beta=beta, A=A)
def half(lr): return 0.5 * float(stats.chi2.sf(lr, 1)) if lr > 0 else 1.0
rows = []
def analyse(lab, groups, note, with_cat=False, minpairs=1, jack=True):
    M = Model(lab, groups, with_cat, minpairs); f = M.fit(); f0 = M.fit(fix_t2=0.0); fs = M.fit(fix_s2=0.0); lr_t = lrt(f["ll"], f0["ll"]); lr_s = lrt(f["ll"], fs["ll"])
    # τ의 프로파일 우도 구간
    grid = np.arange(0, 60.01, 0.5); prof = np.array([M.fit(fix_t2=float(t) ** 2)["ll"] for t in grid]); inside = grid[2 * (f["ll"] - prof) <= 3.841]; lo, hi = (float(inside.min()), float(inside.max())) if len(inside) else (np.nan, np.nan)
    # 조직 묶음 평균이 같다는 Wald 검정
    G = len(M.groups); cov = np.linalg.inv(f["A"]); C = np.zeros((G - 1, len(M.cols))); C[:, 0] = -1; C[np.arange(G - 1), np.arange(1, G)] = 1; d = C @ f["beta"]; Wd = float(d @ np.linalg.solve(C @ cov @ C.T, d))
    se_t2 = se_s2 = np.nan; se_b = np.full(len(M.cols), np.nan)
    if jack:
        dl = [M.fit(y=M.DV[b], start=[max(f["t2"], 1.0), max(f["s2"], 1.0)]) for b in range(M.DV.shape[0])]; se_t2 = jse([x["t2"] for x in dl]); se_s2 = jse([x["s2"] for x in dl]); se_b = np.array([jse([x["beta"][k] for x in dl]) for k in range(len(M.cols))])
    tot = f["t2"] + f["s2"]
    P(f"  [{note}] 쌍 {len(M.pairs)}개, GWAS {len(M.ids)}개, 묶음 {G}개")
    P(f"     GWAS 사이 τ = {np.sqrt(f['t2']):.1f}° (95% 프로파일 구간 {lo:.1f}°–{hi:.1f}°; τ² = {f['t2']:.0f} ± {se_t2:.0f} jackknife); τ² = 0의 우도비 = {lr_t:.1f}, P = {half(lr_t):.1e}")
    P(f"     쌍 고유 σ = {np.sqrt(f['s2']):.1f}° (σ² = {f['s2']:.0f} ± {se_s2:.0f}); σ² = 0의 우도비 = {lr_s:.1f}, P = {half(lr_s):.1e}; 참 분산 가운데 GWAS의 몫 τ²/(τ² + σ²) = {(f['t2'] / tot if tot > 0 else np.nan):.2f}")
    P(f"     조직 묶음 평균이 같다: Wald = {Wd:.1f} (자유도 {G - 1}), P = {stats.chi2.sf(Wd, G - 1):.1e}")
    P("     평균 μ: " + "; ".join(f"{KO[c[1]] if c[0] == 'g' else c[1]} {b:.1f}° ± {s:.1f}°" for c, b, s in zip(M.cols, f["beta"], se_b)))
    for c, b, s in zip(M.cols, f["beta"], se_b): rows.append(dict(set=note, term=c[1], kind={"g": "tissue group", "c": "trait category"}[c[0]], estimate=b, se_model=float(np.sqrt(cov[M.cols.index(c), M.cols.index(c)])), se_jackknife=s))
    rows.append(dict(set=note, term="tau", kind="between-GWAS s.d.", estimate=float(np.sqrt(f["t2"])), ci_lo=lo, ci_hi=hi, lrt_p=half(lr_t), n_pairs=len(M.pairs), n_gwas=len(M.ids)))
    rows.append(dict(set=note, term="sigma", kind="pair-specific s.d.", estimate=float(np.sqrt(f["s2"])), lrt_p=half(lr_s), n_pairs=len(M.pairs), n_gwas=len(M.ids)))
    return M, f
if __name__ == "__main__":
    P("== 가중 모형: y = 조직 묶음 평균 + GWAS 효과(τ) + 쌍 고유 효과(σ) + 추정 오차(jackknife 공분산)")
    M, f = analyse("all", LIN, "전체 트랙, 12묶음")
    A = gmeans("all", LIN); P("     (견줌) w05의 역분산 가중 평균: " + "; ".join(f"{KO[g]} {A[g][0]:.1f}°" for g in M.groups))
    analyse("all", LIN, "전체 트랙, 12묶음 + 형질 범주", with_cat=True)
    analyse("all", [g for g in LIN if g != "stromal"], "전체 트랙, stromal을 뺀 11묶음")
    for lab in ("all", "res", "mat"):
        analyse(lab, FIVE, f"{NAME[lab]}, 다섯 묶음"); analyse(lab, FIVE, f"{NAME[lab]}, 다섯 묶음 + 형질 범주", with_cat=True)
    for lab in ("all", "res", "mat"): analyse(lab, NON, f"{NAME[lab]}, stromal을 뺀 네 묶음", jack=False)
    pd.DataFrame(rows).to_csv(RES/"w09_model.csv", index=False); (RES/"w09_summary.txt").write_text("\n".join(log), encoding="utf-8")
    print(f"가장 큰 투영 기울기 {PG[0]:.1e}", file=sys.stderr)
