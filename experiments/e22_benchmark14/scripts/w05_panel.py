"""실험 3w w05: 106개 GWAS — 농축된 GWAS–조직 묶음 쌍을 단위로 한 CAGE–DNase 균형 (PROTOCOL_exp3w 3절).
 단위: 조직 묶음 모형(3u M1, 정렬 블록)의 τ* > 0, 단측 P < 0.05/12인 쌍. 형질마다 분석 묶음을 고르는 규칙과 예측 묶음은 쓰지 않는다. 형질 범주(3u traits의 group)는 기술에만 쓴다.
 쌍마다: 세 층위 함께 적합한 모형의 균형 φD = atan2(DNase τ*, CAGE τ*) — (가) 전체 트랙(3u T3all), (나) 층위당 세 단위 재추출 20회 평균(3u·3w T3r; 다섯 묶음), (다) 세포형 일치(3w M3; 다섯 묶음).
 조직 묶음 주석마다: 쌍들의 역분산 가중 평균(공동 jackknife)과 중앙값. 주석 사이의 차이. 균형의 분산을 조직 묶음과 GWAS가 설명하는 몫(η²).
 같은 GWAS가 여러 묶음에 농축된 경우: 묶음 사이의 차이(GWAS 안의 비교). 같은 주석 안: 골격 형질과 그 밖의 형질.
실행: python -I -X utf8 w05_panel.py
출력: results/w05_pairs.csv, w05_annotations.csv, w05_summary.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
EXP = Path(__file__).resolve().parents[1]; RES = EXP/"results"; DATA = Path(DATA); X = DATA/"exp3"; Y = DATA/"exp3d"; F3 = DATA/"exp3f"; I3 = DATA/"exp3i"; U = DATA/"exp3u"; W = DATA/"exp3w"
U3 = EXP.parent/"e20_gwas_panel"
LIN = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
FIVE = ["CNS_tissue", "liver", "lymphoid", "myeloid", "stromal"]; LAY = ["CAGE", "DNASE", "H3K27ac"]
M_TOT = sum(float(open(X/"GRCh38"/"baselineLD_v2.2"/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
SD = pd.concat([pd.read_csv(p/"annot_sd_common.csv", index_col=0).sd_common for p in (Y, F3, I3)]); assert not SD.index.duplicated().any()
strip = lambda n: re.sub(r"L2_\d+$", "", n)
def jse(v): v = np.asarray(v, float); B = len(v); return float(np.sqrt((B - 1) / B * np.sum((v - v.mean()) ** 2)))
def load(path):
    z = np.load(path); h2 = float(z["h2_log"]); d = {}
    for k, n in enumerate([strip(n) for n in z["names"]]):
        if n in SD.index: f = SD[n] * M_TOT / h2; d[[m for m in LAY if n.endswith("_" + m)][0]] = (float(z["coef"][k] * f), z["part_delete"][:, k] * f)
    assert sorted(d) == sorted(LAY), path; return d
def avg(ds): return {m: (float(np.mean([d[m][0] for d in ds])), np.mean([d[m][1] for d in ds], axis=0)) for m in LAY}
def phi(Lr):
    c, d = Lr["CAGE"], Lr["DNASE"]; e = float(np.degrees(np.arctan2(d[0], c[0]))); dv = np.degrees(np.arctan2(d[1], c[1])); dv = e + (dv - e + 180) % 360 - 180
    return e, dv, bool(c[0] + d[0] > 0), c[0] - d[0], c[1] - d[1]
def ivw(items):
    w = np.array([1 / jse(d) ** 2 for _, d in items]); w = w / w.sum(); return float(sum(wi * e for wi, (e, _) in zip(w, items))), sum(wi * d for wi, (_, d) in zip(w, items))
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
TR = pd.read_csv(U3/"traits_exp3u.csv"); TR = TR[TR.include].set_index("id"); LG = pd.read_csv(U3/"results"/"u02_lineage.csv", encoding="utf-8-sig")
E = LG[(LG.p1 < .05 / 12) & (LG.tau > 0)].copy(); P(f"농축된 쌍: {len(E)}개 (GWAS {E.id.nunique()}개; 106개 GWAS × 12묶음 = 1,272쌍 가운데)")
CATN = {"S": "skeletal", "BI": "blood and immune", "LL": "lipid", "BR": "brain and behaviour", "U": "other"}
def p_t3r(i, g, s):
    a = U/"runs"/f"{i}__T3r{s:02d}_{g}.npz"; return a if a.exists() else W/"runs"/f"{i}__T3r{s:02d}_{g}.npz"
rows = []; D = {}
for r in E.itertuples():
    rec = dict(id=r.id, trait=TR.trait[r.id], category=CATN[TR.group[r.id]], group=r.lineage, m1_tau=r.tau, m1_p1=r.p1)
    for lab, Lr in (("all", load(U/"runs"/f"{r.id}__T3all_{r.lineage}.npz")),
                    ("res", avg([load(p_t3r(r.id, r.lineage, s)) for s in range(1, 21)]) if r.lineage in FIVE else None),
                    ("mat", load(W/"runs"/f"{r.id}__M3_{r.lineage}.npz") if r.lineage in FIVE else None)):
        if Lr is None: continue
        e, dv, ok, c1, dc = phi(Lr); rec.update({f"{lab}_cage": Lr["CAGE"][0], f"{lab}_dnase": Lr["DNASE"][0], f"{lab}_c1": c1, f"{lab}_c1_se": jse(dc), f"{lab}_phiD": e, f"{lab}_phiD_se": jse(dv), f"{lab}_valid": ok})
        D[(r.id, r.lineage, lab)] = (e, dv, ok)
    rows.append(rec)
PR = pd.DataFrame(rows); PR.to_csv(RES/"w05_pairs.csv", index=False)
P(f"균형이 정의되는 쌍(CAGE + DNase > 0), 전체 트랙: {int(PR.all_valid.sum())}개 / {len(PR)}")
z = PR.all_c1 / PR.all_c1_se; P(f"CAGE > DNase (단측 P < 0.05): {int((z > 1.645).sum())}쌍; DNase > CAGE: {int((z < -1.645).sum())}쌍")
# ---- 조직 묶음 주석마다 ----
def ann_stat(lab, groups, sel=None):
    out = {}
    for g in groups:
        ks = [(i, g, lab) for i in PR.id[PR.group == g] if (i, g, lab) in D and D[(i, g, lab)][2] and (sel is None or sel(i))]
        if not ks: continue
        m, dv = ivw([(D[k][0], D[k][1]) for k in ks]); out[g] = dict(n=len(ks), mean=m, se=jse(dv), median=float(np.median([D[k][0] for k in ks])), dv=dv)
    return out
A_all = ann_stat("all", LIN); A_res = ann_stat("res", FIVE); A_mat = ann_stat("mat", FIVE)
rows = []
for g in LIN:
    rec = dict(group=g, n=A_all[g]["n"], all_mean=A_all[g]["mean"], all_se=A_all[g]["se"], all_median=A_all[g]["median"])
    for lab, A in (("res", A_res), ("mat", A_mat)):
        if g in A: rec.update({f"{lab}_n": A[g]["n"], f"{lab}_mean": A[g]["mean"], f"{lab}_se": A[g]["se"], f"{lab}_median": A[g]["median"]})
    rows.append(rec)
AN = pd.DataFrame(rows).sort_values("all_mean"); AN.to_csv(RES/"w05_annotations.csv", index=False)
P("\n== 조직 묶음 주석마다의 균형 φD: 역분산 가중 평균 ± s.e. (중앙값), 모든 농축 쌍")
for r in AN.itertuples():
    s = f"  {r.group:23s} n {r.n:2d}  전체 트랙 {r.all_mean:5.1f}° ± {r.all_se:.1f}° ({r.all_median:.1f}°)"
    if r.group in FIVE: s += f" | 세 단위 {r.res_mean:5.1f}° ± {r.res_se:.1f}° ({r.res_median:.1f}°; n {int(r.res_n)}) | 세포형 일치 {r.mat_mean:5.1f}° ± {r.mat_se:.1f}° ({r.mat_median:.1f}°; n {int(r.mat_n)})"
    P(s)
P("주석 사이의 차이 (공동 jackknife):")
for lab, A in (("전체 트랙", A_all), ("세 단위", A_res), ("세포형 일치", A_mat)):
    for a, b in (("stromal", "liver"), ("stromal", "CNS_tissue"), ("stromal", "myeloid"), ("stromal", "lymphoid")):
        d = A[a]["mean"] - A[b]["mean"]; se = jse(A[a]["dv"] - A[b]["dv"]); P(f"  [{lab}] {a} − {b}: {d:+.1f}° ± {se:.1f}° (P = {2 * stats.norm.sf(abs(d / se)):.1e})")
# ---- 균형의 분산: 조직 묶음과 GWAS가 설명하는 몫 ----
def eta2(df, col, by):
    g = df.groupby(by)[col]; ssb = float((g.size() * (g.mean() - df[col].mean()) ** 2).sum()); sst = float(((df[col] - df[col].mean()) ** 2).sum()); k = g.ngroups; n = len(df)
    return ssb / sst, 1 - (1 - ssb / sst) * (n - 1) / (n - k)
V = PR[PR.all_valid]; P(f"\n== 균형의 분산 (전체 트랙, 균형이 정의되는 쌍 {len(V)}개)")
for by, lab in (("group", "조직 묶음"), ("id", "GWAS"), ("category", "형질 범주")):
    e, a = eta2(V, "all_phiD", by); P(f"  {lab}: η² = {e:.2f} (보정 R² = {a:.2f}; 수준 {V[by].nunique()}개)")
V5 = PR[PR.group.isin(FIVE)]
for col, lab in (("all_phiD", "전체 트랙"), ("res_phiD", "세 단위"), ("mat_phiD", "세포형 일치")):
    v = V5[V5[col.replace("phiD", "valid")] == True]; e, a = eta2(v, col, "group"); e2, a2 = eta2(v, col, "id"); P(f"  다섯 묶음의 쌍 {len(v)}개 [{lab}]: 조직 묶음 η² = {e:.2f} (보정 {a:.2f}), GWAS η² = {e2:.2f} (보정 {a2:.2f})")
# ---- 같은 GWAS 안: 여러 묶음에 농축된 GWAS에서 묶음 사이의 차이 ----
P("\n== 같은 GWAS가 두 묶음 모두에 농축된 경우의 균형 차이 (전체 트랙; GWAS마다 차이를 내고 역분산 가중 평균, 공동 jackknife)")
def within_gwas(a, b, lab="all"):
    ids = [i for i in set(PR.id[PR.group == a]) & set(PR.id[PR.group == b]) if D[(i, a, lab)][2] and D[(i, b, lab)][2]]
    if len(ids) < 3: return None
    m, dv = ivw([(D[(i, a, lab)][0] - D[(i, b, lab)][0], D[(i, a, lab)][1] - D[(i, b, lab)][1]) for i in ids]); return len(ids), m, jse(dv)
for a, b in (("stromal", "adipose"), ("stromal", "skeletal_muscle"), ("stromal", "cardiovascular"), ("stromal", "liver"), ("reproductive_endocrine", "adipose"), ("adipose", "liver"), ("lymphoid", "myeloid"), ("cardiovascular", "adipose")):
    r = within_gwas(a, b)
    if r: P(f"  {a} − {b}: GWAS {r[0]}개, {r[1]:+.1f}° ± {r[2]:.1f}° (P = {2 * stats.norm.sf(abs(r[1] / r[2])):.1e})")
# ---- 같은 주석 안: 골격 형질과 그 밖 ----
P("\n== 같은 주석 안에서 골격 형질 − 그 밖의 형질 (역분산 가중 평균의 차이, 공동 jackknife)")
SK = set(TR.index[TR.group == "S"])
for lab, name, groups in (("all", "전체 트랙", LIN), ("res", "세 단위", FIVE), ("mat", "세포형 일치", FIVE)):
    items = []
    for g in groups:
        a = ann_stat(lab, [g], lambda i: i in SK); b = ann_stat(lab, [g], lambda i: i not in SK)
        if g in a and g in b and a[g]["n"] >= 2:
            d = a[g]["mean"] - b[g]["mean"]; dv = a[g]["dv"] - b[g]["dv"]; items.append((d, dv)); P(f"  [{name}] {g:23s} 골격 {a[g]['mean']:.1f}° (n {a[g]['n']}) − 그 밖 {b[g]['mean']:.1f}° (n {b[g]['n']}) = {d:+.1f}° ± {jse(dv):.1f}°")
    if len(items) >= 2: m, dv = ivw(items); P(f"  [{name}] 묶음에 걸친 평균: {m:+.1f}° ± {jse(dv):.1f}° (P = {2 * stats.norm.sf(abs(m / jse(dv))):.1e})")
# ---- 전체 트랙과 세 단위·세포형 일치의 일치 ----
v = V5[(V5.all_valid == True) & (V5.res_valid == True) & (V5.mat_valid == True)]
P(f"\n== 다섯 묶음의 쌍 {len(v)}개에서 균형의 상관: 전체 트랙–세 단위 r = {np.corrcoef(v.all_phiD, v.res_phiD)[0, 1]:.2f}, 전체 트랙–세포형 일치 r = {np.corrcoef(v.all_phiD, v.mat_phiD)[0, 1]:.2f}")
# ---- 조직 묶음과 GWAS를 함께 넣은 모형: 두 묶음 이상에 농축된 GWAS의 쌍 (가중하지 않은 최소제곱; 한 요인을 뺐을 때 늘어나는 잔차 제곱합) ----
def fit_ss(df, col, facs):
    y = df[col].to_numpy(float); Xs = [np.ones(len(df))] + [pd.get_dummies(df[f], drop_first=True).to_numpy(float) for f in facs]; Xm = np.column_stack(Xs)
    b, *_ = np.linalg.lstsq(Xm, y, rcond=None); r = y - Xm @ b; return float(r @ r), int(np.linalg.matrix_rank(Xm))
P("\n== 조직 묶음과 GWAS를 함께 넣은 모형 (전체 트랙; 균형이 정의되고 두 묶음 이상에 농축된 GWAS)")
m = V.groupby("id").size(); V2 = V[V.id.isin(m.index[m >= 2])]; n = len(V2); sst = float(((V2.all_phiD - V2.all_phiD.mean()) ** 2).sum())
sf, kf = fit_ss(V2, "all_phiD", ["group", "id"]); sg, kg = fit_ss(V2, "all_phiD", ["id"]); si, ki = fit_ss(V2, "all_phiD", ["group"])
Fg = ((sg - sf) / (kf - kg)) / (sf / (n - kf)); Fi = ((si - sf) / (kf - ki)) / (sf / (n - kf))
P(f"  쌍 {n}개, GWAS {V2.id.nunique()}개, 조직 묶음 {V2.group.nunique()}개; 전체 R² = {1 - sf / sst:.2f}")
P(f"  조직 묶음(GWAS를 넣은 뒤): 부분 R² = {(sg - sf) / sg:.2f}, F({kf - kg}, {n - kf}) = {Fg:.1f}, P = {stats.f.sf(Fg, kf - kg, n - kf):.1e}")
P(f"  GWAS(조직 묶음을 넣은 뒤): 부분 R² = {(si - sf) / si:.2f}, F({kf - ki}, {n - kf}) = {Fi:.1f}, P = {stats.f.sf(Fi, kf - ki, n - kf):.1e}")
P("  GWAS마다 농축된 묶음의 수:", m.value_counts().sort_index().to_dict())
(RES/"w05_summary.txt").write_text("\n".join(log), encoding="utf-8")
