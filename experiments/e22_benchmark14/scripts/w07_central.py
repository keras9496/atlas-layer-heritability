"""실험 3w w07: 균형이 조직 묶음 주석과 GWAS에 따라 얼마나 다른지, 조직 묶음의 차이가 트랙 수와 세포형을 맞춰도 남는지의 점검 (사후 분석).
 새 적합은 없다. w05와 같은 쌍(농축된 GWAS–조직 묶음 쌍 288개)과 같은 적합 결과를 읽는다.
 (가) 조직 묶음 주석 평균 사이의 이질성 Q(공동 jackknife): 전체 트랙 12묶음, 세 단위·세포형 일치 다섯 묶음, 그리고 stromal을 뺀 경우.
 (나) 같은 GWAS 안의 묶음 사이 차이를 세 방식으로(다섯 묶음 사이).
 (다) 조직 묶음 + GWAS 두 요인 모형을 다섯 묶음에서 세 방식으로.
 (라) GWAS의 몫에서 추정 오차를 뺀 값: 조직 묶음 평균을 뺀 잔차의 GWAS 평균, 그 분산에서 jackknife 분산의 평균을 뺀 값(모멘트 추정).
 (마) 형질 범주마다 잔차의 평균(어느 형질이 GWAS의 몫을 만드는가).
실행: python -I -X utf8 w07_central.py
출력: results/w07_summary.txt, w07_gwas_resid.csv
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
    return e, dv, bool(c[0] + d[0] > 0)
def ivw(items):
    w = np.array([1 / jse(d) ** 2 for _, d in items]); w = w / w.sum(); return float(sum(wi * e for wi, (e, _) in zip(w, items))), sum(wi * d for wi, (_, d) in zip(w, items))
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
TR = pd.read_csv(U3/"traits_exp3u.csv"); TR = TR[TR.include].set_index("id"); LG = pd.read_csv(U3/"results"/"u02_lineage.csv", encoding="utf-8-sig")
E = LG[(LG.p1 < .05 / 12) & (LG.tau > 0)].copy(); assert len(E) == 288
CATN = {"S": "skeletal", "BI": "blood and immune", "LL": "lipid", "BR": "brain and behaviour", "U": "other"}
def p_t3r(i, g, s):
    a = U/"runs"/f"{i}__T3r{s:02d}_{g}.npz"; return a if a.exists() else W/"runs"/f"{i}__T3r{s:02d}_{g}.npz"
D = {}; PAIRS = []
for r in E.itertuples():
    PAIRS.append((r.id, r.lineage))
    for lab, Lr in (("all", load(U/"runs"/f"{r.id}__T3all_{r.lineage}.npz")),
                    ("res", avg([load(p_t3r(r.id, r.lineage, s)) for s in range(1, 21)]) if r.lineage in FIVE else None),
                    ("mat", load(W/"runs"/f"{r.id}__M3_{r.lineage}.npz") if r.lineage in FIVE else None)):
        if Lr is not None: D[(r.id, r.lineage, lab)] = phi(Lr)
REF = pd.read_csv(RES/"w05_pairs.csv"); chk = max(abs(D[(r.id, r.group, "all")][0] - r.all_phiD) for r in REF.itertuples()); assert chk < 1e-9, chk
NAME = {"all": "전체 트랙", "res": "세 단위", "mat": "세포형 일치"}
def valid(lab, groups): return [(i, g) for i, g in PAIRS if g in groups and (i, g, lab) in D and D[(i, g, lab)][2]]
def gmeans(lab, groups, pairs=None):
    pairs = valid(lab, groups) if pairs is None else pairs; out = {}
    for g in groups:
        ks = [(i, g) for i, gg in pairs if gg == g]
        if ks: m, dv = ivw([D[(i, g, lab)][:2] for i, g in ks]); out[g] = (m, dv, len(ks))
    return out
def Qtest(means, dvs):
    G = len(means); C = np.hstack([-np.ones((G - 1, 1)), np.eye(G - 1)]); d = C @ np.asarray(means); Dv = np.asarray(dvs) @ C.T; B = len(Dv); Dc = Dv - Dv.mean(0); V = (B - 1) / B * Dc.T @ Dc
    Q = float(d @ np.linalg.solve(V, d)); p = G - 1; nu = B - 1; Fs = Q * (nu - p + 1) / (nu * p); return Q, p, float(stats.chi2.sf(Q, p)), float(stats.f.sf(Fs, p, nu - p + 1))
# ---- (가) 조직 묶음 주석 평균 사이의 이질성 ----
P("== (가) 조직 묶음 주석 평균(역분산 가중) 사이의 이질성 Q — 공동 jackknife, 자유도 = 묶음 수 − 1")
def show_Q(lab, groups, note):
    A = gmeans(lab, groups); gs = [g for g in groups if g in A]; Q, p, pc, pf = Qtest([A[g][0] for g in gs], np.column_stack([A[g][1] for g in gs]))
    rng = (min(A[g][0] for g in gs), max(A[g][0] for g in gs)); P(f"  [{NAME[lab]}] {note}: 묶음 {len(gs)}개, 평균의 범위 {rng[0]:.1f}°–{rng[1]:.1f}°, Q = {Q:.1f} (자유도 {p}), P = {pc:.1e} (F 근사 {pf:.1e})")
show_Q("all", LIN, "12묶음"); show_Q("all", [g for g in LIN if g != "stromal"], "stromal을 뺀 11묶음"); show_Q("all", [g for g in LIN if g not in ("stromal", "reproductive_endocrine")], "stromal, 생식·내분비를 뺀 10묶음")
NON = [g for g in FIVE if g != "stromal"]
for lab in ("all", "res", "mat"): show_Q(lab, FIVE, "다섯 묶음"); show_Q(lab, NON, "stromal을 뺀 네 묶음")
P("  stromal을 뺀 네 묶음 사이의 짝 차이 (공동 jackknife):")
for lab in ("all", "res", "mat"):
    A = gmeans(lab, NON); s = []
    for a in range(4):
        for b in range(a + 1, 4):
            d = A[NON[b]][0] - A[NON[a]][0]; se = jse(A[NON[b]][1] - A[NON[a]][1]); s.append(f"{NON[b]} − {NON[a]} {d:+.1f}° ± {se:.1f}°")
    P(f"  [{NAME[lab]}] " + "; ".join(s))
# ---- (나) 같은 GWAS 안의 차이, 세 방식 ----
P("\n== (나) 같은 GWAS가 두 묶음 모두에 농축된 경우의 균형 차이 — 다섯 묶음 사이, 세 방식 (GWAS마다 차이, 역분산 가중 평균, 공동 jackknife)")
def within(a, b, lab):
    ia = {i for i, g in valid(lab, [a])}; ib = {i for i, g in valid(lab, [b])}; ids = sorted(ia & ib)
    if len(ids) < 3: return None
    m, dv = ivw([(D[(i, a, lab)][0] - D[(i, b, lab)][0], D[(i, a, lab)][1] - D[(i, b, lab)][1]) for i in ids]); return len(ids), m, jse(dv)
for a, b in (("stromal", "liver"), ("stromal", "CNS_tissue"), ("stromal", "lymphoid"), ("stromal", "myeloid"), ("lymphoid", "myeloid"), ("lymphoid", "liver"), ("myeloid", "liver"), ("CNS_tissue", "liver")):
    s = []
    for lab in ("all", "res", "mat"):
        r = within(a, b, lab); s.append(f"{NAME[lab]} {r[1]:+.1f}° ± {r[2]:.1f}° (GWAS {r[0]}개)" if r else f"{NAME[lab]} — (GWAS 3개 미만)")
    P(f"  {a} − {b}: " + " | ".join(s))
# stromal 대 나머지 네 묶음: GWAS마다 (stromal − 그 GWAS의 다른 묶음 평균)
for lab in ("all", "res", "mat"):
    v = valid(lab, FIVE); items = []
    for i in sorted({i for i, g in v if g == "stromal"}):
        oth = [g for j, g in v if j == i and g != "stromal"]
        if oth: items.append((D[(i, "stromal", lab)][0] - np.mean([D[(i, g, lab)][0] for g in oth]), D[(i, "stromal", lab)][1] - np.mean([D[(i, g, lab)][1] for g in oth], axis=0)))
    m, dv = ivw(items); P(f"  stromal − 같은 GWAS의 다른 네 묶음 평균 [{NAME[lab]}]: GWAS {len(items)}개, {m:+.1f}° ± {jse(dv):.1f}°")
# ---- (다) 두 요인 모형, 다섯 묶음, 세 방식 ----
P("\n== (다) 조직 묶음 + GWAS 두 요인 모형 — 다섯 묶음의 쌍, 세 방식 (가중하지 않은 최소제곱; 세 방식 모두 균형이 정의되고 다섯 묶음 가운데 둘 이상에 농축된 GWAS)")
def fit_ss(df, col, facs):
    y = df[col].to_numpy(float); Xm = np.column_stack([np.ones(len(df))] + [pd.get_dummies(df[f], drop_first=True).to_numpy(float) for f in facs]); b, *_ = np.linalg.lstsq(Xm, y, rcond=None); r = y - Xm @ b; return float(r @ r), int(np.linalg.matrix_rank(Xm))
def twoway(df, col, lab):
    n = len(df); sst = float(((df[col] - df[col].mean()) ** 2).sum()); sf, kf = fit_ss(df, col, ["group", "id"]); sg, kg = fit_ss(df, col, ["id"]); si, ki = fit_ss(df, col, ["group"])
    Fg = ((sg - sf) / (kf - kg)) / (sf / (n - kf)); Fi = ((si - sf) / (kf - ki)) / (sf / (n - kf))
    P(f"  [{lab}] 쌍 {n}개, GWAS {df.id.nunique()}개, 묶음 {df.group.nunique()}개; 전체 R² = {1 - sf / sst:.2f}; 조직 묶음 부분 R² = {(sg - sf) / sg:.2f} (F({kf - kg}, {n - kf}) = {Fg:.1f}, P = {stats.f.sf(Fg, kf - kg, n - kf):.1e}); GWAS 부분 R² = {(si - sf) / si:.2f} (F({kf - ki}, {n - kf}) = {Fi:.1f}, P = {stats.f.sf(Fi, kf - ki, n - kf):.1e})")
v3 = sorted(set(valid("all", FIVE)) & set(valid("res", FIVE)) & set(valid("mat", FIVE))); T = pd.DataFrame([dict(id=i, group=g, all=D[(i, g, "all")][0], res=D[(i, g, "res")][0], mat=D[(i, g, "mat")][0]) for i, g in v3])
m = T.groupby("id").size(); T2 = T[T.id.isin(m.index[m >= 2])]
for lab in ("all", "res", "mat"): twoway(T2, lab, NAME[lab])
T12 = pd.DataFrame([dict(id=i, group=g, all=D[(i, g, "all")][0]) for i, g in valid("all", LIN)]); m = T12.groupby("id").size(); twoway(T12[T12.id.isin(m.index[m >= 2])], "all", "전체 트랙, 12묶음 (w05와 같아야 한다)")
twoway(T12[T12.id.isin(m.index[m >= 2]) & (T12.group != "stromal")], "all", "전체 트랙, stromal을 뺀 11묶음")
# ---- (라) GWAS의 몫에서 추정 오차를 뺀 값 ----
P("\n== (라) 조직 묶음 평균을 뺀 잔차의 GWAS 평균: 분산에서 jackknife 분산의 평균을 뺀 값 (GWAS 사이의 참 차이의 모멘트 추정; 가중하지 않음)")
def resid(lab, groups, minpairs):
    v = valid(lab, groups); A = gmeans(lab, groups, v); rows = []
    for i in sorted({i for i, g in v}):
        gs = [g for j, g in v if j == i]
        if len(gs) < minpairs: continue
        e = float(np.mean([D[(i, g, lab)][0] - A[g][0] for g in gs])); dv = np.mean([D[(i, g, lab)][1] - A[g][1] for g in gs], axis=0); rows.append(dict(id=i, trait=TR.trait[i], category=CATN[TR.group[i]], n_groups=len(gs), resid=e, se=jse(dv), lab=lab, dv=dv))
    return pd.DataFrame(rows)
def moments(R, note):
    s2 = float(R.resid.var(ddof=1)); n2 = float((R.se ** 2).mean()); t2 = s2 - n2; chi = float(((R.resid / R.se) ** 2).sum())
    P(f"  [{note}] GWAS {len(R)}개: 잔차 평균의 표준편차 {np.sqrt(s2):.1f}°, 추정 오차(jackknife s.e.의 제곱평균 제곱근) {np.sqrt(n2):.1f}°, 오차를 뺀 GWAS 사이 표준편차 {np.sqrt(max(t2, 0)):.1f}° (분산의 {max(t2, 0) / s2 * 100:.0f}%); Σz² = {chi:.0f} (자유도 약 {len(R)}, P = {stats.chi2.sf(chi, len(R)):.1e}; GWAS 사이의 공분산은 무시)")
keep = []
for lab, groups, note in (("all", LIN, "전체 트랙, 12묶음"), ("all", FIVE, "전체 트랙, 다섯 묶음"), ("res", FIVE, "세 단위, 다섯 묶음"), ("mat", FIVE, "세포형 일치, 다섯 묶음")):
    for mp in (1, 2):
        R = resid(lab, groups, mp); moments(R, f"{note}, 농축된 묶음 {mp}개 이상")
        if mp == 1: R2 = R.drop(columns="dv").copy(); R2["set"] = note; keep.append(R2)
pd.concat(keep).to_csv(RES/"w07_gwas_resid.csv", index=False)
# 쌍 수준: 전체 분산과 추정 오차
for lab, groups in (("all", LIN), ("res", FIVE), ("mat", FIVE)):
    v = valid(lab, groups); e = np.array([D[(i, g, lab)][0] for i, g in v]); s = np.array([jse(D[(i, g, lab)][1]) for i, g in v]); P(f"  쌍 수준 [{NAME[lab]}]: 쌍 {len(v)}개, 균형의 표준편차 {e.std(ddof=1):.1f}°, 쌍의 s.e. 중앙값 {np.median(s):.1f}° (제곱평균 제곱근 {np.sqrt((s ** 2).mean()):.1f}°)")
# ---- (마) 형질 범주마다 잔차 ----
P("\n== (마) 형질 범주마다: 조직 묶음 평균을 뺀 잔차의 역분산 가중 평균 (공동 jackknife)")
for lab, groups, note in (("all", LIN, "전체 트랙, 12묶음"), ("res", FIVE, "세 단위, 다섯 묶음"), ("mat", FIVE, "세포형 일치, 다섯 묶음")):
    v = valid(lab, groups); A = gmeans(lab, groups, v); s = []; cm = {}
    for c in CATN.values():
        items = [(D[(i, g, lab)][0] - A[g][0], D[(i, g, lab)][1] - A[g][1]) for i, g in v if CATN[TR.group[i]] == c]
        if len(items) >= 3: mm, dv = ivw(items); cm[c] = (mm, dv); s.append(f"{c} {mm:+.1f}° ± {jse(dv):.1f}° (쌍 {len(items)}개)")
    cs = list(cm); Q, p, pc, pf = Qtest([cm[c][0] for c in cs], np.column_stack([cm[c][1] for c in cs])); P(f"  [{note}] " + "; ".join(s) + f" — 범주 사이 Q = {Q:.1f} (자유도 {p}), P = {pc:.1e}")
(RES/"w07_summary.txt").write_text("\n".join(log), encoding="utf-8")
