"""실험 3v v02: 형질의 두 지표(T = 조직 묶음의 폭, L = 층위의 폭)와 유형 (PROTOCOL_exp3v 3–4절).
 단계 (python -I -X utf8 v02_types.py main | all):
   main = 주 GWAS 11형질(3d y04_B.csv, 792칸)만. 새 적합이 필요 없다.
   all  = main + 106개 GWAS(exp3v/runs/<id>.npz, 7,632칸)의 분석 A–E.
 정의(프로토콜 3절): 유의한 주석 = 분석 집합의 모든 칸에 걸친 BH FDR < 0.05(단측). K = 유의한 주석 수, T = 유의한 주석이 있는 조직 묶음 수, L = K / T.
   유형: 집중형(T ≤ 2, L < 3), 층위 확장형(T ≤ 2, L ≥ 3), 조직 확장형(T ≥ 3, L < 3), 조직·층위 확장형(T ≥ 3, L ≥ 3), 신호 없음(K = 0).
 출력: results/v02_main11.csv, v02_main11.txt, (all) v02_matrix106.csv, v02_gwas106.csv, v02_pairs106.csv, v02_types.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; RES = EXP/"results"; RES.mkdir(exist_ok=True); ROOT = EXP.parent
DATA = Path(DATA); B = DATA/"exp3"/"GRCh38"; Y = DATA/"exp3d"; V = DATA/"exp3v"
LIN = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
MODS = ["CAGE", "DNASE", "H3K27ac", "H3K4me1", "H3K4me3", "RNA"]
ANN = [f"L5_{g}_{m}" for g in LIN for m in MODS]
ORD = ["EA", "SCZ", "BIP", "NEURO", "BMI", "AD", "MONO", "IBD", "RA", "LDL", "HEIGHT"]
TYPES = ["focal", "layer-spanning", "tissue-spanning", "tissue- and layer-spanning", "none"]
KO = {"focal": "집중형", "layer-spanning": "층위 확장형", "tissue-spanning": "조직 확장형", "tissue- and layer-spanning": "조직·층위 확장형", "none": "신호 없음"}
T_CUT, L_CUT = 3, 3.0
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
def bh(p):
    p = np.asarray(p); o = np.argsort(p); q = p[o] * len(p) / np.arange(1, len(p) + 1); q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty_like(q); out[o] = np.minimum(q, 1); return out
def typ(K, T, L):
    if K == 0: return "none"
    return ("tissue- and layer-spanning" if L >= L_CUT else "tissue-spanning") if T >= T_CUT else ("layer-spanning" if L >= L_CUT else "focal")
def index(D, key, thr):
    """D: 칸 표(key, unit, mod, q). 형질마다 K, T, L, 유형, 묶음별 층위."""
    rows = []
    for t, d in D.groupby(key, sort=False):
        s = d[d.q < thr]; K = len(s); T = s.unit.nunique(); L = K / T if T else np.nan
        per = s.groupby("unit").size(); top = d.loc[d.tau_star.idxmax()]
        rows.append({key: t, "K": K, "T": T, "L": L, "max_layers": int(per.max()) if T else 0, "n_layers_distinct": s["mod"].nunique(), "type": typ(K, T, L),
                     "top_annot": f"{top.unit}|{top['mod']}", "top_tau": float(top.tau_star),
                     "where": "; ".join(f"{g}: {', '.join(m for m in MODS if m in set(s[s.unit == g]['mod']))}" for g in LIN if g in set(s.unit))})
    return pd.DataFrame(rows).set_index(key)

# ------------------------------------------------------------------ 주 GWAS 11형질
D1 = pd.read_csv(ROOT/"e03_single_layer_map"/"results"/"y04_B.csv"); assert len(D1) == 792 and set(D1.trait) == set(ORD)
assert np.allclose(bh(D1.tau_p1), D1.q_bh, rtol=1e-12, atol=0); D1["q"] = D1.q_bh
M = index(D1, "trait", .05).loc[ORD]
for thr, tag in ((.10, "q10"), (.01, "q01")):
    s = index(D1, "trait", thr).loc[ORD]; M[[f"K_{tag}", f"T_{tag}", f"L_{tag}", f"type_{tag}"]] = s[["K", "T", "L", "type"]]
M.to_csv(RES/"v02_main11.csv", encoding="utf-8-sig")
P("=== 주 GWAS 11형질 (792칸, FDR < 0.05) ===")
P(f"{'형질':7s} {'K':>3s} {'T':>3s} {'L':>5s}  {'유형':14s} | FDR<0.10: K T L 유형 | FDR<0.01: K T L 유형 | 어디에")
for t, r in M.iterrows():
    P(f"{t:7s} {r.K:3d} {r['T']:3d} {r.L:5.2f}  {KO[r.type]:14s} | {r.K_q10:2d} {r.T_q10:2d} {r.L_q10:4.2f} {KO[r.type_q10]:10s} | {r.K_q01:2d} {r.T_q01:2d} {r.L_q01:4.2f} {KO[r.type_q01]:10s} | {r['where']}")
P("유형별 형질 수: " + ", ".join(f"{KO[k]} {int((M.type == k).sum())}" for k in TYPES))
P("역치를 바꾸면 유형이 달라지는 형질: FDR<0.10 " + (", ".join(M.index[M.type != M.type_q10]) or "없음") + " / FDR<0.01 " + (", ".join(M.index[M.type != M.type_q01]) or "없음"))
if STAGE == "main":
    (RES/"v02_main11.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit(0)
(RES/"v02_main11.txt").write_text("\n".join(log), encoding="utf-8")

# ------------------------------------------------------------------ 106개 GWAS
TR = pd.read_csv(ROOT/"e20_gwas_panel"/"traits_exp3u.csv"); TR = TR[TR.include].reset_index(drop=True); assert len(TR) == 106
SD = pd.read_csv(Y/"annot_sd_common.csv", index_col=0).sd_common
M_TOT = sum(float(open(B/"baselineLD_v2.2"/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
rows = []; hz = {}
for r in TR.itertuples():
    z = np.load(V/"runs"/f"{r.id}.npz", allow_pickle=False); assert list(z["ann"]) == ANN
    h2 = np.round(z["tot"], 4); f = np.array([SD[a] for a in ANN]) * M_TOT / h2; hz[r.id] = float(np.median(z["tot"] / z["tot_se"]))
    for j, a in enumerate(ANN):
        g, m = LIN[j // 6], MODS[j % 6]; assert a == f"L5_{g}_{m}"
        rows.append(dict(id=r.id, file=r.file, group=r.group, status=r.status, unit=g, mod=m, annot=a, h2=float(h2[j]), tau_star=float(z["coef"][j] * f[j]), tau_star_se=float(z["coef_se"][j] * f[j]),
                         z=float(z["coef"][j] / z["coef_se"][j])))
X = pd.DataFrame(rows); assert len(X) == 7632 and np.isfinite(X.tau_star).all() and (X.tau_star_se > 0).all()
X["p1"] = stats.norm.sf(X.z); X["q"] = bh(X.p1); X.to_csv(RES/"v02_matrix106.csv", index=False, encoding="utf-8-sig")
G = index(X, "id", .05)
for thr, tag in ((.10, "q10"), (.01, "q01")):
    s = index(X, "id", thr); G[[f"K_{tag}", f"T_{tag}", f"L_{tag}", f"type_{tag}"]] = s[["K", "T", "L", "type"]]
G = TR.set_index("id")[["file", "trait", "disease", "group", "status", "discovery_trait", "h2_z"]].rename(columns={"h2_z": "h2_z_source"}).join(G); G["h2_z"] = pd.Series(hz)
G.to_csv(RES/"v02_gwas106.csv", encoding="utf-8-sig")
P(f"\n=== 106개 GWAS (7,632칸): FDR < 0.05인 칸 {int((X.q < .05).sum())}개, 하나 이상 있는 GWAS {int((G.K > 0).sum())}개 ===")
P(f"층위별 유의한 칸: " + ", ".join(f"{m} {int(((X.q < .05) & (X['mod'] == m)).sum())}" for m in MODS))
P(f"조직 묶음별 유의한 칸: " + ", ".join(f"{g} {int(((X.q < .05) & (X.unit == g)).sum())}" for g in LIN))
# A 분포
P("\n--- A. 유형별 GWAS 수 ---")
P(f"{'유형':16s} 전체  " + "  ".join(f"{g:>3s}" for g in ["S", "BI", "LL", "BR", "U"]) + "  새 형질만")
for k in TYPES:
    P(f"{KO[k]:16s} {int((G.type == k).sum()):4d}  " + "  ".join(f"{int(((G.type == k) & (G.group == g)).sum()):3d}" for g in ["S", "BI", "LL", "BR", "U"]) + f"  {int(((G.type == k) & (G.status == 'new')).sum()):4d}")
P("지도의 자리 (T × L 구간; L 구간 = [1,2), [2,3), [3,4), [4,6]) — GWAS 수")
Gk = G[G.K > 0].copy(); Gk["Lbin"] = pd.cut(Gk.L, [1, 2, 3, 4, 6.01], right=False, labels=["1–<2", "2–<3", "3–<4", "4–6"])
tab = pd.crosstab(Gk.Lbin, Gk["T"]).reindex(index=["4–6", "3–<4", "2–<3", "1–<2"], fill_value=0); P(tab.to_string())
for k in TYPES[:4]:
    d = G[G.type == k].sort_values(["T", "L"], ascending=False); P(f"\n[{KO[k]}] {len(d)}개: " + "; ".join(f"{r.file.split('.', 1)[1] if '.' in r.file else r.file} (T {r.T}, L {r.L:.1f}, K {r.K})" for r in d.itertuples()))
# B 검정력
P("\n--- B. 검정력 (h² z = 72개 모형의 총 h² ÷ 표준오차의 중앙값) ---")
for c in ["K", "T"]:
    rho, p = stats.spearmanr(G.h2_z, G[c]); P(f"Spearman(h² z, {c}) = {rho:+.2f} (P = {p:.1e}; 106개)")
rho, p = stats.spearmanr(Gk.h2_z, Gk.L); P(f"Spearman(h² z, L) = {rho:+.2f} (P = {p:.1e}; K > 0인 {len(Gk)}개)")
P("유형별 h² z의 중앙값(범위): " + "; ".join(f"{KO[k]} {G[G.type == k].h2_z.median():.1f} ({G[G.type == k].h2_z.min():.1f}–{G[G.type == k].h2_z.max():.1f}, {int((G.type == k).sum())}개)" for k in TYPES if (G.type == k).any()))
med = G.h2_z.median(); Hh = G[G.h2_z >= med]
P(f"h² z가 중앙값({med:.1f}) 이상인 {len(Hh)}개: " + ", ".join(f"{KO[k]} {int((Hh.type == k).sum())}" for k in TYPES))
Lo = G[G.h2_z < med]; P(f"h² z가 중앙값 미만인 {len(Lo)}개: " + ", ".join(f"{KO[k]} {int((Lo.type == k).sum())}" for k in TYPES))
# C 같은 형질
P("\n--- C. 주 GWAS와 형질이 같은 GWAS ---")
for r in G[G.status != "new"].itertuples():
    m = M.loc[r.discovery_trait] if r.discovery_trait in M.index else None
    P(f"{r.file:48s} K {r.K:2d} T {r.T:2d} L {r.L:4.2f} {KO[r.type]:10s} | 주 GWAS {r.discovery_trait}: " + (f"K {m.K} T {m['T']} L {m.L:.2f} {KO[m.type]}" if m is not None else "행렬에 없는 형질(검증 형질)") + f" | {r.where}")
# D 조직 묶음의 층위 구성
P("\n--- D. 조직 묶음의 층위 구성 (유의한 주석이 하나 이상 있는 형질–조직 묶음 쌍) ---")
S5 = X[X.q < .05]; PR = S5.groupby(["id", "unit"]).agg(n_layers=("mod", "size"), layers=("mod", lambda s: ",".join(m for m in MODS if m in set(s)))).reset_index()
PR = PR.merge(G[["file", "group", "status", "T", "type"]], left_on="id", right_index=True); PR.to_csv(RES/"v02_pairs106.csv", index=False, encoding="utf-8-sig")
P(f"쌍 {len(PR)}개 (GWAS {PR.id.nunique()}개)")
P(f"{'조직 묶음':24s} 쌍  평균 층위 수  " + " ".join(f"{m:>8s}" for m in MODS) + "   (층위별: 그 층위가 유의한 쌍의 비율)")
for g in LIN:
    d = PR[PR.unit == g]
    if len(d): P(f"{g:24s} {len(d):3d}  {d.n_layers.mean():6.2f}       " + " ".join(f"{np.mean([m in s.split(',') for s in d.layers]):8.2f}" for m in MODS))
def eta2(d, by):
    gm = d.n_layers.mean(); sst = ((d.n_layers - gm) ** 2).sum(); ssb = sum(len(x) * (x.n_layers.mean() - gm) ** 2 for _, x in d.groupby(by)); k = d[by].nunique(); n = len(d)
    r2 = ssb / sst; return r2, 1 - (1 - r2) * (n - 1) / (n - k), k, n
for by, nm in (("unit", "조직 묶음"), ("id", "형질")):
    r2, adj, k, n = eta2(PR, by); P(f"층위 수의 분산 가운데 {nm}으로 설명되는 몫: η² = {r2:.2f}, 보정 R² = {adj:.2f} (수준 {k}개, 쌍 {n}개)")
R2 = PR[PR.groupby("id").id.transform("size") >= 2]
for by, nm in (("unit", "조직 묶음"), ("id", "형질")):
    r2, adj, k, n = eta2(R2, by); P(f"  쌍이 둘 이상인 형질만: {nm} η² = {r2:.2f}, 보정 R² = {adj:.2f} (수준 {k}개, 쌍 {n}개)")
# E 민감도
P("\n--- E. 민감도 ---")
for tag, nm in (("q10", "FDR < 0.10"), ("q01", "FDR < 0.01")):
    P(f"{nm}: " + ", ".join(f"{KO[k]} {int((G['type_' + tag] == k).sum())}" for k in TYPES) + f" | 유형이 바뀐 GWAS {int((G.type != G['type_' + tag]).sum())}개" +
      f" (신호 없음과의 오감을 빼면 {int(((G.type != G['type_' + tag]) & (G.type != 'none') & (G['type_' + tag] != 'none')).sum())}개)")
    P("   " + pd.crosstab(G.type, G["type_" + tag]).reindex(index=[k for k in TYPES if (G.type == k).any()], columns=[k for k in TYPES if (G["type_" + tag] == k).any()], fill_value=0).rename(index=KO, columns=KO).to_string().replace("\n", "\n   "))
(RES/"v02_types.txt").write_text("\n".join(log), encoding="utf-8")
