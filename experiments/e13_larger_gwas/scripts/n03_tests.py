"""실험 3n n03: 최신 GWAS 강건성 검정 (PROTOCOL_exp3n_2026-10-21, 2026-10-21 고정).
 계산은 3m m_common(= 3j j03의 τ*·공동 jackknife·Q·결합 대비)과 3f v03의 분류 규칙을 그대로 쓴다.
 갱신 지도 = 지도 형질 11개 중 키·BMI·조현병·양극성장애·LDL만 새 GWAS(exp3n/runs)로 바꾼 것. 나머지 6개는 3j 실행 그대로.
 N-H1 (주) 갱신 지도 Q (자유도 20) P < 0.05 → 유지.
 N-H2 (주) 뇌 형질 4개 대 키, AD 대 키 (c1, 양측, α 0.05/3): 같은 부호 & P < α → 유지 / 같은 부호 & P ≥ α → 약화 / 부호 반대 → 반대. 셋째 대비는 보고만.
 N-H3 새 GWAS 다섯 형질의 분류(3f 규칙), c1·c2. 키: DNase 1위이고 c1 < 0 (단측 P < 0.05) → "키의 DNase 쪽 위치 유지".
 N-H4 원래 대 새 GWAS 비교(기술). N-H5 계통 지도(기술). N-H6 활성도 조건화 EA(기술).
출력: results/n03_tests.txt, n03_layers.csv, n03_lineage.csv, n03_activity.csv, n03_trait_info.csv
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
EXP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXP.parent/"e12_diagnostics"/"scripts"))
import m_common as MC
from m_common import layers, prof, jse, Qtest, comb, CON, SEEDS, LAY, MAP, KO, M_TOT
NR = Path(f"{DATA}/exp3n/runs"); JR = MC.RUNS; C = Path(f"{DATA}/exp3c"); Y = Path(f"{DATA}/exp3d"); L = Path(f"{DATA}/exp3l")
NEW = ["HEIGHT", "BMI", "SCZ", "BIP", "LDL"]
OWN = MC.OWN; OWNL = {t: ("DNASE" if t == "HEIGHT" else "CAGE") for t in MAP}
LIN12 = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
def fp(p): return f"{p:.2g}"
ALPHA = 0.05 / 3
# ---------------- 층위 지도 (원래, 갱신) ----------------
PR_old = {t: prof(t, SEEDS, JR) for t in MAP}
PR_new = {t: (prof(t, SEEDS, NR) if t in NEW else PR_old[t]) for t in MAP}
Qo = Qtest(PR_old, MAP); Qn = Qtest(PR_new, MAP)
assert abs(Qo[0] - 74.3) < .05, Qo   # 주 지도 재현 (3j)
P("=== N-H1 (주): 갱신 지도 형질 간 층위 이질성 ===")
P(f"  원래 지도 Q {Qo[0]:.1f} (자유도 {Qo[1]}, P {fp(Qo[2])}, F 근사 {fp(Qo[3])})")
P(f"  갱신 지도 Q {Qn[0]:.1f} (자유도 {Qn[1]}, P {fp(Qn[2])}, F 근사 {fp(Qn[3])}) → {'형질 간 층위 이질성 유지' if Qn[2] < .05 else '미달'}")
P("\n=== N-H2 (주): 키가 들어간 사전 대비 (양측, α 0.0167) ===")
rows_c = []
for name, groups, k in CON:
    eo, so, po, _ = comb(PR_old, groups, k); en, sn, pn, _ = comb(PR_new, groups, k)
    if "height" in name:
        v = "반대" if np.sign(en) != np.sign(eo) else ("유지" if pn < ALPHA else "약화")
    else: v = "보고만 (바뀐 형질 없음)"
    P(f"  {name}: 원래 {eo:+.2f} ± {so:.2f} (P {fp(po)}) → 갱신 {en:+.2f} ± {sn:.2f} (P {fp(pn)}) → {v}"); rows_c.append(dict(contrast=name, old=eo, old_se=so, old_p=po, new=en, new_se=sn, new_p=pn, verdict=v))
# ---------------- N-H3·N-H4 형질별 ----------------
def seed_tops(t, runs): return pd.Series([max((MC.load(t, s, runs)[m][0], m) for m in LAY)[1] for s in SEEDS]).value_counts()
def classify(t, runs, PRd):
    Ld = layers(t, SEEDS, runs); est = {m: Ld[m][0] for m in LAY}; order = sorted(LAY, key=lambda m: -est[m]); X1, X2 = order[0], order[1]
    d = Ld[X1][1] - Ld[X2][1]; cx = est[X1] - est[X2]; px = float(stats.norm.sf(cx / jse(d))); top = seed_tops(t, runs); ntop = int(top.get(X1, 0))
    conf = (ntop >= .75 * len(SEEDS)) and (px < .05)
    c1, c2 = PRd[t][0]; c1se, c2se = jse(PRd[t][1][:, 0]), jse(PRd[t][1][:, 1])
    return dict(top=X1, n_top=ntop, second=X2, top_vs_second=cx, top_vs_second_p1=px, cls=(f"확정 {X1}형" if conf else f"미정 (최빈 {top.index[0]} {int(top.iloc[0])}/{len(SEEDS)})"),
                **{f"{m}_tau": est[m] for m in LAY}, **{f"{m}_se": jse(Ld[m][1]) for m in LAY}, c1=c1, c1_se=c1se, c2=c2, c2_se=c2se)
P("\n=== N-H3·N-H4: 새 GWAS 다섯 형질의 층위 (재추출 20회 평균; 원래 → 새) ===")
rows = []
for t in NEW:
    o = classify(t, JR, PR_old); n = classify(t, NR, PR_new)
    for k, r in [("old", o), ("new", n)]: rows.append(dict(trait=t, gwas=k, lineage=OWN[t], **r))
    P(f"  {KO[t]:6s} ({OWN[t]}): " + " | ".join(f"{m} {o[m + '_tau']:+.2f} → {n[m + '_tau']:+.2f} ± {n[m + '_se']:.2f}" for m in LAY))
    P(f"           c1 {o['c1']:+.2f} ± {o['c1_se']:.2f} → {n['c1']:+.2f} ± {n['c1_se']:.2f} | c2 {o['c2']:+.2f} → {n['c2']:+.2f} ± {n['c2_se']:.2f} | 분류 {o['cls']} → {n['cls']}")
h = [r for r in rows if r["trait"] == "HEIGHT" and r["gwas"] == "new"][0]; hp = float(stats.norm.cdf(h["c1"] / h["c1_se"]))
P(f"  키 판정: 1위 {h['top']}, c1 {h['c1']:+.2f} ± {h['c1_se']:.2f} (단측 P(c1 < 0) {fp(hp)}) → {'키의 DNase 쪽 위치 유지' if h['top'] == 'DNASE' and hp < .05 else '미달'}")
pd.DataFrame(rows).to_csv(EXP/"results"/"n03_layers.csv", index=False)
# ---------------- 형질 정보 (T3all 로그) ----------------
def info(t, runs):
    s = open(runs/f"{t}__T3all.log", encoding="utf-8").read()
    g = lambda pat: re.search(pat, s)
    h2 = g(r"Total Observed scale h2: ([\d.eE-]+) \(([\d.eE-]+)\)"); ic = g(r"Intercept: ([\d.eE-]+) \(([\d.eE-]+)\)"); mc = g(r"Mean Chi\^2: ([\d.eE-]+)"); lg = g(r"Lambda GC: ([\d.eE-]+)")
    return dict(h2=float(h2.group(1)), h2_se=float(h2.group(2)), Z=float(h2.group(1)) / float(h2.group(2)), mean_chi2=float(mc.group(1)), lambda_gc=float(lg.group(1)), intercept=float(ic.group(1)), intercept_se=float(ic.group(2)),
                n_snp=int(np.load(runs/f"{t}__T3all.npz")["n_snp"]))
I = pd.DataFrame([dict(trait=t, gwas=k, **info(t, r)) for t in NEW for k, r in [("old", JR), ("new", NR)]]); I.to_csv(EXP/"results"/"n03_trait_info.csv", index=False)
P("\n=== 형질 정보 (T3all 모형 로그; 원래 → 새) ===")
for t in NEW:
    o, n = I[(I.trait == t) & (I.gwas == "old")].iloc[0], I[(I.trait == t) & (I.gwas == "new")].iloc[0]
    P(f"  {KO[t]:6s} h2 {o.h2:.3f} (Z {o.Z:.1f}) → {n.h2:.3f} (Z {n.Z:.1f}) | 평균 χ² {o.mean_chi2:.2f} → {n.mean_chi2:.2f} | 절편 {o.intercept:.3f} → {n.intercept:.3f} ({n.intercept_se:.3f}) | 회귀 SNP {o.n_snp:,} → {n.n_snp:,}")
# ---------------- N-H5 계통 지도 ----------------
SDC = pd.read_csv(C/"annot_sd_common.csv", index_col=0).sd_common
def tau1(t, model, runs, name):
    z = np.load(runs/f"{t}__{model}.npz"); h2 = float(z["h2_log"]); names = [re.sub(r"L2_\d+$", "", x) for x in z["names"]]; k = names.index(name); f = SDC[name] * M_TOT / h2
    return z["coef"][k] * f, jse(z["part_delete"][:, k] * f)
S3C = pd.read_csv(EXP.parent/"e02_lineage_map"/"results"/"z04_sldsc.csv"); S3C = S3C[S3C.model == "M1"]
P("\n=== N-H5: 계통 지도 (단일 계통 주석; 고정 계통의 τ*와 12개 중 순위; 원래 → 새) ===")
lr = []
for t in NEW:
    for g in LIN12:
        v, s = tau1(t, f"M1_{g}", NR, f"L5_{g}"); lr.append(dict(trait=t, lineage=g, tau=v, se=s, p1=float(stats.norm.sf(v / s))))
    D = pd.DataFrame([r for r in lr if r["trait"] == t]).sort_values("tau", ascending=False).reset_index(drop=True); rk = int(D.index[D.lineage == OWN[t]][0]) + 1; r = D[D.lineage == OWN[t]].iloc[0]
    o = S3C[(S3C.trait == t)].sort_values("tau_star", ascending=False).reset_index(drop=True); ork = int(o.index[o.lineage == OWN[t]][0]) + 1; oo = o[o.lineage == OWN[t]].iloc[0]
    P(f"  {KO[t]:6s} {OWN[t]}: τ* {oo.tau_star:+.2f} ({ork}위) → {r.tau:+.2f} ± {r.se:.2f} (단측 P {fp(r.p1)}, {rk}위) | 새 1위 {D.lineage[0]} {D.tau[0]:+.2f}")
pd.DataFrame(lr).to_csv(EXP/"results"/"n03_lineage.csv", index=False)
# ---------------- N-H6 활성도 조건화 ----------------
SDA = pd.concat([pd.read_csv(Y/"annot_sd_common.csv", index_col=0).sd_common, pd.read_csv(L/"annot_sd_common.csv", index_col=0).sd_common]); SDA = SDA[~SDA.index.duplicated()]
def load_ea(t, runs, model="EA"):
    z = np.load(runs/f"{t}__{model}.npz"); h2 = float(z["h2_log"]); res = {}
    for k, n in enumerate(re.sub(r"L2_\d+$", "", x) for x in z["names"]):
        if n not in SDA.index: continue
        kind = "E" if n.startswith("L5_") else "A"; m = [x for x in LAY if n.endswith("_" + x)][0]; f = SDA[n] * M_TOT / h2
        res[(kind, m)] = (z["coef"][k] * f, jse(z["part_delete"][:, k] * f))
    return res
V3L = pd.read_csv(EXP.parent/"e11_activity"/"results"/"l05_values.csv")
P("\n=== N-H6: 활성도 조건화 (EA, 전체 트랙; 자기 층위 효과 τ*; 원래 → 새) ===")
ar = []
for t in NEW:
    m = OWNL[t]; E0 = load_ea(t, NR, "T3all")[("E", m)]; EA = load_ea(t, NR)
    e, se = EA[("E", m)]; p1 = float(stats.norm.sf(e / se))
    o = V3L[(V3L.trait == t) & (V3L.model == "EA") & (V3L.kind == "E") & (V3L.layer == m)].iloc[0]; o0 = V3L[(V3L.trait == t) & (V3L.model == "E") & (V3L.kind == "E") & (V3L.layer == m)].iloc[0]
    ar.append(dict(trait=t, layer=m, E=E0[0], EA=e, EA_se=se, p1=p1, retention=e / E0[0], old_E=o0.tau, old_EA=o.tau, old_EA_se=o.se, old_retention=o.tau / o0.tau,
                   **{f"act_{x}": EA[("A", x)][0] for x in LAY}, **{f"act_{x}_se": EA[("A", x)][1] for x in LAY}))
    P(f"  {KO[t]:6s} {m}: 원래 E {o0.tau:+.2f} → EA {o.tau:+.2f} ± {o.se:.2f} ({o.tau / o0.tau * 100:.0f} %) | 새 E {E0[0]:+.2f} → EA {e:+.2f} ± {se:.2f} (단측 P {fp(p1)}, {e / E0[0] * 100:.0f} %)"
      f" | 새 활성도 " + ", ".join(f"{x} {EA[('A', x)][0]:+.2f}" for x in LAY))
pd.DataFrame(ar).to_csv(EXP/"results"/"n03_activity.csv", index=False)
pd.DataFrame(rows_c).to_csv(EXP/"results"/"n03_contrasts.csv", index=False)
(EXP/"results"/"n03_tests.txt").write_text("\n".join(log), encoding="utf-8")
