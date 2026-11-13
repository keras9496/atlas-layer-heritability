"""실험 3w w04: 벤치마크 14형질의 표 (PROTOCOL_exp3w 2절, 5절). 적합은 하지 않고, 3w에서 새로 적합한 것과 기존 실행 파일(3j, 3k, 3l, 3n, 3o, 3q, 3r, 3t, 3d)을 읽어 한 기준으로 정리한다.
 τ* = 계수 × sd(공통 SNP) × M_5_50 / h²(로그의 반올림 값). 정렬 블록 실행의 표준오차는 delete 값 200개의 jackknife. 72칸 지도만 ldsc 기본 블록(3d, 3v와 같다).
 여러 형질에 걸친 유의성은 14형질의 Benjamini–Hochberg (72칸 지도는 14 × 72칸, 조직 묶음 지도는 14 × 12칸).
 균형 φD = atan2(DNase τ*, CAGE τ*) (도). delete 값의 각도는 추정값 둘레 ±180° 안에 놓는다(3u u03과 같다). CAGE + DNase ≤ 0이면 정의하지 않는다.
실행: python -I -X utf8 w04_bench.py
출력: results/w04_*.csv, w04_summary.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
EXP = Path(__file__).resolve().parents[1]; RES = EXP/"results"; DATA = Path(DATA)
X = DATA/"exp3"; C3 = DATA/"exp3c"; Y = DATA/"exp3d"; E3 = DATA/"exp3e"; F3 = DATA/"exp3f"; I3 = DATA/"exp3i"; J3 = DATA/"exp3j"; K3 = DATA/"exp3k"; L3 = DATA/"exp3l"
N3 = DATA/"exp3n"; O3 = DATA/"exp3o"; Q3 = DATA/"exp3q"; R3 = DATA/"exp3r"; T3 = DATA/"exp3t"; W = DATA/"exp3w"
LIN = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
MODS6 = ["CAGE", "DNASE", "H3K27ac", "H3K4me1", "H3K4me3", "RNA"]; LAY = ["CAGE", "DNASE", "H3K27ac"]
OWN = {"EA": "CNS_tissue", "N_SCZ": "CNS_tissue", "N_BIP": "CNS_tissue", "NEURO": "CNS_tissue", "N_BMI": "CNS_tissue", "N_AD": "myeloid", "MONO": "myeloid", "NEUT": "myeloid",
       "IBD": "lymphoid", "RA": "lymphoid", "N_LDL": "liver", "TC": "liver", "N_HEIGHT": "stromal", "BMD": "stromal"}
B14 = list(OWN); NEW6 = [t for t in B14 if t.startswith("N_")]
CAT = {"EA": "brain", "N_SCZ": "brain", "N_BIP": "brain", "NEURO": "brain", "N_BMI": "brain", "N_AD": "myeloid", "MONO": "myeloid", "NEUT": "myeloid", "IBD": "lymphoid", "RA": "lymphoid",
       "N_LDL": "liver", "TC": "liver", "N_HEIGHT": "skeletal", "BMD": "skeletal"}
EN = {"EA": "Educational attainment", "N_SCZ": "Schizophrenia", "N_BIP": "Bipolar disorder", "NEURO": "Neuroticism", "N_BMI": "BMI", "N_AD": "Alzheimer's disease", "MONO": "Monocyte count",
      "NEUT": "Neutrophil count", "IBD": "IBD", "RA": "Rheumatoid arthritis", "N_LDL": "LDL cholesterol", "TC": "Total cholesterol", "N_HEIGHT": "Height", "BMD": "Bone mineral density"}
GWAS = {"EA": "Okbay 2016", "N_SCZ": "Trubetskoy 2022", "N_BIP": "O'Connell 2025", "NEURO": "Okbay 2016", "N_BMI": "Yengo 2018", "N_AD": "2026 meta-analysis", "MONO": "Vuckovic 2020",
        "NEUT": "Vuckovic 2020", "IBD": "Jostins 2012", "RA": "Okada 2014", "N_LDL": "Graham 2021", "TC": "Zoodsma 2025", "N_HEIGHT": "Yengo 2022", "BMD": "Morris 2019"}
OWNL = {t: ("DNASE" if t in ("N_HEIGHT", "BMD") else "CAGE") for t in B14}                    # 재추출 지도에서 τ*가 가장 큰 층위
M_TOT = sum(float(open(X/"GRCh38"/"baselineLD_v2.2"/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
def sd_of(*dirs_or_files):
    s = pd.concat([pd.read_csv(p if str(p).endswith(".csv") else p/"annot_sd_common.csv", index_col=0).sd_common for p in dirs_or_files]); assert not s.index.duplicated().any(); return s
SD_C = sd_of(C3); SD_D = sd_of(Y); SD_T3 = sd_of(Y, F3, I3); SD_PD = sd_of(K3); SD_M = sd_of(Y, E3/"measured_sd_common.csv"); SD_A = sd_of(Y, L3); SD_Q = sd_of(Q3); SD_O = sd_of(O3)
strip = lambda n: re.sub(r"L2_\d+$", "", n)
def jse(v): v = np.asarray(v, float); B = len(v); return float(np.sqrt((B - 1) / B * np.sum((v - v.mean()) ** 2)))
def jcov(Z): B = Z.shape[0]; z = Z - Z.mean(0); return (B - 1) / B * z.T @ z
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n); m = 1.0
    for r in range(n - 1, -1, -1): m = min(m, p[o[r]] * n / (r + 1)); q[o[r]] = m
    return q
def load(path, SD):
    """주석 이름 → (τ*, delete 값). 정렬 블록 npz."""
    z = np.load(path); h2 = float(z["h2_log"]); out = {}
    for k, n in enumerate([strip(n) for n in z["names"]]):
        if n in SD.index: f = SD[n] * M_TOT / h2; out[n] = (float(z["coef"][k] * f), z["part_delete"][:, k] * f)
    return out
def nk(t): return t[2:] if t.startswith("N_") else t                                          # 3n의 실행 파일 이름(N_ 없음)
def runs_main(t): return R3/"runs" if t == "N_AD" else (N3/"runs" if t in NEW6 else J3/"runs")
def key_main(t): return "N_AD" if t == "N_AD" else nk(t)
def p_m1(t, g): return (W/"runs"/f"{t}__M1_{g}.npz") if t not in NEW6 else runs_main(t)/f"{key_main(t)}__M1_{g}.npz"
def p_t3r(t, s): return runs_main(t)/f"{key_main(t)}__T3_{s:02d}.npz"
def p_t3all(t, g):
    if t in NEW6: return W/"runs"/f"{t}__T3all_{g}.npz"
    if g == OWN[t]: return J3/"runs"/f"{t}__T3all.npz"
    a = J3/"runs"/f"{t}__T3all_{g}.npz"; return a if a.exists() else T3/"runs"/f"{t}__T3all_{g}.npz"
def layers3(d):
    """세 층위의 (τ*, delete) — 이름 끝의 층위로 고른다(층위마다 주석 하나여야 한다)."""
    out = {}
    for m in LAY:
        ks = [n for n in d if n.endswith("_" + m)]; assert len(ks) == 1, (m, list(d)); out[m] = d[ks[0]]
    return out
def avg(ds):
    return {m: (float(np.mean([d[m][0] for d in ds])), np.mean([d[m][1] for d in ds], axis=0)) for m in LAY}
def phi(Lr, a="DNASE"):
    c, d = Lr["CAGE"], Lr[a]; e = float(np.degrees(np.arctan2(d[0], c[0]))); dv = np.degrees(np.arctan2(d[1], c[1])); dv = e + (dv - e + 180) % 360 - 180
    return e, dv, bool(c[0] + d[0] > 0)
def summ(Lr):
    c1 = Lr["CAGE"][0] - Lr["DNASE"][0]; d1 = Lr["CAGE"][1] - Lr["DNASE"][1]; c2 = Lr["CAGE"][0] - Lr["H3K27ac"][0]; d2 = Lr["CAGE"][1] - Lr["H3K27ac"][1]; e, dv, ok = phi(Lr)
    r = {m.lower(): Lr[m][0] for m in LAY}; r.update({m.lower() + "_se": jse(Lr[m][1]) for m in LAY})
    r.update(c1=c1, c1_se=jse(d1), c2=c2, c2_se=jse(d2), phiD=e, phiD_se=jse(dv), phiD_valid=ok); return r, dict(c1=d1, c2=d2, phi=dv)
def Qmv(est, dels):
    """공동 jackknife 다변량 Q: est = 형질마다의 벡터(길이 p), dels = (200 × p) 배열들. 공통 평균 GLS."""
    K = len(est); p = len(est[0]); d = np.concatenate(est); S = jcov(np.hstack(dels)); Si = np.linalg.inv(S); Xm = np.tile(np.eye(p), (K, 1))
    mu = np.linalg.solve(Xm.T @ Si @ Xm, Xm.T @ Si @ d); r = d - Xm @ mu; Q = float(r @ Si @ r); df = p * (K - 1); nu = dels[0].shape[0] - 1; Fs = (nu - df + 1) / (nu * df) * Q
    return Q, df, float(stats.chi2.sf(Q, df)), float(stats.f.sf(Fs, df, nu - df + 1)), mu
def ivw(items):
    """items = [(추정값, delete 벡터)] → 역분산 가중 평균과 그 delete 벡터(공동 jackknife)."""
    w = np.array([1 / jse(d) ** 2 for _, d in items]); w = w / w.sum(); return float(sum(wi * e for wi, (e, _) in zip(w, items))), sum(wi * d for wi, (_, d) in zip(w, items))
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
def f2(x, n=2): return f"{x:+.{n}f}".replace("-", "−") if x == x else "–"

# ================= 1. 조직 묶음 지도 14 × 12 (M1; 정렬 블록), 전체 묶음 보정(M2) =================
rows = []
for t in B14:
    for g in LIN:
        a = load(p_m1(t, g), SD_C)[f"L5_{g}"]; b = load(W/"runs"/f"{t}__M2_{g}.npz", SD_C)[f"L5_{g}"]
        rows.append(dict(trait=t, group=g, tau=a[0], se=jse(a[1]), tau_m2=b[0], se_m2=jse(b[1])))
TG = pd.DataFrame(rows); TG["p1"] = stats.norm.sf(TG.tau / TG.se); TG["q"] = bh(TG.p1); TG["p1_m2"] = stats.norm.sf(TG.tau_m2 / TG.se_m2); TG["q_m2"] = bh(TG.p1_m2)
TG["rank"] = TG.groupby("trait").tau.rank(ascending=False).astype(int); TG["own"] = [OWN[t] == g for t, g in zip(TG.trait, TG.group)]; TG["enriched"] = (TG.p1 < .05 / 12) & (TG.tau > 0)
TG.to_csv(RES/"w04_tissue.csv", index=False)
P("== 1. 조직 묶음 지도 (14 × 12; 단측 P, FDR은 168칸)"); P(f"FDR < 0.05인 칸: {int((TG.q < .05).sum())}개 / 168; 전체 묶음 보정 모형(M2)에서: {int((TG.q_m2 < .05).sum())}개")
for t in B14:
    d = TG[TG.trait == t].set_index("group"); o = d.loc[OWN[t]]; top = d.tau.idxmax()
    P(f"  {EN[t]:24s} 범주 묶음 {OWN[t]:11s} τ* {o.tau:.2f} ± {o.se:.2f} (P {o.p1:.1e}, 순위 {int(o['rank'])}) | 1위 {top} {d.tau[top]:.2f} | 농축 묶음(P<0.05/12) {int(d.enriched.sum())}개 | M2 {o.tau_m2:.2f} ± {o.se_m2:.2f} (P {o.p1_m2:.1e})")
P(f"범주 묶음이 1위인 형질: {int(sum(TG[(TG.trait == t) & TG.own]['rank'].iloc[0] == 1 for t in B14))}개 / 14")

# ================= 2. 72칸 지도 (층위 하나씩; ldsc 기본 블록) =================
ANN = [f"L5_{g}_{m}" for g in LIN for m in MODS6]; Y4 = pd.read_csv(EXP.parent/"e03_single_layer_map"/"results"/"y04_B.csv"); rows = []
for t in B14:
    f = W/"runs72"/f"{t}.npz"
    if f.exists():
        z = np.load(f); assert list(z["ann"]) == ANN; h2 = np.round(z["tot"], 4); fac = np.array([SD_D[a] for a in ANN]) * M_TOT / h2; ts, se = z["coef"] * fac, z["coef_se"] * fac; src = "3w"
    else:
        r = Y4[(Y4.trait == t) & (Y4.model == "D1")].set_index("annot").loc[ANN]; ts, se = r.tau_star.to_numpy(), r.tau_star_se.to_numpy(); src = "3d"
    for a, x, s in zip(ANN, ts, se): g = [g for g in LIN if a.startswith(f"L5_{g}_")][0]; rows.append(dict(trait=t, annot=a, group=g, layer=a[len(f"L5_{g}_"):], tau=x, se=s, source=src))
M72 = pd.DataFrame(rows); M72["p1"] = stats.norm.sf(M72.tau / M72.se); M72["q"] = bh(M72.p1); M72.to_csv(RES/"w04_map72.csv", index=False)
def typ(K, T, L): return "none" if K == 0 else (("tissue- and layer-spanning" if L >= 3 else "tissue-spanning") if T >= 3 else ("layer-spanning" if L >= 3 else "focal"))
rows = []
for t in B14:
    d = M72[M72.trait == t]; rec = dict(trait=t)
    for lab, thr in (("", .05), ("_q10", .10), ("_q01", .01)):
        s = d[d.q < thr]; K = len(s); T = s.group.nunique(); L = K / T if T else float("nan"); rec.update({"K" + lab: K, "T" + lab: T, "L" + lab: L, "type" + lab: typ(K, T, L)})
    s = d[d.q < .05]; rec["where"] = "; ".join(f"{g}: {', '.join(s[s.group == g].layer)}" for g in LIN if (s.group == g).any()); top = d.loc[d.tau.idxmax()]; rec.update(top_annot=top.annot, top_tau=top.tau)
    rows.append(rec)
TY = pd.DataFrame(rows); TY.to_csv(RES/"w04_types.csv", index=False)
P(f"\n== 2. 72칸 지도 (FDR은 14 × 72 = 1,008칸): 유의 {int((M72.q < .05).sum())}칸")
for r in TY.itertuples(): P(f"  {EN[r.trait]:24s} K {r.K:2d} T {r.T:2d} L {r.L:.2f} {r.type:26s} (FDR 0.10: {r.type_q10}; 0.01: {r.type_q01}) | 가장 큰 τ*: {r.top_annot} {r.top_tau:.2f}")
P("  유형별 수:", TY.type.value_counts().to_dict(), "| RNA-seq가 유의한 칸:", int(((M72.layer == "RNA") & (M72.q < .05)).sum()))

# ================= 3. 층위 함께 (재추출 20회 평균; 정렬 블록) =================
LR = {}; rows = []; DEL = {}
for t in B14:
    ds = [layers3(load(p_t3r(t, s), SD_T3)) for s in range(1, 21)]; LR[t] = avg(ds); r, DEL[t] = summ(LR[t]); first = [max(LAY, key=lambda m: d[m][0]) for d in ds]
    r.update(trait=t, name=EN[t], gwas=GWAS[t], group=OWN[t], top=max(LAY, key=lambda m: LR[t][m][0])); r["n_first"] = first.count(r["top"])
    sec = sorted(LAY, key=lambda m: -LR[t][m][0])[1]; dd = LR[t][r["top"]][1] - LR[t][sec][1]; r["top_vs_second_p1"] = float(stats.norm.sf((LR[t][r["top"]][0] - LR[t][sec][0]) / jse(dd)))
    r["resolved"] = bool(r["n_first"] >= 15 and r["top_vs_second_p1"] < .05)
    a = summ(layers3(load(p_t3all(t, OWN[t]), SD_T3)))[0]; r.update(all_c1=a["c1"], all_c1_se=a["c1_se"], all_phiD=a["phiD"], all_phiD_se=a["phiD_se"]); rows.append(r)
LY = pd.DataFrame(rows).set_index("trait"); LY.to_csv(RES/"w04_layers.csv")
P("\n== 3. 층위 함께, 재추출 20회 평균 (범주 묶음)")
for t in B14:
    r = LY.loc[t]; P(f"  {EN[t]:24s} CAGE {r.cage:.2f} ± {r.cage_se:.2f}  DNase {r.dnase:.2f} ± {r.dnase_se:.2f}  H3K27ac {r.h3k27ac:.2f} ± {r.h3k27ac_se:.2f} | C−D {f2(r.c1)} ± {r.c1_se:.2f} | φD {r.phiD:.1f}° ± {r.phiD_se:.1f}° | 1위 {r.top} {int(r.n_first)}/20, 둘째 대비 P {r.top_vs_second_p1:.3f}{' *' if r.resolved else ''} | 전체 트랙 C−D {f2(r.all_c1)} ± {r.all_c1_se:.2f}, φD {r.all_phiD:.1f}°")
def Q_of(ts, what):
    if what == "contrast": return Qmv([np.array([LY.c1[t], LY.c2[t]]) for t in ts], [np.column_stack([DEL[t]["c1"], DEL[t]["c2"]]) for t in ts])
    return Qmv([np.array([LY.phiD[t]]) for t in ts], [DEL[t]["phi"].reshape(-1, 1) for t in ts])
NOSK = [t for t in B14 if CAT[t] != "skeletal"]; QR = {}
for lab, ts in (("14형질", B14), ("골격 둘을 뺀 12형질", NOSK)):
    for what in ("contrast", "phi"):
        Q, df, p, pf, mu = Q_of(ts, what); QR[(lab, what)] = (Q, df, p, pf, mu)
        P(f"  이질성 [{lab}; {'(C−D, C−H)' if what == 'contrast' else 'φD'}]: Q = {Q:.1f}, 자유도 {df}, P = {p:.1e} (F 근사 {pf:.1e})" + (f"; 공통 φD {mu[0]:.1f}°" if what == "phi" else ""))
def contrast(A, Bs, key):
    ea, da = ivw([(LY[key][t] if key != "phi" else LY.phiD[t], DEL[t][key if key != "phi" else "phi"]) for t in A]); eb, db = ivw([(LY[key][t] if key != "phi" else LY.phiD[t], DEL[t][key if key != "phi" else "phi"]) for t in Bs])
    d = ea - eb; se = jse(da - db); return ea, jse(da), eb, jse(db), d, se, float(2 * stats.norm.sf(abs(d / se)))
CT = []
for lab, A, Bs in (("골격(키, BMD) − 뇌·행동 다섯", ["N_HEIGHT", "BMD"], [t for t in B14 if CAT[t] == "brain"]), ("골격 − 나머지 12", ["N_HEIGHT", "BMD"], NOSK), ("키 − 뇌·행동 다섯", ["N_HEIGHT"], [t for t in B14 if CAT[t] == "brain"]),
                   ("BMD − 뇌·행동 다섯", ["BMD"], [t for t in B14 if CAT[t] == "brain"]), ("키 − AD", ["N_HEIGHT"], ["N_AD"])):
    for key in ("c1", "phi"):
        ea, sa, eb, sb, d, se, p = contrast(A, Bs, key); CT.append(dict(contrast=lab, stat=key, a=ea, a_se=sa, b=eb, b_se=sb, diff=d, diff_se=se, p2=p))
        P(f"  대비 [{lab}; {'C−D' if key == 'c1' else 'φD'}]: {ea:.2f} ± {sa:.2f} 대 {eb:.2f} ± {sb:.2f}, 차이 {f2(d)} ± {se:.2f}, 양측 P = {p:.1e}")
pd.DataFrame(CT).to_csv(RES/"w04_contrasts.csv", index=False)
P("  규칙(20회 중 15회 이상 1위, 둘째 층위 대비 단측 P < 0.05)을 채운 형질:", ", ".join(f"{EN[t]}({LY.top[t]})" for t in B14 if LY.resolved[t]), "| 못 채운 형질:", ", ".join(EN[t] for t in B14 if not LY.resolved[t]))

# ================= 4. 농축된 모든 묶음의 층위 (전체 트랙) =================
rows = []
for t in B14:
    for g in LIN:
        r, _ = summ(layers3(load(p_t3all(t, g), SD_T3))); m = TG[(TG.trait == t) & (TG.group == g)].iloc[0]; r.update(trait=t, group=g, own=OWN[t] == g, m1_tau=m.tau, m1_p1=m.p1, enriched=bool(m.enriched), m1_rank=int(m["rank"])); rows.append(r)
AG = pd.DataFrame(rows); AG.to_csv(RES/"w04_allgroups.csv", index=False)
P("\n== 4. 농축된 모든 묶음의 층위(전체 트랙): 형질마다 농축 묶음(P < 0.05/12)의 φD")
for t in B14:
    d = AG[(AG.trait == t) & AG.enriched].sort_values("m1_rank"); P(f"  {EN[t]:24s} " + "; ".join(f"{r.group}{'*' if r.own else ''} {r.phiD:.0f}° ± {r.phiD_se:.0f}°" + ("" if r.phiD_valid else " (정의 안 됨)") for r in d.itertuples()))

# ================= 5. 세포형 일치, 원위·근위 =================
rows = []
for t in B14:
    g = OWN[t]; f = (W if t in NEW6 else J3)/"runs"/f"{t}__M3_{g}.npz"; r, dm = summ(layers3(load(f, SD_T3))); r.update(trait=t, model="matched"); rows.append((r, dm))
    f = (W if t in NEW6 else K3)/"runs"/f"{t}__PD_{g}.npz"; d = load(f, SD_PD)
    for cls, lab in (("D", "distal"), ("P", "proximal")):
        r, dm = summ({m: d[f"{cls}_{g}_{m}"] for m in LAY}); r.update(trait=t, model=lab); rows.append((r, dm))
CTRL = pd.DataFrame([r for r, _ in rows]); CTRL.to_csv(RES/"w04_controls.csv", index=False); CD = {(r["trait"], r["model"]): dm for r, dm in rows}
P("\n== 5. 세포형 일치 층위와 원위(> 1 kb)·근위 SNP의 층위: C−D")
for t in B14:
    d = CTRL[CTRL.trait == t].set_index("model"); P(f"  {EN[t]:24s} 재추출 {f2(LY.c1[t])} | 세포형 일치 {f2(d.c1['matched'])} ± {d.c1_se['matched']:.2f} | 원위 {f2(d.c1['distal'])} ± {d.c1_se['distal']:.2f} | 근위 {f2(d.c1['proximal'])} ± {d.c1_se['proximal']:.2f}")
for lab in ("matched", "distal", "proximal"):
    d = CTRL[CTRL.model == lab].set_index("trait"); Q, df, p, pf, _ = Qmv([np.array([d.c1[t], d.c2[t]]) for t in B14], [np.column_stack([CD[(t, lab)]["c1"], CD[(t, lab)]["c2"]]) for t in B14])
    same = sum((d.c1[t] > 0) == (LY.c1[t] > 0) for t in B14); P(f"  [{lab}] Q = {Q:.1f} (자유도 {df}, P = {p:.1e}); C−D의 부호가 재추출 지도와 같은 형질 {same}/14")

# ================= 6. 실측 주석 조건화 =================
rows = []; BMD_ = {}
for t in B14:
    g = OWN[t]; m = OWNL[t]; jm = "J1D" if m == "DNASE" else "J1C"; d = load(W/"runs"/f"{t}__{jm}_{g}.npz", SD_M); a = d[f"L5_{g}_{m}"]
    alone = M72[(M72.trait == t) & (M72.annot == f"L5_{g}_{m}")].tau.iloc[0]; rec = dict(trait=t, group=g, layer=m, alone=alone, cond=a[0], cond_se=jse(a[1]), retained=a[0] / alone)
    for mm in LAY: x = d[f"M_{g}_{mm}"]; rec[f"meas_{mm}"] = x[0]; rec[f"meas_{mm}_se"] = jse(x[1])
    c = load(W/"runs"/f"{t}__J1C_{g}.npz", SD_M)[f"L5_{g}_CAGE"]; rec.update(cage_alone=M72[(M72.trait == t) & (M72.annot == f"L5_{g}_CAGE")].tau.iloc[0], cage_cond=c[0], cage_cond_se=jse(c[1]))
    bm = layers3({k: v for k, v in load(W/"runs"/f"{t}__BM_{g}.npz", SD_M).items() if k.startswith("M_")}); sr, sdm = summ(bm); rec.update(bm_c1=sr["c1"], bm_c1_se=sr["c1_se"], bm_c2=sr["c2"]); BMD_[t] = sdm; rows.append(rec)
ME = pd.DataFrame(rows); ME["p1"] = stats.norm.sf(ME.cond / ME.cond_se); ME["q"] = bh(ME.p1); ME["cage_p1"] = stats.norm.sf(ME.cage_cond / ME.cage_cond_se); ME["cage_q"] = bh(ME.cage_p1)
Qb = Qmv([np.array([r.bm_c1, r.bm_c2]) for r in ME.itertuples()], [np.column_stack([BMD_[t]["c1"], BMD_[t]["c2"]]) for t in ME.trait]); ME.to_csv(RES/"w04_measured.csv", index=False)
P("\n== 6. 실측 CAGE·DNase·H3K27ac 주석과 함께 (형질의 층위: 키·BMD는 DNase, 나머지는 CAGE)")
for r in ME.itertuples(): P(f"  {EN[r.trait]:24s} {r.layer:5s} 단독 {r.alone:.2f} → 조건부 {r.cond:.2f} ± {r.cond_se:.2f} ({r.retained * 100:.0f}%), 단측 P {r.p1:.1e}, q {r.q:.3f}{' *' if r.q < .05 else ''}")
P(f"  FDR < 0.05: {int((ME.q < .05).sum())}/14; 유지율 중앙값 {ME.retained.median() * 100:.0f}% (범위 {ME.retained.min() * 100:.0f}–{ME.retained.max() * 100:.0f}%) | 모든 형질을 예측 CAGE로: FDR < 0.05 {int((ME.cage_q < .05).sum())}/14")
P(f"  실측 주석만의 층위 이질성(14형질): Q = {Qb[0]:.1f} (자유도 {Qb[1]}, P = {Qb[2]:.2g}; F 근사 {Qb[3]:.2g})")

# ================= 7. 활성도 조건화(상위 5%)와 활성 요소 안의 효과 크기 =================
def p_ea(t): return R3/"runs"/"N_AD__EA.npz" if t == "N_AD" else (N3/"runs"/f"{nk(t)}__EA.npz" if t in NEW6 else L3/"runs"/f"{t}__EA.npz")
def p_hl(t): return R3/"runs"/"N_AD__HL10.npz" if t == "N_AD" else Q3/"runs"/f"{t}__HL10.npz"
rows = []
for t in B14:
    g = OWN[t]; m = OWNL[t]; d = load(p_ea(t), SD_A); e = d[f"L5_{g}_{m}"]; a = d[f"A5_{g}_{m}"]; e0 = layers3(load(p_t3all(t, g), SD_T3))[m]
    h = load(p_hl(t), SD_Q); hi, lo = h[f"H10_{g}_{m}"], h[f"Lo10_{g}_{m}"]
    rows.append(dict(trait=t, group=g, layer=m, effect_only=e0[0], effect=e[0], effect_se=jse(e[1]), retained=e[0] / e0[0], activity=a[0], activity_se=jse(a[1]),
                     high=hi[0], high_se=jse(hi[1]), low=lo[0], low_se=jse(lo[1]), delta=hi[0] - lo[0], delta_se=jse(hi[1] - lo[1])))
AC = pd.DataFrame(rows); AC["p1"] = stats.norm.sf(AC.effect / AC.effect_se); AC["q"] = bh(AC.p1); AC["delta_p1"] = stats.norm.sf(AC.delta / AC.delta_se); AC["delta_q"] = bh(AC.delta_p1); AC.to_csv(RES/"w04_activity.csv", index=False)
P("\n== 7. 활성도 상위 5%를 공변량으로 (효과 주석의 τ*), 활성도 상위 10% 안의 효과 큰 절반 − 작은 절반")
for r in AC.itertuples():
    P(f"  {EN[r.trait]:24s} 효과만 {r.effect_only:.2f} → 활성도와 함께 {r.effect:.2f} ± {r.effect_se:.2f} ({r.retained * 100:.0f}%), q {r.q:.3f}{' *' if r.q < .05 else ''}; 활성도 τ* {r.activity:.2f} ± {r.activity_se:.2f} | Δ {f2(r.delta)} ± {r.delta_se:.2f}, q {r.delta_q:.3f}{' *' if r.delta_q < .05 else ''}")
P(f"  효과 주석이 남는 형질(FDR < 0.05): {int((AC.q < .05).sum())}/14 — " + ", ".join(EN[t] for t in AC.trait[AC.q < .05]) + f" | 활성 요소 안 Δ > 0 (FDR < 0.05): {int((AC.delta_q < .05).sum())}/14 — " + ", ".join(EN[t] for t in AC.trait[AC.delta_q < .05]))

# ================= 8. 효과 크기 구간 (겹치지 않는 넷) =================
BINS = {"B1": "0–0.1%", "B2": "0.1–0.5%", "B3": "0.5–1%", "B4": "1–5%"}; rows = []
for t in B14:
    g = OWN[t]; f = W/"runs"/"N_AD__BIN_myeloid.npz" if t == "N_AD" else O3/"runs"/f"{t}__BIN.npz"; d = load(f, SD_O)
    for b, lab in BINS.items():
        c, dn = d[f"{b}_{g}_CAGE"], d[f"{b}_{g}_DNASE"]; rows.append(dict(trait=t, bin=lab, cage=c[0], cage_se=jse(c[1]), dnase=dn[0], dnase_se=jse(dn[1]), c1=c[0] - dn[0], c1_se=jse(c[1] - dn[1])))
BN = pd.DataFrame(rows); BN.to_csv(RES/"w04_bins.csv", index=False)
P("\n== 8. 효과 크기 구간별 C−D (12개 구간 주석을 함께)")
for t in B14: d = BN[BN.trait == t]; P(f"  {EN[t]:24s} " + " | ".join(f"{r.bin} {f2(r.c1)} ± {r.c1_se:.2f}" for r in d.itertuples()))
for b in BINS.values():
    d = BN[BN.bin == b].set_index("trait"); P(f"  [{b}] C−D의 부호가 재추출 지도와 같은 형질: {sum((d.c1[t] > 0) == (LY.c1[t] > 0) for t in B14)}/14")

# ================= 9. 영향이 큰 블록 (전체 트랙, 범주 묶음; delete 값으로 계산) =================
rows = []
for t in B14:
    Lr = layers3(load(p_t3all(t, OWN[t]), SD_T3)); est = Lr["CAGE"][0] - Lr["DNASE"][0]; dv = Lr["CAGE"][1] - Lr["DNASE"][1]; infl = est - dv; o = np.argsort(-infl * np.sign(est))
    rec = dict(trait=t, c1=est, c1_se=jse(dv), top1_share=float((infl[o[0]]) / est))
    for k in (5, 10): rec[f"c1_without_top{k}"] = float(est - infl[o[:k]].sum()); rec[f"top{k}_blocks"] = " ".join(str(int(b)) for b in o[:k])
    rows.append(rec)
BL = pd.DataFrame(rows); BL.to_csv(RES/"w04_blocks.csv", index=False)
P("\n== 9. C−D에 영향이 큰 블록(전체 트랙): 추정값과, 대비를 가장 크게 만드는 블록 5·10개의 기여를 뺀 값(일차 근사; 다시 적합한 값이 아니다)")
for r in BL.itertuples(): P(f"  {EN[r.trait]:24s} {f2(r.c1)} ± {r.c1_se:.2f} → 5개 빼면 {f2(r.c1_without_top5)}, 10개 빼면 {f2(r.c1_without_top10)} (가장 큰 블록 하나의 몫 {r.top1_share * 100:.0f}%)")
(RES/"w04_summary.txt").write_text("\n".join(log), encoding="utf-8")
