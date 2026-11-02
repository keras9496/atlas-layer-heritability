"""실험 3o o03: 연속·구간 기준의 층위 대비 (PROTOCOL_exp3o_2026-10-21, 2026-10-21 고정).
 τ* = 계수 × sd(공통 SNP) × M_5_50 / h²(로그 값), delete 값도 같은 배율(3j와 같다). Q·결합 대비는 3m m_common(= 3j j03)과 같다.
 기준: CONT(연속 백분위 3층위 동시), B1 0–0.1 %, B2 0.1–0.5 %, B3 0.5–1 %, B4 1–5 % (BIN 모형에서 12개 동시).
 O-H1 (주) CONT Q (지도 형질 11개, 자유도 20) P < 0.05 → "연속 주석에서도 형질 간 층위 이질성".
 O-H2 (주) 형질마다 c1 부호가 전체 트랙 5 % 모형(T3all)과 CONT·B1·B2·B3·B4 다섯 기준 모두 같으면 "방향 안정", 아니면 "기준 민감".
 O-H3 구간별 Q. O-H4 사전 대비 세 개와 검증 대비. O-H5 새 GWAS 다섯 형질(기술).
 기술: CONT·BIN 주석 LD 점수의 기준선 보정 VIF (계통별). BIN에서 VIF > 5인 계통의 구간 대비는 해석하지 않는다.
출력: results/o03_tests.txt, o03_contrasts_long.csv (포함 관계 역치 3m M2와 나란히), o03_stability.csv, o03_vif.csv
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
from m_common import jse, jcov, Qtest, comb, CON, LAY, MAP, VAL, ALL, KO, M_TOT, OWN
X = Path(f"{DATA}/exp3"); O = Path(f"{DATA}/exp3o"); OR = O/"runs"; NR = Path(f"{DATA}/exp3n/runs")
NEW = ["HEIGHT", "BMI", "SCZ", "BIP", "LDL"]; CRIT = ["CONT", "B1", "B2", "B3", "B4"]
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
def fp(p): return f"{p:.2g}"
SD = pd.read_csv(O/"annot_sd_common.csv", index_col=0).sd_common
def load(key, model):
    """기준 → 층위 → (τ*, delete 값)"""
    z = np.load(OR/f"{key}__{model}.npz"); h2 = float(z["h2_log"]); res = {}
    for k, n in enumerate(re.sub(r"L2_\d+$", "", x) for x in z["names"]):
        if n not in SD.index: continue
        crit = "CONT" if n.startswith("C_") else n.split("_")[0]; m = [x for x in LAY if n.endswith("_" + x)][0]; f = SD[n] * M_TOT / h2
        res.setdefault(crit, {})[m] = (z["coef"][k] * f, z["part_delete"][:, k] * f)
    assert sorted(res) == (["CONT"] if model == "CONT" else ["B1", "B2", "B3", "B4"]) and all(len(v) == 3 for v in res.values()), (key, model, list(res))
    return res
def prof_of(Lc): return (np.array([Lc["CAGE"][0] - Lc["DNASE"][0], Lc["CAGE"][0] - Lc["H3K27ac"][0]]), np.column_stack([Lc["CAGE"][1] - Lc["DNASE"][1], Lc["CAGE"][1] - Lc["H3K27ac"][1]]))
def profiles(keys):
    out = {c: {} for c in CRIT}
    for t, key in keys.items():
        c0 = load(key, "CONT"); b0 = load(key, "BIN"); out["CONT"][t] = prof_of(c0["CONT"])
        for b in CRIT[1:]: out[b][t] = prof_of(b0[b])
    return out
PR = profiles({t: t for t in ALL})
REF = {t: MC.prof(t, ["T3all"]) for t in ALL}   # 전체 트랙 5 % (3j)
# ---------------- O-H1 ----------------
P("=== O-H1 (주): 연속 백분위 주석의 형질 간 층위 이질성 (지도 형질 11개) ===")
Q = Qtest(PR["CONT"], MAP); P(f"  CONT Q {Q[0]:.1f} (자유도 {Q[1]}, P {fp(Q[2])}, F 근사 {fp(Q[3])}) → {'연속 주석에서도 형질 간 층위 이질성' if Q[2] < .05 else '미달'}")
Qr = Qtest(REF, MAP); P(f"  (참고: 전체 트랙 5 % 모형 Q {Qr[0]:.1f}, P {fp(Qr[2])})")
P("\n=== O-H3: 구간별 Q (BIN 모형, 11형질) ===")
for b in CRIT[1:]: q = Qtest(PR[b], MAP); P(f"  {b}: Q {q[0]:.1f}, P {fp(q[2])} (F 근사 {fp(q[3])})")
# ---------------- O-H2 ----------------
P("\n=== O-H2 (주): 형질별 c1 (CAGE − DNase) 부호 안정성 — 기준 T3all 부호와 CONT·B1–B4 비교 ===")
rows, longr = [], []
M2 = pd.read_csv(EXP.parent/"e12_diagnostics"/"results"/"m02_threshold_contrasts.csv")
def add_long(t, gw, crit, pr):
    c, d = pr; longr.append(dict(trait=t, gwas=gw, criterion=crit, c1=c[0], c1_se=jse(d[:, 0]), c2=c[1], c2_se=jse(d[:, 1])))
for t in ALL:
    s0 = np.sign(REF[t][0][0]); same = 0; excl = 0; txt = []
    for c in CRIT:
        e, d = PR[c][t]; se = jse(d[:, 0]); same += int(np.sign(e[0]) == s0); excl += int(abs(e[0]) > 1.96 * se); txt.append(f"{c} {e[0]:+.2f}±{se:.2f}")
        add_long(t, "original", c, PR[c][t])
    add_long(t, "original", "T3all_5%", REF[t])
    for _, r in M2[M2.trait == t].iterrows():
        if r.threshold != "5%": longr.append(dict(trait=t, gwas="original", criterion=f"nested_{r.threshold}", c1=r.c1, c1_se=r.c1_se, c2=r.c2, c2_se=r.c2_se))
    lab = "방향 안정" if same == len(CRIT) else "기준 민감"
    rows.append(dict(trait=t, gwas="original", ref_c1=REF[t][0][0], n_same_sign=same, n_ci_excl0=excl, label=lab))
    P(f"  {KO[t]:7s} T3all {REF[t][0][0]:+.2f} | " + " | ".join(txt) + f" → 같은 부호 {same}/5, 0 벗어남 {excl}/5 → {lab}")
P("  방향 안정: " + ", ".join(KO[r["trait"]] for r in rows if r["label"] == "방향 안정"))
P("  기준 민감: " + ", ".join(KO[r["trait"]] for r in rows if r["label"] == "기준 민감"))
# ---------------- O-H4 ----------------
P("\n=== O-H4: 사전 대비와 검증 대비 (기준별) ===")
VALC = ("Validation: NEUT, TC vs BMD (c1, one-sided)", [(1, ["NEUT", "TC"]), (-1, ["BMD"])], 0)
cr = []
for name, groups, k in CON + [VALC]:
    parts = []
    for c in ["T3all_5%"] + CRIT:
        pr = REF if c == "T3all_5%" else PR[c]; est, se, p2, _ = comb(pr, groups, k); p = float(stats.norm.sf(est / se)) if name.startswith("Validation") else p2
        parts.append(f"{c} {est:+.2f}±{se:.2f} (P {fp(p)})"); cr.append(dict(contrast=name, criterion=c, est=est, se=se, p=p))
    P(f"  {name}:\n     " + "\n     ".join(parts))
# ---------------- O-H5 새 GWAS ----------------
P("\n=== O-H5: 새 GWAS 다섯 형질 (기술; 기준은 새 GWAS의 T3all, 3n) ===")
PRN = profiles({t: f"N_{t}" for t in NEW})
for t in NEW:
    rn = MC.prof(t, ["T3all"], NR); s0 = np.sign(rn[0][0]); same = 0; excl = 0; txt = []
    for c in CRIT:
        e, d = PRN[c][t]; se = jse(d[:, 0]); same += int(np.sign(e[0]) == s0); excl += int(abs(e[0]) > 1.96 * se); txt.append(f"{c} {e[0]:+.2f}±{se:.2f}"); add_long(t, "new", c, PRN[c][t])
    add_long(t, "new", "T3all_5%", rn)
    lab = "방향 안정" if same == len(CRIT) else "기준 민감"; rows.append(dict(trait=t, gwas="new", ref_c1=rn[0][0], n_same_sign=same, n_ci_excl0=excl, label=lab))
    P(f"  {KO[t]:7s} T3all {rn[0][0]:+.2f} | " + " | ".join(txt) + f" → 같은 부호 {same}/5, 0 벗어남 {excl}/5 → {lab}")
UPD = {c: {t: (PRN[c][t] if t in NEW else PR[c][t]) for t in MAP} for c in CRIT}
q = Qtest(UPD["CONT"], MAP); P(f"  갱신 지도(새 GWAS 다섯 개) CONT Q {q[0]:.1f}, P {fp(q[2])}")
for b in CRIT[1:]: q = Qtest(UPD[b], MAP); P(f"  갱신 지도 {b} Q {q[0]:.1f}, P {fp(q[2])}")
# ---------------- 기술: VIF ----------------
P("\n=== 기술: 주석 LD 점수의 기준선 보정 VIF (회귀 SNP) ===")
BL = X/"GRCh38"/"baselineLD_v2.2"; Wt = X/"GRCh38"/"weights"
wsnp = set(pd.concat([pd.read_csv(Wt/f"weights.hm3_noMHC.{c}.l2.ldscore.gz", sep="\t", usecols=["SNP"]) for c in range(1, 23)]).SNP)
BLD = pd.concat([pd.read_csv(BL/f"baselineLD.{c}.l2.ldscore.gz", sep="\t") for c in range(1, 23)]).drop_duplicates("SNP"); BLD = BLD[BLD.SNP.isin(wsnp)].reset_index(drop=True)
Xb = np.column_stack([np.ones(len(BLD)), BLD[[c for c in BLD.columns if c not in ("CHR", "SNP", "BP")]].to_numpy(np.float64)])
def ldmat(grp):
    Lx = pd.concat([pd.read_csv(O/"annot"/grp/f"{grp}.{c}.l2.ldscore.gz", sep="\t") for c in range(1, 23)]).drop_duplicates("SNP").set_index("SNP")
    return Lx.loc[BLD.SNP].drop(columns=["CHR", "BP"]).to_numpy(np.float64), [re.sub(r"L2$", "", x) for x in Lx.columns if x not in ("CHR", "BP")]
vr = []
for g in ["myeloid", "CNS_tissue", "lymphoid", "liver", "stromal"]:
    for grp in [f"C3_{g}", f"B12_{g}"]:
        Lm, names = ldmat(grp); Rz = Lm - Xb @ np.linalg.lstsq(Xb, Lm, rcond=None)[0]; Cm = np.corrcoef(Rz, rowvar=False); vif = np.diag(np.linalg.inv(Cm))
        vr += [dict(group=grp, annot=n, vif=v) for n, v in zip(names, vif)]
        P(f"  {grp:15s} VIF {vif.min():.2f}–{vif.max():.2f}" + (" → 해석 제외" if grp.startswith("B12") and vif.max() > 5 else ""))
pd.DataFrame(vr).to_csv(EXP/"results"/"o03_vif.csv", index=False)
pd.DataFrame(longr).to_csv(EXP/"results"/"o03_contrasts_long.csv", index=False); pd.DataFrame(rows).to_csv(EXP/"results"/"o03_stability.csv", index=False)
pd.DataFrame(cr).to_csv(EXP/"results"/"o03_pairwise.csv", index=False)
(EXP/"results"/"o03_tests.txt").write_text("\n".join(log), encoding="utf-8")
