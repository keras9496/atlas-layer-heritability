"""실험 3l l04: 활성도 조건화 검정 (PROTOCOL_exp3l_2026-10-17, 2026-10-18 고정).
 τ* = 계수 × sd(공통 SNP) × M_5_50 / h²(로그의 반올림 값), delete 값도 같은 배율 (3j와 같다). 정렬 블록 200개, 공동 jackknife.
 주석 종류: 효과 L5_ (sd: 3d), 활성도 A5_ (sd: 3l l02), 실측 M_ (sd: 3e). pretest에서는 3i 대역 M5_ 를 활성도로 읽는다.
 자기 층위: 키 DNase, 나머지 지도 형질 10개 CAGE.
 A-H1 (주) EA에서 자기 층위 효과 τ* > 0, 단측, 지도 형질 11개 BH → FDR < 0.05가 6개 이상 "활성도 조건부로 유지", 3–5 "부분", ≤ 2 "대부분 활성도로 설명".
 A-H2 EA 효과 대비 Q(11), 사전 대비 3개(양측 α 0.05/3), 검증 대비(단측). A-H3 유지 비율 EA/E, 활성도 τ*. A-H4 A 모형 활성도 대비 Q, 1위 층위. A-H5 EAM 자기 층위 효과 τ* BH.
 기술: 효과·활성도 6개 LD 점수(기준선 97개를 뺀 잔차)의 VIF와 같은 층위 편상관. 자기 층위 효과의 VIF > 5인 형질은 A-H1에서 해석하지 않는다(판정 계산에는 넣는다).
 사용: python l04_tests.py pretest | dry | main
출력: results/l04_<단계>_tests.txt, l04_<단계>_summary.csv, l04_main_vif.csv
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]
X = Path(f"{DATA}/exp3"); Y = Path(f"{DATA}/exp3d"); I = Path(f"{DATA}/exp3i"); W = Path(f"{DATA}/exp3e")
J = Path(f"{DATA}/exp3j"); L = Path(f"{DATA}/exp3l"); RUN = L/"test" if STAGE in ("pretest", "dry") else L/"runs"
log, SUM = [], []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LAY = ["CAGE", "DNASE", "H3K27ac"]
SDp = [Y/"annot_sd_common.csv", W/"measured_sd_common.csv"] + ([I/"annot_sd_common.csv"] if STAGE == "pretest" else [L/"annot_sd_common.csv"])
SD = pd.concat([pd.read_csv(p, index_col=0).sd_common for p in SDp]); SD = SD[~SD.index.duplicated()]
B0 = X/"GRCh38"/"baselineLD_v2.2"; M_TOT = sum(float(open(B0/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
OWN = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
       "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal", "BMD": "stromal", "NEUT": "myeloid", "TC": "liver"}
MAP = ["AD", "MONO", "EA", "SCZ", "BIP", "NEURO", "BMI", "IBD", "RA", "LDL", "HEIGHT"]; VAL = ["BMD", "NEUT", "TC"]; ALL = MAP + VAL
OWNL = {t: ("DNASE" if t == "HEIGHT" else "CAGE") for t in MAP}
KO = {"AD": "AD", "MONO": "단핵구 수", "EA": "교육 연수", "SCZ": "조현병", "BIP": "양극성장애", "NEURO": "신경증", "BMI": "BMI", "IBD": "IBD", "RA": "RA", "LDL": "LDL", "HEIGHT": "키",
      "BMD": "골밀도", "NEUT": "호중구 수", "TC": "총콜레스테롤", "FAKE1": "가짜", "FAKEpre": "가짜(대역)"}
KIND = {"L5_": "E", "A5_": "A", "M_": "M"} if STAGE != "pretest" else {"L5_": "E", "M5_": "A", "M_": "M"}
strip = lambda n: re.sub(r"L2_\d+$", "", n)
def load(t, model, runs):
    """종류(E/A/M) → 층위 → (τ*, delete 값)"""
    z = np.load(runs/f"{t}__{model}.npz"); h2 = float(z["h2_log"]); res = {}
    for k, n in enumerate(strip(x) for x in z["names"]):
        if n not in SD.index: continue
        kind = [v for p, v in KIND.items() if n.startswith(p)][0]; m = [x for x in LAY if n.endswith("_" + x)][0]; f = SD[n] * M_TOT / h2
        res.setdefault(kind, {})[m] = (z["coef"][k] * f, z["part_delete"][:, k] * f)
    return res
def jcov(Z): B = Z.shape[0]; z = Z - Z.mean(0); return (B - 1) / B * z.T @ z
def jse(v): return float(np.sqrt(jcov(np.asarray(v).reshape(-1, 1))[0, 0]))
def prof(L3):
    return (np.array([L3["CAGE"][0] - L3["DNASE"][0], L3["CAGE"][0] - L3["H3K27ac"][0]]), np.column_stack([L3["CAGE"][1] - L3["DNASE"][1], L3["CAGE"][1] - L3["H3K27ac"][1]]))
def Qtest(PR, ts):
    d = np.concatenate([PR[t][0] for t in ts]); S = jcov(np.hstack([PR[t][1] for t in ts])); K = len(ts)
    Si = np.linalg.inv(S); Xm = np.tile(np.eye(2), (K, 1)); mu = np.linalg.solve(Xm.T @ Si @ Xm, Xm.T @ Si @ d); r = d - Xm @ mu
    Q = float(r @ Si @ r); df = 2 * (K - 1); nu = PR[ts[0]][1].shape[0] - 1; Fs = (nu - df + 1) / (nu * df) * Q
    return Q, df, float(stats.chi2.sf(Q, df)), float(stats.f.sf(Fs, df, nu - df + 1))
def comb(PR, groups, k, side):
    a = {}
    for sgn, ts in groups:
        w = np.array([1 / jse(PR[t][1][:, k]) ** 2 for t in ts])
        for wi, t in zip(w, ts): a[t] = a.get(t, 0) + sgn * wi / w.sum()
    est = float(sum(ai * PR[t][0][k] for t, ai in a.items())); se = jse(sum(ai * PR[t][1][:, k] for t, ai in a.items()))
    return est, se, float(stats.norm.sf(est / se) if side == 1 else 2 * stats.norm.sf(abs(est / se)))
def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p); q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]; return np.minimum(q, 1)
if STAGE in ("pretest", "dry"):
    tag = "FAKEpre" if STAGE == "pretest" else "FAKE1"
    for model in (["EA", "A", "EAM"] if STAGE == "pretest" else ["EA"]):
        R = load(tag, model, RUN); P(f"[{STAGE}] {tag} {model}: " + " | ".join(f"{k} " + ", ".join(f"{m} {v[0]:+.2f} ± {jse(v[1]):.2f}" for m, v in R[k].items()) for k in sorted(R)))
        assert sorted(R) == {"EA": ["A", "E"], "A": ["A"], "EAM": ["A", "E", "M"]}[model] and all(len(R[k]) == 3 for k in R), "주석 종류·층위 수가 맞지 않다"
    P("형식 확인 통과 (가짜 형질 — 결과 아님)"); (EXP/"results"/f"l04_{STAGE}_tests.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit()
# ---------------- 본 검정 ----------------
E = {t: {"E": load(t, "T3all", J/"runs")["E"]} for t in ALL}           # 기준 E = 3j T3all (효과 3층위만)
EA = {t: load(t, "EA", RUN) for t in ALL}; A = {t: load(t, "A", RUN) for t in ALL}; EAM = {t: load(t, "EAM", RUN) for t in ALL}
assert abs(E["AD"]["E"]["CAGE"][0] - 2.91) < .006, E["AD"]["E"]["CAGE"][0]   # 기준 E가 3m M3의 AD T3all CAGE τ*(+2.91)와 같은 계산인지
# --- 기술: 효과·활성도 LD 점수 VIF (계통별, 기준선 잔차) ---
P("=== 기술: 효과·활성도 주석의 LD 점수 (회귀 SNP, 기준선 97개를 뺀 잔차) ===")
BL = X/"GRCh38"/"baselineLD_v2.2"; Wt = X/"GRCh38"/"weights"
wsnp = set(pd.concat([pd.read_csv(Wt/f"weights.hm3_noMHC.{c}.l2.ldscore.gz", sep="\t", usecols=["SNP"]) for c in range(1, 23)]).SNP)
BLD = pd.concat([pd.read_csv(BL/f"baselineLD.{c}.l2.ldscore.gz", sep="\t") for c in range(1, 23)]).drop_duplicates("SNP"); BLD = BLD[BLD.SNP.isin(wsnp)].reset_index(drop=True)
Xb1 = np.column_stack([np.ones(len(BLD)), BLD[[c for c in BLD.columns if c not in ("CHR", "SNP", "BP")]].to_numpy(np.float64)])
def ldmat(d, stem, names):
    Lx = pd.concat([pd.read_csv(d/f"{stem}.{c}.l2.ldscore.gz", sep="\t") for c in range(1, 23)]).drop_duplicates("SNP").set_index("SNP")
    return Lx.loc[BLD.SNP, [f"{n}L2" for n in names]].to_numpy(np.float64)
VIF = {}; vrows = []
for g in sorted(set(OWN.values())):
    Le = np.hstack([ldmat(Y/"annot"/f"L5_{g}_{m}", f"L5_{g}_{m}", [f"L5_{g}_{m}"]) for m in LAY]); La = ldmat(L/"annot"/f"A5_{g}", f"A5_{g}", [f"A5_{g}_{m}" for m in LAY])
    Lm = np.hstack([Le, La]); Rz = Lm - Xb1 @ np.linalg.lstsq(Xb1, Lm, rcond=None)[0]; pr = np.corrcoef(Rz, rowvar=False); v = np.diag(np.linalg.inv(pr))
    VIF[g] = dict(zip([f"E_{m}" for m in LAY] + [f"A_{m}" for m in LAY], v))
    vrows.append(dict(lineage=g, **VIF[g], **{f"pr_EA_{m}": pr[i, i + 3] for i, m in enumerate(LAY)}))
    P(f"  {g:11s} VIF " + " ".join(f"{k} {x:.2f}" for k, x in VIF[g].items()) + " | 같은 층위 효과–활성도 편상관 " + " ".join(f"{m} {pr[i, i + 3]:.2f}" for i, m in enumerate(LAY)))
pd.DataFrame(vrows).to_csv(EXP/"results"/"l04_main_vif.csv", index=False)
# --- A-H1 (주) ---
P("\n=== A-H1 (주): EA에서 자기 층위 효과 τ* > 0 (단측, 지도 형질 11개 BH) ===")
p1 = []
for t in MAP:
    m = OWNL[t]; e, d = EA[t]["E"][m]; s = jse(d); p1.append(stats.norm.sf(e / s))
q1 = bh(p1); n1 = int((q1 < .05).sum())
for t, p, q in zip(MAP, p1, q1):
    m = OWNL[t]; e, d = EA[t]["E"][m]; e0 = E[t]["E"][m][0]; vf = VIF[OWN[t]][f"E_{m}"]
    P(f"  {KO[t]:6s} {m:7s} E {e0:+.2f} → EA {e:+.2f} ± {jse(d):.2f} (P {p:.2g}, q {q:.2g}) | 유지 비율 {e / e0:.0%} | VIF {vf:.2f}{' — 해석하지 않음' if vf > 5 else ''}")
    SUM.append(dict(test="A-H1", trait=t, layer=m, tau_E=e0, tau_EA=e, se_EA=jse(d), p=p, q=q, vif=vf))
v1 = "활성도 조건부로 유지" if n1 >= 6 else ("부분" if n1 >= 3 else "대부분 활성도로 설명")
P(f"  → FDR < 0.05: {n1}/11 → A-H1 판정: {v1}")
# --- A-H2 ---
P("\n=== A-H2: EA 효과 대비 ===")
PRe = {t: prof(EA[t]["E"]) for t in ALL}; Q, df, p, pf = Qtest(PRe, MAP)
P(f"  다변량 Q (지도 형질 11개) {Q:.1f}, 자유도 {df}, P {p:.2g} (F 근사 {pf:.2g})"); SUM.append(dict(test="A-H2 Q", tau_EA=Q, p=p, q=pf))
CON = [("AD 대 키 (c1)", [(1, ["AD"]), (-1, ["HEIGHT"])], 0, 2), ("뇌 형질 4개 대 키 (c1)", [(1, ["EA", "SCZ", "BIP", "NEURO"]), (-1, ["HEIGHT"])], 0, 2),
       ("IBD·RA 대 AD·단핵구 수 (c2)", [(1, ["IBD", "RA"]), (-1, ["AD", "MONO"])], 1, 2), ("검증 대비: 호중구 수·총콜레스테롤 대 골밀도 (c1)", [(1, ["NEUT", "TC"]), (-1, ["BMD"])], 0, 1)]
for lab, gr, k, side in CON:
    e, s, pp = comb(PRe, gr, k, side); P(f"  {lab}: {e:+.2f} ± {s:.2f}, P {pp:.2g} [{'단측' if side == 1 else '양측, α 0.0167'}]"); SUM.append(dict(test=f"A-H2 {lab}", tau_EA=e, se_EA=s, p=pp))
P("  형질별 EA 효과 τ* C/D/H 와 1위: " + "; ".join(f"{KO[t]} " + "/".join(f"{EA[t]['E'][m][0]:+.2f}" for m in LAY) + f" {max(LAY, key=lambda m: EA[t]['E'][m][0])}" for t in ALL))
# --- A-H3 ---
P("\n=== A-H3: 유지 비율과 활성도 τ* (EA) ===")
for t in ALL:
    P(f"  {KO[t]:6s} 효과 E → EA: " + ", ".join(f"{m} {E[t]['E'][m][0]:+.2f} → {EA[t]['E'][m][0]:+.2f}" for m in LAY) + " | 활성도 " + ", ".join(f"{m} {EA[t]['A'][m][0]:+.2f} ± {jse(EA[t]['A'][m][1]):.2f}" for m in LAY))
rr = [EA[t]["E"][OWNL[t]][0] / E[t]["E"][OWNL[t]][0] for t in MAP]; P(f"  자기 층위 유지 비율 중앙값 {np.median(rr):.0%} (범위 {min(rr):.0%}–{max(rr):.0%})")
# --- A-H4 ---
P("\n=== A-H4: 활성도만의 층위 지도 (A 모형) ===")
PRa = {t: prof(A[t]["A"]) for t in ALL}; Q, df, p, pf = Qtest(PRa, MAP)
P(f"  활성도 대비 Q (지도 형질 11개) {Q:.1f}, P {p:.2g} (F 근사 {pf:.2g}) | 1위: " + ", ".join(f"{KO[t]} {max(LAY, key=lambda m: A[t]['A'][m][0])}" for t in ALL))
SUM.append(dict(test="A-H4 Q", tau_EA=Q, p=p, q=pf))
# --- A-H5 ---
P("\n=== A-H5: EAM(실측 주석까지)에서 자기 층위 효과 τ* ===")
p5 = [stats.norm.sf(EAM[t]["E"][OWNL[t]][0] / jse(EAM[t]["E"][OWNL[t]][1])) for t in MAP]; q5 = bh(p5)
for t, p, q in zip(MAP, p5, q5): e, d = EAM[t]["E"][OWNL[t]]; P(f"  {KO[t]:6s} {OWNL[t]:7s} {e:+.2f} ± {jse(d):.2f} (P {p:.2g}, q {q:.2g})"); SUM.append(dict(test="A-H5", trait=t, layer=OWNL[t], tau_EA=e, se_EA=jse(d), p=p, q=q))
P(f"  → FDR < 0.05: {int((q5 < .05).sum())}/11")
pd.DataFrame(SUM).to_csv(EXP/"results"/"l04_main_summary.csv", index=False)
(EXP/"results"/"l04_main_tests.txt").write_text("\n".join(log), encoding="utf-8")
