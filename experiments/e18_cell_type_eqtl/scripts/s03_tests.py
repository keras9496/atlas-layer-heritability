"""실험 3s s03: 세포형 eQTL MaxCPP 조건부 검정 (PROTOCOL_exp3s_2026-10-25, 2026-10-25 고정).
 τ* = 계수 × sd(공통 SNP) × M_5_50 / h²(로그 값), 공동 jackknife(정렬 블록 200개) — 3l·3p와 같은 식(3p p_common의 jcov·jse·prof·Qtest·comb·bh).
 S-H1 (주) AD + 뇌 조직 형질 5개(EA, SCZ, BIP, NEURO, BMI): EQ에서 자기 층위(CAGE) 효과 τ* > 0, 단측, 6개 BH.
        FDR < 0.05 형질 수 5–6 → "세포형 eQTL을 넘어 유지", 3–4 → "부분", 0–2 → "대부분 세포형 eQTL로 설명". 유지 비율(EQ/E) 중앙값을 함께 적는다.
 S-H2 (주, AD) EAQ에서 AD 골수계 CAGE 활성도 τ* > 0 (단측 α 0.05).
 S-H3 지도 11 + 검증 3: EQ·EAQ 자기 층위 효과 τ*와 유지 비율, EQ 층위 지도 Q(지도 11개)·사전 대비 3개·검증 대비. S-H4 Q5 각 주석의 τ*. S-H5 더 큰 GWAS 6개.
 기술: 계통별 기준선 보정 LD 점수 VIF(효과 3 + 활성도 3 + Q5).
 사용: python s03_tests.py dry | main
출력: results/s03_<단계>_tests.txt, s03_summary.csv, s03_values.csv, s03_vif.csv
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"e15_activity_constructions"/"scripts"))
from p_common import X, Y, L, J, N3, LAY, M_TOT, OWN, MAP, VAL, ALL, NEW, OWNL, KO, CON, jcov, jse, prof, Qtest, comb, bh
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; S3 = Path(f"{DATA}/exp3s"); R3 = Path(f"{DATA}/exp3r")
RUN = S3/"test" if STAGE == "dry" else S3/"runs"
QN = ["Q_MIC", "Q_NEU", "Q_GLIA", "Q_BULK", "Q_MONO"]
SD = pd.concat([pd.read_csv(p, index_col=0).sd_common for p in [Y/"annot_sd_common.csv", L/"annot_sd_common.csv", S3/"annot_sd_common.csv"]]); SD = SD[~SD.index.duplicated()]
OWN = dict(OWN, FAKE1="myeloid", **{f"N_{t}": OWN[t] for t in NEW + ["AD"]}); OWNL = dict(OWNL, **{f"N_{t}": OWNL[t] for t in NEW + ["AD"]}); KO = dict(KO, **{f"N_{t}": f"{KO[t]}(새)" for t in NEW + ["AD"]})
log, SUM, VALS = [], [], []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
def fp(p): return f"{p:.2g}"
strip = lambda n: re.sub(r"L2_\d+$", "", n)
def load(path):
    """종류(E/A/Q) → 층위 또는 Q 이름 → (τ*, delete 값)"""
    z = np.load(path); h2 = float(z["h2_log"]); res = {}
    for k, n in enumerate(strip(x) for x in z["names"]):
        if n not in SD.index: continue
        if n.startswith("Q_"): kind, m = "Q", n
        else: kind = "E" if n.startswith("L5_") else "A"; m = [x for x in LAY if n.endswith("_" + x)][0]
        f = SD[n] * M_TOT / h2; assert m not in res.get(kind, {}); res.setdefault(kind, {})[m] = (z["coef"][k] * f, z["part_delete"][:, k] * f)
    return res
def ref_path(t, model):
    """비교 대상: E = T3all, EA = EA (원래 3j·3l, 새 GWAS 3n, 새 AD 3r)"""
    if t == "N_AD": return R3/"runs"/f"N_AD__{ {'E': 'T3all', 'EA': 'EA'}[model] }.npz"
    if t.startswith("N_"): return N3/"runs"/f"{t[2:]}__{ {'E': 'T3all', 'EA': 'EA'}[model] }.npz"
    return {"E": J/"runs"/f"{t}__T3all.npz", "EA": L/"runs"/f"{t}__EA.npz"}[model]
if STAGE == "dry":
    for model in ["EQ", "EAQ"]:
        R = load(RUN/f"FAKE1__{model}.npz"); P(f"[dry] FAKE1 {model}: " + " | ".join(f"{k}: " + ", ".join(R[k]) for k in R))
        assert sorted(R["Q"]) == sorted(QN) and len(R["E"]) == 3 and (model == "EQ" or len(R["A"]) == 3)
    P("형식 확인 통과 (가짜 형질 — 결과 아님)"); (EXP/"results"/"s03_dry_tests.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit()
TR = ALL + [f"N_{t}" for t in NEW] + ["N_AD"]
RES = {(t, m): load(RUN/f"{t}__{m}.npz") for t in TR for m in ["EQ", "EAQ"]}
REF = {(t, m): load(ref_path(t, m)) for t in TR for m in ["E", "EA"]}
for (t, m), R in list(RES.items()) + list(REF.items()):
    for kind in R:
        for x, (tau, d) in R[kind].items(): VALS.append(dict(trait=t, model=m, kind=kind, layer=x, tau=float(tau), se=jse(d)))
pd.DataFrame(VALS).to_csv(EXP/"results"/"s03_values.csv", index=False)
# ---------------- 기술: VIF ----------------
P("=== 기술: 기준선 보정 LD 점수 VIF (효과 3 + 활성도 3 + Q5; 회귀 SNP, 기준선 97개를 뺀 잔차) ===")
BL = X/"GRCh38"/"baselineLD_v2.2"; Wt = X/"GRCh38"/"weights"
wsnp = set(pd.concat([pd.read_csv(Wt/f"weights.hm3_noMHC.{c}.l2.ldscore.gz", sep="\t", usecols=["SNP"]) for c in range(1, 23)]).SNP)
BLD = pd.concat([pd.read_csv(BL/f"baselineLD.{c}.l2.ldscore.gz", sep="\t") for c in range(1, 23)]).drop_duplicates("SNP"); BLD = BLD[BLD.SNP.isin(wsnp)].reset_index(drop=True)
Xb1 = np.column_stack([np.ones(len(BLD)), BLD[[c for c in BLD.columns if c not in ("CHR", "SNP", "BP")]].to_numpy(np.float64)])
def ldtab(d, gname):
    T_ = pd.concat([pd.read_csv(d/f"{gname}.{c}.l2.ldscore.gz", sep="\t") for c in range(1, 23)]).drop_duplicates("SNP").set_index("SNP")
    return T_[[c for c in T_.columns if c.endswith("L2")]]   # CHR·BP 열 제외
QL = ldtab(S3/"annot"/"Q5", "Q5"); vrows = []; VIFMAX = {}
for g in sorted(set(OWN[t] for t in ALL)):
    parts = [ldtab(Y/"annot"/f"L5_{g}_{m}", f"L5_{g}_{m}")[[f"L5_{g}_{m}L2"]] for m in LAY] + [ldtab(L/"annot"/f"A5_{g}", f"A5_{g}"), QL]
    Lx = pd.concat(parts, axis=1).loc[BLD.SNP]; names = list(Lx.columns); Lx = Lx.to_numpy(np.float64)
    Rz = Lx - Xb1 @ np.linalg.lstsq(Xb1, Lx, rcond=None)[0]; pr = np.corrcoef(Rz, rowvar=False); v = np.diag(np.linalg.inv(pr)); VIFMAX[g] = float(v.max())
    vrows.append(dict(lineage=g, **dict(zip(names, v)))); P(f"  {g:11s} VIF " + " ".join(f"{n.replace('L2', '')} {x:.2f}" for n, x in zip(names, v)))
    qi = [i for i, n in enumerate(names) if n.startswith("Q_")]; ei = [i for i, n in enumerate(names) if n.startswith("L5_")]
    P(f"  {'':11s} 효과–Q5 편상관 최대 |r| {np.abs(pr[np.ix_(ei, qi)]).max():.2f}")
pd.DataFrame(vrows).to_csv(EXP/"results"/"s03_vif.csv", index=False)
def own(t, R, kind="E"): v, d = R[kind][OWNL[t]]; return v, d
# ---------------- S-H1 (주) ----------------
H1 = ["AD", "EA", "SCZ", "BIP", "NEURO", "BMI"]
P("\n=== S-H1 (주): EQ에서 자기 층위(CAGE) 효과 τ* > 0 (단측, 6개 BH); 비교 E(T3all) ===")
pv, ret = [], []
for t in H1:
    v, d = own(t, RES[(t, "EQ")]); s = jse(d); pv.append(float(stats.norm.sf(v / s))); e0 = own(t, REF[(t, "E")])[0]; ret.append(v / e0)
q = bh(pv); n = int((q < .05).sum())
for t, v_p, q_, r_ in zip(H1, pv, q, ret):
    v, d = own(t, RES[(t, "EQ")]); e0, d0 = own(t, REF[(t, "E")])
    P(f"  {KO[t]:6s} ({OWN[t]}) E {e0:+.2f} ± {jse(d0):.2f} → EQ {v:+.2f} ± {jse(d):.2f} (유지 {r_ * 100:.0f} %, P {fp(v_p)}, q {fp(q_)})" + (f" — VIF 최대 {VIFMAX[OWN[t]]:.1f}" if VIFMAX[OWN[t]] > 5 else ""))
    SUM.append(dict(test="S-H1", trait=t, layer=OWNL[t], E=e0, tau=v, se=jse(d), p=v_p, q=q_, retention=r_))
v1 = "세포형 eQTL을 넘어 유지" if n >= 5 else ("부분" if n >= 3 else "대부분 세포형 eQTL로 설명")
P(f"  → FDR < 0.05: {n}/6, 유지 비율 중앙값 {np.median(ret) * 100:.0f} % → 판정: {v1}"); SUM.append(dict(test="S-H1 판정", tau=n, retention=float(np.median(ret)), trait=v1))
# ---------------- S-H2 (주, AD) ----------------
P("\n=== S-H2 (주, AD): EAQ에서 골수계 CAGE 활성도 τ* > 0 (단측 α 0.05); 비교 EA ===")
a, d = RES[("AD", "EAQ")]["A"]["CAGE"]; s = jse(d); p = float(stats.norm.sf(a / s)); a0, d0 = REF[("AD", "EA")]["A"]["CAGE"]
v2 = "AD의 골수계 전사 개시 활성 위치 정보는 미세아교·단핵구·뇌 세포형 eQTL을 넘는다" if p < .05 else "세포형 eQTL과 구분되지 않는다"
P(f"  AD CAGE 활성도 EA {a0:+.2f} ± {jse(d0):.2f} → EAQ {a:+.2f} ± {s:.2f} (유지 {a / a0 * 100:.0f} %, 단측 P {fp(p)}) → {v2}")
eff, de = RES[("AD", "EAQ")]["E"]["CAGE"]; P(f"  (참고) 같은 모형의 AD CAGE 효과 τ* {eff:+.2f} ± {jse(de):.2f}; 활성도 DNase {RES[('AD', 'EAQ')]['A']['DNASE'][0]:+.2f}, H3K27ac {RES[('AD', 'EAQ')]['A']['H3K27ac'][0]:+.2f}")
SUM.append(dict(test="S-H2", trait="AD", layer="CAGE activity", E=a0, tau=a, se=s, p=p, retention=a / a0, verdict=v2))
# ---------------- S-H3 ----------------
P("\n=== S-H3: 자기 층위 효과 τ* (지도 11 + 검증 3); E → EQ, EA → EAQ ===")
for t in ALL:
    e0 = own(t, REF[(t, "E")])[0]; v, d = own(t, RES[(t, "EQ")]); ea0 = own(t, REF[(t, "EA")])[0]; w, dw = own(t, RES[(t, "EAQ")])
    pe, pw = float(stats.norm.sf(v / jse(d))), float(stats.norm.sf(w / jse(dw)))
    P(f"  {KO[t]:6s} {OWNL[t]:6s} E {e0:+.2f} → EQ {v:+.2f} ± {jse(d):.2f} ({v / e0 * 100:.0f} %, P {fp(pe)}) | EA {ea0:+.2f} → EAQ {w:+.2f} ± {jse(dw):.2f} ({w / ea0 * 100 if abs(ea0) > 1e-9 else np.nan:.0f} %, P {fp(pw)})")
    SUM.append(dict(test="S-H3", trait=t, layer=OWNL[t], E=e0, tau=v, se=jse(d), p=pe, EA=ea0, tau_EAQ=w, se_EAQ=jse(dw), p_EAQ=pw))
PRq = {t: prof(RES[(t, "EQ")]["E"]) for t in ALL}; PRe = {t: prof(REF[(t, "E")]["E"]) for t in ALL}
Qe, Qq = Qtest(PRe, MAP), Qtest(PRq, MAP); assert abs(Qe[0] - 70.1) < .05, Qe   # 3p E(T3all) 재현
P(f"  층위 지도 Q (지도 11개, 자유도 20): E {Qe[0]:.1f} (P {fp(Qe[2])}) → EQ {Qq[0]:.1f} (P {fp(Qq[2])}, F 근사 {fp(Qq[3])})"); SUM.append(dict(test="S-H3 Q", E=Qe[0], tau=Qq[0], p=Qq[2]))
for lab, gr, k, side in CON:
    e1, s1, p1 = comb(PRe, gr, k, side); e2, s2, p2 = comb(PRq, gr, k, side)
    P(f"  {lab}: E {e1:+.2f} ± {s1:.2f} (P {fp(p1)}) → EQ {e2:+.2f} ± {s2:.2f} (P {fp(p2)}) [{'단측' if side == 1 else '양측, α 0.0167'}]"); SUM.append(dict(test=f"S-H3 {lab}", E=e1, tau=e2, se=s2, p=p2))
# ---------------- S-H4 ----------------
P("\n=== S-H4: 세포형 MaxCPP 주석 τ* (EQ; 기준선 GTEx·BLUEPRINT MaxCPP와 예측 효과 3층위를 함께 넣은 조건부 값; 단측 P) ===")
for t in TR:
    R = RES[(t, "EQ")]["Q"]; P(f"  {KO[t]:9s} " + " | ".join(f"{x[2:]} {R[x][0]:+.2f} ± {jse(R[x][1]):.2f} (P {fp(stats.norm.sf(R[x][0] / jse(R[x][1])))})" for x in QN))
    for x in QN: SUM.append(dict(test="S-H4", trait=t, layer=x, tau=R[x][0], se=jse(R[x][1]), p=float(stats.norm.sf(R[x][0] / jse(R[x][1])))))
# ---------------- S-H5 ----------------
P("\n=== S-H5: 더 큰 GWAS 6개 (기술; S-H1·S-H2와 같은 값) ===")
for t in [f"N_{x}" for x in NEW] + ["N_AD"]:
    e0 = own(t, REF[(t, "E")])[0]; v, d = own(t, RES[(t, "EQ")]); pe = float(stats.norm.sf(v / jse(d)))
    extra = ""
    if t == "N_AD":
        a, da = RES[(t, "EAQ")]["A"]["CAGE"]; a0 = REF[(t, "EA")]["A"]["CAGE"][0]; extra = f" | CAGE 활성도 EA {a0:+.2f} → EAQ {a:+.2f} ± {jse(da):.2f} (P {fp(stats.norm.sf(a / jse(da)))})"
    P(f"  {KO[t]:9s} {OWNL[t]:6s} E {e0:+.2f} → EQ {v:+.2f} ± {jse(d):.2f} ({v / e0 * 100:.0f} %, P {fp(pe)})" + extra)
    SUM.append(dict(test="S-H5", trait=t, layer=OWNL[t], E=e0, tau=v, se=jse(d), p=pe))
pd.DataFrame(SUM).to_csv(EXP/"results"/"s03_summary.csv", index=False)
(EXP/"results"/"s03_main_tests.txt").write_text("\n".join(log), encoding="utf-8")
