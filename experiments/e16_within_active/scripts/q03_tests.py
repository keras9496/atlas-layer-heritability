"""실험 3q q03: 활성 요소 안 효과 정보 검정 (PROTOCOL_exp3q_2026-10-23, 2026-10-23 고정; 사후 분석).
 τ*와 공동 jackknife는 3l l04·3p p03과 같은 식(3p p_common의 함수). Δ_m = τ*(H_m) − τ*(Lo_m), delete 값도 같은 차이.
 Q-H1 (주) HL10에서 자기 층위 Δ > 0, 단측, 지도 형질 11개 BH → FDR < 0.05가 6개 이상 "활성 요소 안에서 효과 정보 있음", 3–5 "부분", ≤ 2 "검출 안 됨".
 Q-H2 3l 유지군 대 설명군 자기 층위 Δ 역분산 가중 평균 차이(양측). Q-H3 Δ의 층위 대비(Δ_CAGE − Δ_DNASE, Δ_CAGE − Δ_H3K27ac) Q(11)·사전 대비·검증 대비.
 Q-H4 HL20 같은 검정. Q-H5 새 GWAS 다섯 형질. 기술: 계통별 기준선 보정 LD 점수 VIF(VIF > 5 계통의 형질은 Q-H1에서 해석하지 않음).
 사용: python q03_tests.py dry | main
출력: results/q03_<단계>_tests.txt, q03_main_summary.csv, q03_main_vif.csv, q03_values.csv
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"e15_activity_constructions"/"scripts"))
from p_common import X, Y, LAY, M_TOT, OWN, MAP, VAL, ALL, NEW, OWNL, KO, CON, jcov, jse, prof, Qtest, comb, bh
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; Qd = Path(f"{DATA}/exp3q")
RUN = Qd/"test" if STAGE == "dry" else Qd/"runs"
SDQ = pd.read_csv(Qd/"annot_sd_common.csv", index_col=0).sd_common
log, SUM, VAL_ROWS = [], [], []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
strip = lambda n: re.sub(r"L2_\d+$", "", n)
def load(name, model):
    """H/Lo → 층위 → (τ*, delete 값)"""
    z = np.load(RUN/f"{name}__{model}.npz"); h2 = float(z["h2_log"]); res = {}
    for k, n in enumerate(strip(x) for x in z["names"]):
        if n not in SDQ.index: continue
        kind = "H" if n.startswith("H") else "Lo"; m = [x for x in LAY if n.endswith("_" + x)][0]; f = SDQ[n] * M_TOT / h2
        assert m not in res.get(kind, {}); res.setdefault(kind, {})[m] = (z["coef"][k] * f, z["part_delete"][:, k] * f)
    assert {k: len(v) for k, v in res.items()} == {"H": 3, "Lo": 3}, (name, model, {k: len(v) for k, v in res.items()})
    return res
def delta(R): return {m: (R["H"][m][0] - R["Lo"][m][0], R["H"][m][1] - R["Lo"][m][1]) for m in LAY}
if STAGE == "dry":
    for model in ["HL10", "HL20"]:
        R = load("FAKE1", model); P(f"[dry] FAKE1 {model}: " + " | ".join(f"{k} " + ", ".join(f"{m} {v[0]:+.2f}" for m, v in R[k].items()) for k in R))
    P("형식 확인 통과 (가짜 형질 — 결과 아님)"); (EXP/"results"/"q03_dry_tests.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit()
RES = {(t, m): load(t, m) for t in ALL for m in ["HL10", "HL20"]}; RESN = {(t, m): load(f"N_{t}", m) for t in NEW for m in ["HL10", "HL20"]}
for (t, m), R in list(RES.items()) + [((f"N_{t}", m), R) for (t, m), R in RESN.items()]:
    for kind in R:
        for x, (tau, d) in R[kind].items(): VAL_ROWS.append(dict(trait=t, model=m, kind=kind, layer=x, tau=float(tau), se=jse(d)))
    for x, (dv, dd) in delta(R).items(): VAL_ROWS.append(dict(trait=t, model=m, kind="Delta", layer=x, tau=float(dv), se=jse(dd)))
pd.DataFrame(VAL_ROWS).to_csv(EXP/"results"/"q03_values.csv", index=False)
# ---------------- 기술: VIF ----------------
P("=== 기술: 기준선 보정 LD 점수 VIF (H 3 + Lo 3; 회귀 SNP, 기준선 97개를 뺀 잔차) ===")
BL = X/"GRCh38"/"baselineLD_v2.2"; Wt = X/"GRCh38"/"weights"
wsnp = set(pd.concat([pd.read_csv(Wt/f"weights.hm3_noMHC.{c}.l2.ldscore.gz", sep="\t", usecols=["SNP"]) for c in range(1, 23)]).SNP)
BLD = pd.concat([pd.read_csv(BL/f"baselineLD.{c}.l2.ldscore.gz", sep="\t") for c in range(1, 23)]).drop_duplicates("SNP"); BLD = BLD[BLD.SNP.isin(wsnp)].reset_index(drop=True)
Xb1 = np.column_stack([np.ones(len(BLD)), BLD[[c for c in BLD.columns if c not in ("CHR", "SNP", "BP")]].to_numpy(np.float64)])
VIF = {}; vrows = []
for g in sorted(set(OWN.values())):
    for k in ["10", "20"]:
        names = [f"H{k}_{g}_{m}" for m in LAY] + [f"Lo{k}_{g}_{m}" for m in LAY]; d = Qd/"annot"/f"HL{k}_{g}"
        Lx = pd.concat([pd.read_csv(d/f"HL{k}_{g}.{c}.l2.ldscore.gz", sep="\t") for c in range(1, 23)]).drop_duplicates("SNP").set_index("SNP").loc[BLD.SNP, [f"{n}L2" for n in names]].to_numpy(np.float64)
        Rz = Lx - Xb1 @ np.linalg.lstsq(Xb1, Lx, rcond=None)[0]; pr = np.corrcoef(Rz, rowvar=False); v = np.diag(np.linalg.inv(pr))
        VIF[(g, k)] = float(v.max()); vrows.append(dict(lineage=g, cut=k, **dict(zip(names, v)), **{f"pr_HLo_{m}": pr[i, i + 3] for i, m in enumerate(LAY)}))
        P(f"  {g:11s} {k}% VIF " + " ".join(f"{x:.2f}" for x in v) + " | 같은 층위 H–Lo 편상관 " + " ".join(f"{pr[i, i + 3]:.2f}" for i in range(3)))
pd.DataFrame(vrows).to_csv(EXP/"results"/"q03_main_vif.csv", index=False)
# ---------------- Q-H1 (주), Q-H4 ----------------
def h1(model, tag):
    P(f"\n=== {tag}: {model}에서 자기 층위 Δ = τ*(H) − τ*(Lo) > 0 (단측, 지도 형질 11개 BH) ===")
    pv = []
    for t in MAP: dv, dd = delta(RES[(t, model)])[OWNL[t]]; pv.append(float(stats.norm.sf(dv / jse(dd))))
    q = bh(pv); n = int((q < .05).sum()); k = model[2:]
    for t, p_, q_ in zip(MAP, pv, q):
        m = OWNL[t]; R = RES[(t, model)]; dv, dd = delta(R)[m]; vf = VIF[(OWN[t], k)]
        P(f"  {KO[t]:6s} {m:7s} H {R['H'][m][0]:+.2f} ± {jse(R['H'][m][1]):.2f}, Lo {R['Lo'][m][0]:+.2f} ± {jse(R['Lo'][m][1]):.2f} → Δ {dv:+.2f} ± {jse(dd):.2f} (P {p_:.2g}, q {q_:.2g})"
          + (" — VIF > 5, 해석하지 않음" if vf > 5 else ""))
        SUM.append(dict(test=f"{tag} {model}", trait=t, layer=m, tau_H=R["H"][m][0], tau_Lo=R["Lo"][m][0], delta=dv, se=jse(dd), p=p_, q=q_))
    v = "활성 요소 안에서 효과 정보 있음" if n >= 6 else ("부분" if n >= 3 else "검출 안 됨")
    P(f"  → FDR < 0.05: {n}/11 → 판정: {v}"); SUM.append(dict(test=f"{tag} {model} 판정", delta=n, trait=v)); return dict(zip(MAP, q))
q10 = h1("HL10", "Q-H1 (주)")
# ---------------- Q-H2 ----------------
P("\n=== Q-H2: 3l 유지군 대 설명군 (자기 층위 Δ, 역분산 가중 평균 차이, 양측) ===")
RG = ["EA", "SCZ", "BIP", "LDL", "HEIGHT"]; XG = ["AD", "MONO", "IBD", "RA"]
for model in ["HL10", "HL20"]:
    PRd = {t: (np.array([delta(RES[(t, model)])[OWNL[t]][0]]), delta(RES[(t, model)])[OWNL[t]][1].reshape(-1, 1)) for t in MAP}
    e1, s1, _ = comb(PRd, [(1, RG)], 0, 1); e2, s2, _ = comb(PRd, [(1, XG)], 0, 1); e, s, p = comb(PRd, [(1, RG), (-1, XG)], 0, 2)
    P(f"  {model}: 유지군 Δ {e1:+.2f} ± {s1:.2f}, 설명군 Δ {e2:+.2f} ± {s2:.2f}, 차이 {e:+.2f} ± {s:.2f} (P {p:.2g})"); SUM.append(dict(test=f"Q-H2 {model}", delta=e, se=s, p=p))
# ---------------- Q-H3 ----------------
P("\n=== Q-H3: 활성 요소 안 효과의 층위 대비 (Δ_CAGE − Δ_DNASE, Δ_CAGE − Δ_H3K27ac; 공동 jackknife) ===")
for model in ["HL10", "HL20"]:
    PR = {t: prof(delta(RES[(t, model)])) for t in ALL}; Q, df, p, pf = Qtest(PR, MAP)
    P(f"  [{model}] Q (지도 형질 11개) {Q:.1f}, 자유도 {df}, P {p:.2g} (F 근사 {pf:.2g})"); SUM.append(dict(test=f"Q-H3 {model} Q", delta=Q, p=p, q=pf))
    for lab, gr, k, side in CON:
        e, s, pp = comb(PR, gr, k, side); P(f"      {lab}: {e:+.2f} ± {s:.2f}, P {pp:.2g} [{'단측' if side == 1 else '양측, α 0.0167'}]"); SUM.append(dict(test=f"Q-H3 {model} {lab}", delta=e, se=s, p=pp))
    P("      형질별 Δ C/D/H: " + "; ".join(f"{KO[t]} " + "/".join(f"{delta(RES[(t, model)])[m][0]:+.2f}" for m in LAY) for t in ALL))
# ---------------- Q-H4 ----------------
h1("HL20", "Q-H4 (민감도)")
# ---------------- Q-H5 ----------------
P("\n=== Q-H5: 새 GWAS 다섯 형질 (자기 층위 Δ, 단측 P) ===")
for t in NEW:
    m = OWNL[t]; row = []
    for model in ["HL10", "HL20"]:
        R = RESN[(t, model)]; dv, dd = delta(R)[m]; s = jse(dd); p_ = float(stats.norm.sf(dv / s)); row.append(f"{model} H {R['H'][m][0]:+.2f} Lo {R['Lo'][m][0]:+.2f} Δ {dv:+.2f} ± {s:.2f} (P {p_:.2g})")
        SUM.append(dict(test=f"Q-H5 {model}", trait=t, layer=m, tau_H=R["H"][m][0], tau_Lo=R["Lo"][m][0], delta=dv, se=s, p=p_))
    P(f"  {KO[t]:6s} {m}: " + " | ".join(row))
pd.DataFrame(SUM).to_csv(EXP/"results"/"q03_main_summary.csv", index=False)
(EXP/"results"/"q03_main_tests.txt").write_text("\n".join(log), encoding="utf-8")
