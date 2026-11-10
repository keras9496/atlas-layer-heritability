"""실험 3p p03: 활성도 구성별 검정 (PROTOCOL_exp3p_2026-10-23, 2026-10-23 고정).
 모형: E(3j T3all), EA·EAM(3l), EAc·S·SA(3p). 새 GWAS 다섯 형질: E·EA(3n), EAc·S·SA(3p, N_ 접두어).
 P-H1 (주) EAc·S 각각에서 자기 층위 효과 τ* > 0, 단측, 지도 형질 11개 BH → FDR < 0.05가 6개 이상 "활성도 조건부로 유지", 3–5 "부분", ≤ 2 "대부분 활성도로 설명".
 P-H2 (주) 3l 양상 재현: 유지군 R = {EA, SCZ, BIP, LDL, HEIGHT}의 4/5 이상 FDR < 0.05, 설명군 X = {AD, MONO, IBD, RA}의 1/4 이하 FDR < 0.05 → "재현".
       EAc·S 둘 다 재현 → "활성도 구성에 강건", 하나 → "부분", 없음 → "구성에 의존".
 P-H3 유지 비율(자기 층위 효과 τ*, 모형 k / E). P-H4 모형별 대비 c1·c2, Q(11), 사전 대비 3개(양측 α 0.05/3), 검증 대비(단측). P-H5 새 GWAS 다섯 형질.
 기술: 계통별 기준선 보정 LD 점수 VIF (EAc 효과 3 + 연속 활성도 3, S 층화 효과 3, SA 층화 효과 3 + 활성도 3). VIF > 5인 계통의 형질은 P-H1·P-H2에서 해석하지 않는다(판정 계산에는 넣는다).
 사용: python p03_tests.py dry | main
출력: results/p03_<단계>_tests.txt, p03_main_summary.csv, p03_main_vif.csv
"""
import sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
sys.path.insert(0, str(Path(__file__).resolve().parent))
from p_common import *
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]
log, SUM = [], []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
if STAGE == "dry":
    for model, want in [("EAc", {"E": 3, "A": 3}), ("S", {"E": 3}), ("SA", {"E": 3, "A": 3})]:
        R = load("FAKE1", model, runs=O/"test"); P(f"[dry] FAKE1 {model}: " + " | ".join(f"{k} " + ", ".join(f"{m} {v[0]:+.2f} ± {jse(v[1]):.2f}" for m, v in R[k].items()) for k in sorted(R)))
        assert {k: len(v) for k, v in R.items()} == want, "주석 종류·층위 수가 맞지 않다"
    P("형식 확인 통과 (가짜 형질 — 결과 아님)"); (EXP/"results"/"p03_dry_tests.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit()
# ---------------- 불러오기와 기준값 확인 ----------------
RES = {(t, m): load(t, m) for t in ALL for m in MODELS}
RESN = {(t, m): load(t, m, "new") for t in NEW for m in ["E", "EA", "EAc", "S", "SA"]}
assert abs(RES[("AD", "E")]["E"]["CAGE"][0] - 2.91) < .006, RES[("AD", "E")]["E"]["CAGE"][0]          # 3l l04 기준
assert abs(RES[("AD", "EA")]["A"]["CAGE"][0] - 1.91) < .006
for t in ALL:
    for m in MODELS:
        want = {"E": {"E"}, "EA": {"E", "A"}, "EAM": {"E", "A", "M"}, "EAc": {"E", "A"}, "S": {"E"}, "SA": {"E", "A"}}[m]
        assert set(RES[(t, m)]) == want and all(len(v) == 3 for v in RES[(t, m)].values()), (t, m)
QE = Qtest({t: prof(RES[(t, "E")]["E"]) for t in ALL}, MAP)[0]; QEA = Qtest({t: prof(RES[(t, "EA")]["E"]) for t in ALL}, MAP)[0]
P(f"기준값 확인: E(전체 트랙) Q {QE:.1f} (원고 70.1), EA Q {QEA:.1f} (3l 53.1)"); assert abs(QE - 70.1) < .1 and abs(QEA - 53.1) < .1
# ---------------- 기술: VIF ----------------
P("=== 기술: 기준선 보정 LD 점수 VIF (회귀 SNP, 기준선 97개를 뺀 잔차) ===")
BL = X/"GRCh38"/"baselineLD_v2.2"; Wt = X/"GRCh38"/"weights"
wsnp = set(pd.concat([pd.read_csv(Wt/f"weights.hm3_noMHC.{c}.l2.ldscore.gz", sep="\t", usecols=["SNP"]) for c in range(1, 23)]).SNP)
BLD = pd.concat([pd.read_csv(BL/f"baselineLD.{c}.l2.ldscore.gz", sep="\t") for c in range(1, 23)]).drop_duplicates("SNP"); BLD = BLD[BLD.SNP.isin(wsnp)].reset_index(drop=True)
Xb1 = np.column_stack([np.ones(len(BLD)), BLD[[c for c in BLD.columns if c not in ("CHR", "SNP", "BP")]].to_numpy(np.float64)])
def ldmat(d, stem, names):
    Lx = pd.concat([pd.read_csv(d/f"{stem}.{c}.l2.ldscore.gz", sep="\t") for c in range(1, 23)]).drop_duplicates("SNP").set_index("SNP")
    return Lx.loc[BLD.SNP, [f"{n}L2" for n in names]].to_numpy(np.float64)
def vif(Lm):
    Rz = Lm - Xb1 @ np.linalg.lstsq(Xb1, Lm, rcond=None)[0]; pr = np.corrcoef(Rz, rowvar=False); return np.diag(np.linalg.inv(pr)), pr
VIF = {}; vrows = []
for g in sorted(set(OWN.values())):
    Le = np.hstack([ldmat(Y/"annot"/f"L5_{g}_{m}", f"L5_{g}_{m}", [f"L5_{g}_{m}"]) for m in LAY])
    Lac = ldmat(O/"annot"/f"Ac3_{g}", f"Ac3_{g}", [f"Ac_{g}_{m}" for m in LAY]); Ls = ldmat(O/"annot"/f"S5_{g}", f"S5_{g}", [f"S5_{g}_{m}" for m in LAY])
    La = ldmat(L/"annot"/f"A5_{g}", f"A5_{g}", [f"A5_{g}_{m}" for m in LAY])
    v1, p1 = vif(np.hstack([Le, Lac])); v2, _ = vif(Ls); v3, p3 = vif(np.hstack([Ls, La]))
    VIF[g] = {"EAc": v1.max(), "S": v2.max(), "SA": v3.max()}
    vrows.append(dict(lineage=g, **{f"EAc_{k}": x for k, x in zip([f"E_{m}" for m in LAY] + [f"Ac_{m}" for m in LAY], v1)}, **{f"S_{m}": x for m, x in zip(LAY, v2)},
                      **{f"SA_{k}": x for k, x in zip([f"S_{m}" for m in LAY] + [f"A_{m}" for m in LAY], v3)}, **{f"pr_EAc_{m}": p1[i, i + 3] for i, m in enumerate(LAY)}, **{f"pr_SA_{m}": p3[i, i + 3] for i, m in enumerate(LAY)}))
    P(f"  {g:11s} EAc VIF " + " ".join(f"{x:.2f}" for x in v1) + f" | S VIF " + " ".join(f"{x:.2f}" for x in v2) + f" | SA VIF " + " ".join(f"{x:.2f}" for x in v3)
      + " | 같은 층위 편상관 EAc " + " ".join(f"{p1[i, i + 3]:.2f}" for i in range(3)) + " / SA " + " ".join(f"{p3[i, i + 3]:.2f}" for i in range(3)))
pd.DataFrame(vrows).to_csv(EXP/"results"/"p03_main_vif.csv", index=False)
# ---------------- P-H1 (주) ----------------
RG = ["EA", "SCZ", "BIP", "LDL", "HEIGHT"]; XG = ["AD", "MONO", "IBD", "RA"]
QS = {}
for model in ["EA", "EAc", "S", "SA", "EAM"]:
    tag = "주" if model in ("EAc", "S") else ("참고, 3l" if model in ("EA", "EAM") else "기술")
    P(f"\n=== P-H1 ({tag}): {model}에서 자기 층위 효과 τ* > 0 (단측, 지도 형질 11개 BH) ===")
    pv = []
    for t in MAP: e, d = RES[(t, model)]["E"][OWNL[t]]; pv.append(float(stats.norm.sf(e / jse(d))))
    q = bh(pv); QS[model] = dict(zip(MAP, q)); n = int((q < .05).sum())
    for t, p_, q_ in zip(MAP, pv, q):
        m = OWNL[t]; e, d = RES[(t, model)]["E"][m]; e0 = RES[(t, "E")]["E"][m][0]; vf = VIF[OWN[t]].get(model, np.nan)
        P(f"  {KO[t]:6s} {m:7s} E {e0:+.2f} → {model} {e:+.2f} ± {jse(d):.2f} (P {p_:.2g}, q {q_:.2g}) | 유지 비율 {e / e0:.0%}" + (f" | VIF 최대 {vf:.2f}" + (" — 해석하지 않음" if vf > 5 else "") if model in VIF[OWN[t]] else ""))
        SUM.append(dict(test=f"P-H1 {model}", trait=t, layer=m, tau_E=e0, tau=e, se=jse(d), p=p_, q=q_, retention=e / e0))
    v = "활성도 조건부로 유지" if n >= 6 else ("부분" if n >= 3 else "대부분 활성도로 설명")
    P(f"  → FDR < 0.05: {n}/11 → 등급: {v}"); SUM.append(dict(test=f"P-H1 {model} 판정", tau=n, trait=v))
# ---------------- P-H2 (주) ----------------
P("\n=== P-H2 (주): 3l 계통 양상 재현 (유지군 4/5 이상, 설명군 1/4 이하가 FDR < 0.05) ===")
rep = {}
for model in ["EAc", "S", "SA", "EA"]:
    nr = sum(QS[model][t] < .05 for t in RG); nx = sum(QS[model][t] < .05 for t in XG); ok = nr >= 4 and nx <= 1; rep[model] = ok
    P(f"  {model:4s} 유지군 {nr}/5 ({', '.join(KO[t] for t in RG if QS[model][t] < .05)}) | 설명군 {nx}/4 ({', '.join(KO[t] for t in XG if QS[model][t] < .05) or '없음'}) → {'재현' if ok else '재현 안 됨'}"
      + (" (참고: 3l 자체)" if model == "EA" else " (기술)" if model == "SA" else ""))
    SUM.append(dict(test=f"P-H2 {model}", tau=nr, se=nx, trait="재현" if ok else "재현 안 됨"))
nrep = rep["EAc"] + rep["S"]; v2 = "활성도 구성에 강건" if nrep == 2 else ("부분" if nrep == 1 else "구성에 의존")
P(f"  → P-H2 판정 (EAc·S): {v2}"); SUM.append(dict(test="P-H2 판정", trait=v2))
# ---------------- P-H3 ----------------
P("\n=== P-H3: 자기 층위 효과의 유지 비율 (모형 / E) ===")
for model in ["EA", "EAc", "S", "SA", "EAM"]:
    rr = [RES[(t, model)]["E"][OWNL[t]][0] / RES[(t, "E")]["E"][OWNL[t]][0] for t in MAP]
    P(f"  {model:4s} 중앙값 {np.median(rr):.0%} | 유지군 " + ", ".join(f"{KO[t]} {RES[(t, model)]['E'][OWNL[t]][0] / RES[(t, 'E')]['E'][OWNL[t]][0]:.0%}" for t in RG)
      + " | 설명군 " + ", ".join(f"{KO[t]} {RES[(t, model)]['E'][OWNL[t]][0] / RES[(t, 'E')]['E'][OWNL[t]][0]:.0%}" for t in XG))
P("  활성도 τ* (EAc 연속 활성도; EA 이진 활성도 참고) — " + "; ".join(f"{KO[t]} " + "/".join(f"{RES[(t, 'EAc')]['A'][m][0]:+.2f}" for m in LAY) + " (EA " + "/".join(f"{RES[(t, 'EA')]['A'][m][0]:+.2f}" for m in LAY) + ")" for t in MAP))
# ---------------- P-H4 ----------------
P("\n=== P-H4: 모형별 효과 층위 대비 (c1 = CAGE − DNase, c2 = CAGE − H3K27ac; 공동 jackknife) ===")
for model in MODELS:
    PR = {t: prof(RES[(t, model)]["E"]) for t in ALL}; Q, df, p, pf = Qtest(PR, MAP)
    P(f"  [{model}] Q (지도 형질 11개) {Q:.1f}, 자유도 {df}, P {p:.2g} (F 근사 {pf:.2g})"); SUM.append(dict(test=f"P-H4 {model} Q", tau=Q, p=p, q=pf))
    for lab, gr, k, side in CON:
        e, s, pp = comb(PR, gr, k, side); P(f"      {lab}: {e:+.2f} ± {s:.2f}, P {pp:.2g} [{'단측' if side == 1 else '양측, α 0.0167'}]"); SUM.append(dict(test=f"P-H4 {model} {lab}", tau=e, se=s, p=pp))
    P("      형질별 c1 (± se) / 1위 효과 층위: " + "; ".join(f"{KO[t]} {PR[t][0][0]:+.2f} ± {jse(PR[t][1][:, 0]):.2f} {max(LAY, key=lambda m: RES[(t, model)]['E'][m][0])}" for t in ALL))
# ---------------- P-H5 ----------------
P("\n=== P-H5: 새 GWAS 다섯 형질 (자기 층위 효과 τ*, 단측 P, 유지 비율; c1) ===")
for t in NEW:
    m = OWNL[t]; e0 = RESN[(t, "E")]["E"][m][0]; row = []
    for model in ["E", "EA", "EAc", "S", "SA"]:
        e, d = RESN[(t, model)]["E"][m]; s = jse(d); p_ = float(stats.norm.sf(e / s)); c1 = prof(RESN[(t, model)]["E"]); row.append(f"{model} {e:+.2f} ± {s:.2f} (P {p_:.2g}, {e / e0:.0%}; c1 {c1[0][0]:+.2f})")
        SUM.append(dict(test=f"P-H5 {model}", trait=t, layer=m, tau_E=e0, tau=e, se=s, p=p_, retention=e / e0))
    P(f"  {KO[t]:6s} {m}: " + " | ".join(row))
pd.DataFrame(SUM).to_csv(EXP/"results"/"p03_main_summary.csv", index=False)
(EXP/"results"/"p03_main_tests.txt").write_text("\n".join(log), encoding="utf-8")
