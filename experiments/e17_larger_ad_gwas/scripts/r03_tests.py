"""실험 3r r03: 새 AD GWAS(N_AD) 검정 (PROTOCOL_exp3r_2026-10-25, 2026-10-25 고정).
 τ*와 공동 jackknife는 3j·3l·3p·3q와 같은 식(3m m_common, 3p p_common). 원래 AD(Bellenguez) 값은 기존 실행에서 읽어 나란히 적는다(차이 검정은 하지 않음).
 R-H1 (주) HL10 골수계 Δ = τ*(H) − τ*(Lo), CAGE·DNase 각각 단측 α 0.025:
        P < 0.025 → "검출" / 아니고 95 % 구간 상한(Δ + 1.96 SE) < +2.28(3q 키 DNase Δ) → "키 크기의 효과는 배제" / 그 밖 → "미해결(검정력 부족)"
 R-H2 자기 층위(CAGE) 효과 τ*: EA, EAc, S, SA, EAM (단측 P, 유지 비율 = τ*/E). R-H3 HL20 Δ (R-H1과 같은 세 범주로 기술).
 R-H4 지도 위치: 12계통 τ*·순위, 층위 τ*와 c1·c2(재추출 평균), AD 대 키 c1(원래 키·3n 새 키), AD를 바꾼 지도의 Q(주 지도판, 3n 갱신 지도판). R-H0 정보량(T3all 로그).
출력: results/r03_tests.txt, r03_summary.csv
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
EXP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXP.parent/"e12_diagnostics"/"scripts")); sys.path.insert(0, str(EXP.parent/"e15_activity_constructions"/"scripts"))
import m_common as MC
import p_common as PC
from m_common import SEEDS, LAY, MAP, M_TOT, jse
RR = Path(f"{DATA}/exp3r/runs"); JR = MC.RUNS; NR = Path(f"{DATA}/exp3n/runs"); Qd = Path(f"{DATA}/exp3q"); C = Path(f"{DATA}/exp3c")
LIN12 = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
NEW5 = ["HEIGHT", "BMI", "SCZ", "BIP", "LDL"]; HEIGHT_D = 2.28   # 3q 키 DNase Δ (HL10)
log, SUM = [], []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
def fp(p): return f"{p:.2g}"
strip = lambda n: re.sub(r"L2_\d+$", "", n)
# ---------------- R-H0 정보량 ----------------
def info(path_log, path_npz):
    s = open(path_log, encoding="utf-8").read(); g = lambda pat: re.search(pat, s)
    h2 = g(r"Total Observed scale h2: ([\d.eE-]+) \(([\d.eE-]+)\)"); ic = g(r"Intercept: ([\d.eE-]+) \(([\d.eE-]+)\)"); mc = g(r"Mean Chi\^2: ([\d.eE-]+)"); rt = g(r"Ratio: ([\d.eE-]+) \(([\d.eE-]+)\)")
    return dict(h2=float(h2.group(1)), h2_se=float(h2.group(2)), Z=float(h2.group(1)) / float(h2.group(2)), mean_chi2=float(mc.group(1)), intercept=float(ic.group(1)), intercept_se=float(ic.group(2)),
                ratio=(float(rt.group(1)) if rt else np.nan), n_snp=int(np.load(path_npz)["n_snp"]))
io_, in_ = info(JR/"AD__T3all.log", JR/"AD__T3all.npz"), info(RR/"N_AD__T3all.log", RR/"N_AD__T3all.npz")
P("=== R-H0 정보량 (T3all 모형 로그; 원래 Bellenguez → 새 GWAS) ===")
P(f"  h² {io_['h2']:.4f} ({io_['h2_se']:.4f}), Z {io_['Z']:.1f} → {in_['h2']:.4f} ({in_['h2_se']:.4f}), Z {in_['Z']:.1f} | 평균 χ² {io_['mean_chi2']:.3f} → {in_['mean_chi2']:.3f}"
  f" | 절편 {io_['intercept']:.3f} ({io_['intercept_se']:.3f}) → {in_['intercept']:.3f} ({in_['intercept_se']:.3f}) | 비율 {io_['ratio']:.3f} → {in_['ratio']:.3f} | 회귀 SNP {io_['n_snp']:,} → {in_['n_snp']:,}")
SUM += [dict(test="R-H0", gwas=k, **v) for k, v in [("original", io_), ("new", in_)]]
# ---------------- R-H1 (주), R-H3 ----------------
SDQ = pd.read_csv(Qd/"annot_sd_common.csv", index_col=0).sd_common
def load_hl(path):
    z = np.load(path); h2 = float(z["h2_log"]); res = {}
    for k, n in enumerate(strip(x) for x in z["names"]):
        if n not in SDQ.index: continue
        kind = "H" if n.startswith("H") else "Lo"; m = [x for x in LAY if n.endswith("_" + x)][0]; f = SDQ[n] * M_TOT / h2
        assert m not in res.get(kind, {}); res.setdefault(kind, {})[m] = (z["coef"][k] * f, z["part_delete"][:, k] * f)
    assert {k: len(v) for k, v in res.items()} == {"H": 3, "Lo": 3}; return res
def verdict(d, s, p):
    if p < .025: return "검출"
    return "키 크기의 효과는 배제" if d + 1.96 * s < HEIGHT_D else "미해결(검정력 부족)"
for model, tag in [("HL10", "R-H1 (주)"), ("HL20", "R-H3 (민감도)")]:
    Ro, Rn = load_hl(Qd/"runs"/f"AD__{model}.npz"), load_hl(RR/f"N_AD__{model}.npz")
    P(f"\n=== {tag}: {model} 골수계 Δ = τ*(H) − τ*(Lo) (단측; 판정 α 0.025, 키 크기 기준 +{HEIGHT_D}) ===")
    for m in LAY:
        do, sdo = Ro["H"][m][0] - Ro["Lo"][m][0], jse(Ro["H"][m][1] - Ro["Lo"][m][1])
        dn, sdn = Rn["H"][m][0] - Rn["Lo"][m][0], jse(Rn["H"][m][1] - Rn["Lo"][m][1]); pn = float(stats.norm.sf(dn / sdn))
        v = verdict(dn, sdn, pn) if m != "H3K27ac" else "(기술)"
        P(f"  {m:7s} 원래 Δ {do:+.2f} ± {sdo:.2f} (P {fp(stats.norm.sf(do / sdo))}) → 새 H {Rn['H'][m][0]:+.2f} ± {jse(Rn['H'][m][1]):.2f}, Lo {Rn['Lo'][m][0]:+.2f} ± {jse(Rn['Lo'][m][1]):.2f}"
          f" → Δ {dn:+.2f} ± {sdn:.2f} (95 % 구간 {dn - 1.96 * sdn:+.2f} ~ {dn + 1.96 * sdn:+.2f}, 단측 P {fp(pn)}) → {v}")
        SUM.append(dict(test=f"{tag} {model}", layer=m, old=do, old_se=sdo, new=dn, new_se=sdn, p=pn, lo95=dn - 1.96 * sdn, hi95=dn + 1.96 * sdn, verdict=v))
# ---------------- R-H2 ----------------
P("\n=== R-H2: 자기 층위(CAGE) 효과 τ* — 활성도 조건화 방식별 (단측 P, 유지 비율 = τ*/E; 원래 → 새) ===")
Eo = PC.load("AD", "E")["E"]["CAGE"][0]; En = PC.load("N_AD", "T3all", runs=RR)["E"]["CAGE"]
P(f"  E   : 원래 {Eo:+.2f} → 새 {En[0]:+.2f} ± {jse(En[1]):.2f} (P {fp(stats.norm.sf(En[0] / jse(En[1])))})")
SUM.append(dict(test="R-H2", model="E", old=Eo, new=En[0], new_se=jse(En[1]), p=float(stats.norm.sf(En[0] / jse(En[1])))))
for model in ["EA", "EAc", "S", "SA", "EAM"]:
    o = PC.load("AD", model)["E"]["CAGE"]; R = PC.load("N_AD", model, runs=RR); n = R["E"]["CAGE"]; s = jse(n[1]); p = float(stats.norm.sf(n[0] / s))
    act = R.get("A", {}); acts = (" | 활성도 " + ", ".join(f"{m} {act[m][0]:+.2f} ± {jse(act[m][1]):.2f}" for m in LAY if m in act)) if act else ""
    P(f"  {model:4s}: 원래 {o[0]:+.2f} ± {jse(o[1]):.2f} ({o[0] / Eo * 100:.0f} %) → 새 {n[0]:+.2f} ± {s:.2f} ({n[0] / En[0] * 100:.0f} %, 단측 P {fp(p)}){acts}")
    SUM.append(dict(test="R-H2", model=model, old=o[0], old_se=jse(o[1]), new=n[0], new_se=s, p=p, retention=n[0] / En[0], **{f"act_{m}": act[m][0] for m in act}, **{f"act_{m}_se": jse(act[m][1]) for m in act}))
# ---------------- R-H4 지도 위치 ----------------
P("\n=== R-H4: 지도 위치 (기술) ===")
SDC = pd.read_csv(C/"annot_sd_common.csv", index_col=0).sd_common
def tau1(path, name):
    z = np.load(path); h2 = float(z["h2_log"]); names = [strip(x) for x in z["names"]]; k = names.index(name); f = SDC[name] * M_TOT / h2
    return z["coef"][k] * f, jse(z["part_delete"][:, k] * f)
LR = pd.DataFrame([dict(lineage=g, **dict(zip(["tau", "se"], tau1(RR/f"N_AD__M1_{g}.npz", f"L5_{g}")))) for g in LIN12]).sort_values("tau", ascending=False).reset_index(drop=True)
S3C = pd.read_csv(EXP.parent/"e02_lineage_map"/"results"/"z04_sldsc.csv"); S3C = S3C[(S3C.model == "M1") & (S3C.trait == "AD")].sort_values("tau_star", ascending=False).reset_index(drop=True)
P("  계통 지도 τ* (새, 높은 순): " + "; ".join(f"{r.lineage} {r.tau:+.2f} ± {r.se:.2f}" for r in LR.itertuples()))
P("  원래 상위 3: " + "; ".join(f"{r.lineage} {r.tau_star:+.2f}" for r in S3C.head(3).itertuples()) + f" | 골수계 순위 원래 {int(S3C.index[S3C.lineage == 'myeloid'][0]) + 1} → 새 {int(LR.index[LR.lineage == 'myeloid'][0]) + 1}")
for r in LR.itertuples(): SUM.append(dict(test="R-H4 계통", model=r.lineage, new=r.tau, new_se=r.se, p=float(stats.norm.sf(r.tau / r.se))))
Lo_, Ln_ = MC.layers("AD", SEEDS, JR), MC.layers("N_AD", SEEDS, RR)
PRo, PRn = MC.prof("AD", SEEDS, JR), MC.prof("N_AD", SEEDS, RR)
P("  층위 τ* (재추출 평균): " + " | ".join(f"{m} {Lo_[m][0]:+.2f} → {Ln_[m][0]:+.2f} ± {jse(Ln_[m][1]):.2f}" for m in LAY))
P(f"  c1 (CAGE − DNase) {PRo[0][0]:+.2f} ± {jse(PRo[1][:, 0]):.2f} → {PRn[0][0]:+.2f} ± {jse(PRn[1][:, 0]):.2f} | c2 (CAGE − H3K27ac) {PRo[0][1]:+.2f} → {PRn[0][1]:+.2f} ± {jse(PRn[1][:, 1]):.2f}")
for m in LAY: SUM.append(dict(test="R-H4 층위", model=m, old=Lo_[m][0], new=Ln_[m][0], new_se=jse(Ln_[m][1])))
SUM += [dict(test="R-H4 c1", old=PRo[0][0], new=PRn[0][0], new_se=jse(PRn[1][:, 0])), dict(test="R-H4 c2", old=PRo[0][1], new=PRn[0][1], new_se=jse(PRn[1][:, 1]))]
PR_main = {t: MC.prof(t, SEEDS, JR) for t in MAP}; PR_h3n = MC.prof("HEIGHT", SEEDS, NR)
for lab, pa, ph in [("원래 AD 대 원래 키", PR_main["AD"], PR_main["HEIGHT"]), ("새 AD 대 원래 키", PRn, PR_main["HEIGHT"]), ("새 AD 대 3n 새 키", PRn, PR_h3n)]:
    e, s, p, _ = MC.comb({"AD": pa, "HEIGHT": ph}, [(1, ["AD"]), (-1, ["HEIGHT"])], 0)
    P(f"  {lab} (c1, 공동 jackknife): {e:+.2f} ± {s:.2f}, 양측 P {fp(p)}"); SUM.append(dict(test="R-H4 AD 대 키", model=lab, new=e, new_se=s, p=p))
Q0 = MC.Qtest(PR_main, MAP); assert abs(Q0[0] - 74.3) < .05, Q0
PR_a = dict(PR_main, AD=PRn); PR_b = {t: (MC.prof(t, SEEDS, NR) if t in NEW5 else PR_main[t]) for t in MAP}; PR_b["AD"] = PRn
for lab, PRx in [("주 지도 (원래 14 중 지도 11)", PR_main), ("AD만 새 GWAS", PR_a), ("3n 갱신 지도 + 새 AD", PR_b)]:
    Q, df, p, pf = MC.Qtest(PRx, MAP); P(f"  Q {lab}: {Q:.1f} (자유도 {df}, P {fp(p)}, F 근사 {fp(pf)})"); SUM.append(dict(test="R-H4 Q", model=lab, new=Q, p=p, p_F=pf))
pd.DataFrame(SUM).to_csv(EXP/"results"/"r03_summary.csv", index=False)
(EXP/"results"/"r03_tests.txt").write_text("\n".join(log), encoding="utf-8")
