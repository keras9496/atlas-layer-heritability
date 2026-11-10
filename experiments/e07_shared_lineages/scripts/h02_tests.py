"""실험 3h h02: 계통 고정 비교 검정 (PROTOCOL_exp3h_2026-10-13).
 형질 t × 계통 g의 T3all(전체 트랙, 단일 실행): τ*와 200 블록 삭제값 → c1 = CAGE − DNase, c2 = CAGE − H3K27ac, jackknife SE (3f·3g와 같다).
 H-H1 (주): DNase형과 CAGE형 형질이 함께 농축된 계통마다 Δ_g = IVW(CAGE형 c1) − IVW(DNase형 c1); 계통 결합 Δ̄ (IVW), 보수적 SE(완전 상관) 단측 p < 0.05 그리고 Δ_g > 0 계통 ≥ 4/6.
 H-H2: 형질별 농축 계통에서 c1 부호가 유형과 맞는 비율(기술). H-H3: 농축 형질 3개 이상인 계통의 (c1, c2) 다변량 Q(기술).
 사용: python h02_tests.py dry | main
출력: results/h02_<단계>_profiles.csv, h02_<단계>_tests.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re, sys, glob
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; R0 = Path(__file__).resolve().parents[2]
X = Path(f"{DATA}/exp3"); Y = Path(f"{DATA}/exp3d"); F = Path(f"{DATA}/exp3f"); H = Path(f"{DATA}/exp3h")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LAY = ["CAGE", "DNASE", "H3K27ac"]
SD = pd.read_csv(Y/"annot_sd_common.csv", index_col=0).sd_common
B0 = X/"GRCh38"/"baselineLD_v2.2"; M_TOT = sum(float(open(B0/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
OWN = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
       "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal", "BMD": "stromal", "NEUT": "myeloid", "TC": "liver"}
DTYPE = ["HEIGHT", "BMD"]; CTYPE = ["AD", "MONO", "EA", "SCZ", "BIP", "BMI", "LDL", "TC"]   # 프로토콜 유형 (미정: NEURO, IBD, RA, NEUT)
KO = {"AD": "AD", "MONO": "단핵구 수", "EA": "교육 연수", "SCZ": "조현병", "BIP": "양극성장애", "NEURO": "신경증", "BMI": "BMI", "IBD": "IBD", "RA": "RA", "LDL": "LDL", "HEIGHT": "키",
      "BMD": "골밀도", "NEUT": "호중구 수", "TC": "총콜레스테롤", "FAKE": "가짜"}
def path(t, g):
    if STAGE == "dry": return H/"test"/"sldsc"/f"{t}__T3all_{g}"
    return F/"sldsc"/f"{t}__T3all" if g == OWN[t] else H/"sldsc"/f"{t}__T3all_{g}"
def load(t, g):
    out = path(t, g); R = pd.read_csv(f"{out}.results", sep="\t"); R["name"] = R.Category.str.replace(r"L2_\d+$", "", regex=True)
    h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log", encoding="utf-8").read()).group(1))
    D = pd.read_csv(f"{out}.part_delete", sep=r"\s+", header=None).to_numpy(); res = {}
    for m in LAY:
        n = f"L5_{g}_{m}"; k = int(R.index[R.name == n][0]); f = SD[n] * M_TOT / h2; res[m] = (R.Coefficient[k] * f, D[:, k] * f)
    return res
def jse(dl): B = len(dl); return np.sqrt((B - 1) / B * ((dl - dl.mean()) ** 2).sum())
def jcov(dls): B = dls.shape[0]; z = dls - dls.mean(0); return (B - 1) / B * z.T @ z
def enriched():
    m = pd.read_csv(R0/"e02_lineage_map"/"results"/"z04_sldsc.csv"); m = m[m.model == "M1"]
    E = {t: sorted(g.lineage[g.tau_p1 < .05]) for t, g in m.groupby("trait")}
    for t in ["BMD", "NEUT", "TC"]:
        E[t] = sorted(f.split("__M1_")[1][:-8] for f in glob.glob(str(F/"sldsc"/f"{t}__M1_*.results"))
                      if stats.norm.sf(pd.read_csv(f, sep="\t")["Coefficient_z-score"].iloc[-1]) < .05)
    return E
E = {"FAKE": ["adipose"]} if STAGE == "dry" else enriched()
PR, rows = {}, []
for t, gs in E.items():
    for g in gs:
        r = load(t, g); est = {m: r[m][0] for m in LAY}; dl = {m: r[m][1] for m in LAY}
        c1, d1 = est["CAGE"] - est["DNASE"], dl["CAGE"] - dl["DNASE"]; c2, d2 = est["CAGE"] - est["H3K27ac"], dl["CAGE"] - dl["H3K27ac"]
        PR[(t, g)] = dict(d=np.array([c1, c2]), S=jcov(np.column_stack([d1, d2])))
        rows.append(dict(trait=t, lineage=g, own=(g == OWN.get(t)), **{f"{m}_tau": est[m] for m in LAY}, **{f"{m}_se": jse(dl[m]) for m in LAY},
                         top=max(est, key=est.get), c1=c1, c1_se=jse(d1), c2=c2, c2_se=jse(d2)))
L = pd.DataFrame(rows); L.to_csv(EXP/"results"/f"h02_{STAGE}_profiles.csv", index=False)
P("=== 형질 × 농축 계통의 3-층위 프로필 (T3all, 전체 트랙, 단일 실행; τ* ± SE; * = 자기 계통) ===")
for r in L.sort_values(["trait", "lineage"]).itertuples():
    P(f"  {KO[r.trait]:7s} {r.lineage:22s}{'*' if r.own else ' '} " + "  ".join(f"{m} {getattr(r, m + '_tau'):+5.2f} ± {getattr(r, m + '_se'):.2f}" for m in LAY)
      + f"  → 1위 {r.top:7s} | c1 {r.c1:+5.2f} ± {r.c1_se:.2f}")
if STAGE == "dry":
    (EXP/"results"/"h02_dry_tests.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit()
I = L.set_index(["trait", "lineage"])
def ivw(keys):
    w = np.array([1 / I.loc[k, "c1_se"] ** 2 for k in keys]); v = np.array([I.loc[k, "c1"] for k in keys]); return (w * v).sum() / w.sum(), np.sqrt(1 / w.sum())

P("\n=== H-H1 (주): 같은 계통에서 CAGE형 c1 − DNase형 c1 (Δ_g > 0 예측) ===")
LG = sorted({g for g in set(L.lineage) if any(g in E[t] for t in DTYPE) and any(g in E[t] for t in CTYPE)})
D_ = []
for g in LG:
    kd = [(t, g) for t in DTYPE if g in E[t]]; kc = [(t, g) for t in CTYPE if g in E[t]]
    md, sd_ = ivw(kd); mc, sc = ivw(kc); dg, sg = mc - md, np.hypot(sd_, sc); D_.append((g, dg, sg))
    P(f"  {g:22s} CAGE형 [{', '.join(KO[t] for t, _ in kc)}] {mc:+5.2f} ± {sc:.2f}  대  DNase형 [{', '.join(KO[t] for t, _ in kd)}] {md:+5.2f} ± {sd_:.2f}"
      f"  → Δ {dg:+5.2f} ± {sg:.2f} (단측 p {stats.norm.sf(dg / sg):.3g})")
assert len(LG) == 6, LG   # 프로토콜 표와 같아야 한다
w = np.array([1 / s ** 2 for _, _, s in D_]); dv = np.array([d for _, d, _ in D_]); se_i = np.array([s for _, _, s in D_])
Db = (w * dv).sum() / w.sum(); se_ind = np.sqrt(1 / w.sum()); se_con = (w * se_i).sum() / w.sum()
p_con, p_ind = stats.norm.sf(Db / se_con), stats.norm.sf(Db / se_ind); npos = int((dv > 0).sum())
P(f"  결합 Δ̄ = {Db:+.2f}; 보수적 SE {se_con:.2f} → 단측 p {p_con:.3g} (주 판정) | 독립 가정 SE {se_ind:.2f} → 단측 p {p_ind:.3g} (기술)")
P(f"  Δ_g > 0 계통 {npos}/6")
if p_con < .05 and npos >= 4: verdict = "유지"
elif npos >= 4: verdict = "방향만 유지"
else: verdict = "유지 안 됨"
P(f"  → 판정: {verdict}")

P("\n=== H-H2 (기술): 형질별 농축 계통에서 c1 부호가 유형과 맞는가 (CAGE형 c1 > 0, DNase형 c1 < 0) ===")
for t in DTYPE + CTYPE + ["NEURO", "IBD", "RA", "NEUT"]:
    sub = L[L.trait == t]
    if t in DTYPE: ok = int((sub.c1 < 0).sum()); lab = "DNase형"
    elif t in CTYPE: ok = int((sub.c1 > 0).sum()); lab = "CAGE형"
    else: ok = None; lab = "미정"
    tops = ", ".join(f"{r.lineage}:{r.top}({r.c1:+.2f})" for r in sub.itertuples())
    P(f"  {KO[t]:7s} {lab:6s} " + (f"일치 {ok}/{len(sub)}" if ok is not None else f"c1>0 {int((sub.c1 > 0).sum())}/{len(sub)}") + f"  | {tops}")

P("\n=== H-H3 (기술): 계통 안 이질성 — 농축 형질 3개 이상인 계통, (c1, c2) 다변량 Cochran Q ===")
for g in sorted(set(L.lineage)):
    ts = sorted(t for t in E if g in E[t])
    if len(ts) < 3: continue
    d = {t: PR[(t, g)]["d"] for t in ts}; W = {t: np.linalg.inv(PR[(t, g)]["S"]) for t in ts}
    dbar = np.linalg.solve(sum(W.values()), sum(W[t] @ d[t] for t in ts)); Q = float(sum((d[t] - dbar) @ W[t] @ (d[t] - dbar) for t in ts)); df = 2 * (len(ts) - 1)
    P(f"  {g:22s} 형질 {len(ts)}개 [{', '.join(KO[t] for t in ts)}]  Q {Q:5.1f}, 자유도 {df}, p {stats.chi2.sf(Q, df):.3g}")
(EXP/"results"/f"h02_{STAGE}_tests.txt").write_text("\n".join(log), encoding="utf-8")
