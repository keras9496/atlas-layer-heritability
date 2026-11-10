"""실험 3i i03: 세포형 일치 층위 지도 검정 (PROTOCOL_exp3i_2026-10-13).
 형질 × 자기 계통의 M3 (단일 실행): τ*, 200 블록 삭제값 → c1 = CAGE − DNase, c2 = CAGE − H3K27ac, jackknife SE (3f·3g·3h와 같다).
 I-H1 (주): 주 지도 확정 8개 형질(CAGE형 7 + 키)의 1위 층위가 같은가 — 8개 중 7개 이상 그리고 키 DNase·AD CAGE.
 I-H2: 지도 형질 11개 (c1, c2) 다변량 Q. I-H3: 3f 사전 대비 3개(양측, α 0.05/3 기준 보고). I-H4: 검증 형질 1위 층위와 검증 대비(단측). 모두 기술.
 사용: python i03_tests.py dry | main
출력: results/i03_<단계>_layers.csv, i03_<단계>_tests.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); I = Path(f"{DATA}/exp3i")
SL = I/"test"/"sldsc" if STAGE == "dry" else I/"sldsc"
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LAY = ["CAGE", "DNASE", "H3K27ac"]
SD = pd.read_csv(I/"annot_sd_common.csv", index_col=0).sd_common
B0 = X/"GRCh38"/"baselineLD_v2.2"; M_TOT = sum(float(open(B0/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
OWN = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
       "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal", "BMD": "stromal", "NEUT": "myeloid", "TC": "liver"}
MAP = ["AD", "MONO", "EA", "SCZ", "BIP", "NEURO", "BMI", "IBD", "RA", "LDL", "HEIGHT"]; VAL = ["BMD", "NEUT", "TC"]
MAIN_TOP = {t: "CAGE" for t in MAP}; MAIN_TOP["HEIGHT"] = "DNASE"                       # 3f 시드 평균 1위 (11개 모두)
CONF = ["AD", "MONO", "EA", "SCZ", "BIP", "BMI", "LDL", "HEIGHT"]                          # 3f 확정 8개
PRED = {"BMD": "DNASE", "NEUT": "CAGE", "TC": "CAGE"}
KO = {"AD": "AD", "MONO": "단핵구 수", "EA": "교육 연수", "SCZ": "조현병", "BIP": "양극성장애", "NEURO": "신경증", "BMI": "BMI", "IBD": "IBD", "RA": "RA", "LDL": "LDL", "HEIGHT": "키",
      "BMD": "골밀도", "NEUT": "호중구 수", "TC": "총콜레스테롤", "FAKE": "가짜"}
TR = ["FAKE"] if STAGE == "dry" else MAP + VAL
if STAGE == "dry": OWN = {"FAKE": "myeloid"}
def jse(dl): B = len(dl); return np.sqrt((B - 1) / B * ((dl - dl.mean()) ** 2).sum())
def jcov(dls): B = dls.shape[0]; z = dls - dls.mean(0); return (B - 1) / B * z.T @ z
PR, rows = {}, []
for t in TR:
    g = OWN[t]; out = SL/f"{t}__M3_{g}"; R = pd.read_csv(f"{out}.results", sep="\t"); R["name"] = R.Category.str.replace(r"L2_\d+$", "", regex=True)
    h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log", encoding="utf-8").read()).group(1))
    D = pd.read_csv(f"{out}.part_delete", sep=r"\s+", header=None).to_numpy(); est, dl = {}, {}
    for m in LAY:
        n = f"M5_{g}_{m}"; k = int(R.index[R.name == n][0]); f = SD[n] * M_TOT / h2; est[m] = R.Coefficient[k] * f; dl[m] = D[:, k] * f
    order = sorted(LAY, key=lambda m: -est[m]); x1, x2 = order[0], order[1]; dx = dl[x1] - dl[x2]; p12 = stats.norm.sf((est[x1] - est[x2]) / jse(dx))
    c1, d1 = est["CAGE"] - est["DNASE"], dl["CAGE"] - dl["DNASE"]; c2, d2 = est["CAGE"] - est["H3K27ac"], dl["CAGE"] - dl["H3K27ac"]
    PR[t] = dict(d=np.array([c1, c2]), S=jcov(np.column_stack([d1, d2])), d1=d1, d2=d2)
    rows.append(dict(trait=t, lineage=g, **{f"{m}_tau": est[m] for m in LAY}, **{f"{m}_se": jse(dl[m]) for m in LAY}, top=x1, second=x2, top_vs_second_p1=p12,
                     main_top=MAIN_TOP.get(t, PRED.get(t, "")), c1=c1, c1_se=jse(d1), c2=c2, c2_se=jse(d2)))
L = pd.DataFrame(rows); L.to_csv(EXP/"results"/f"i03_{STAGE}_layers.csv", index=False)
P("=== 세포형 일치 3-층위 프로필 (M3, 단일 실행; τ* ± SE; * = 1위가 2위보다 단측 p < 0.05) ===")
for r in L.itertuples():
    P(f"  {KO[r.trait]:7s} {r.lineage:11s} " + "  ".join(f"{m} {getattr(r, m + '_tau'):+5.2f} ± {getattr(r, m + '_se'):.2f}" for m in LAY)
      + f"  → 1위 {r.top}{'*' if r.top_vs_second_p1 < .05 else ' '} (주 지도/예측 {r.main_top}) | c1 {r.c1:+5.2f} ± {r.c1_se:.2f}, c2 {r.c2:+5.2f} ± {r.c2_se:.2f}")
if STAGE == "dry":
    (EXP/"results"/"i03_dry_tests.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit()
T = L.set_index("trait")
same = [t for t in CONF if T.loc[t, "top"] == MAIN_TOP[t]]; diff = [f"{KO[t]} {MAIN_TOP[t]}→{T.loc[t, 'top']}" for t in CONF if t not in same]
h_ok, a_ok = T.loc["HEIGHT", "top"] == "DNASE", T.loc["AD", "top"] == "CAGE"
verdict = "유지" if (len(same) >= 7 and h_ok and a_ok) else "유지 안 됨"
P(f"\n=== I-H1 (주): 확정 8개 형질의 1위 층위 유지 {len(same)}/8, 키 {T.loc['HEIGHT', 'top']}, AD {T.loc['AD', 'top']} → 판정: {verdict}")
if diff: P("  바뀐 형질: " + ", ".join(diff))
P("  미정 3개(판정 밖): " + ", ".join(f"{KO[t]} {T.loc[t, 'top']}" for t in ["NEURO", "IBD", "RA"]))
P(f"  1위가 2위보다 단측 p < 0.05인 지도 형질: {int((T.loc[MAP, 'top_vs_second_p1'] < .05).sum())}/11")
d = {t: PR[t]["d"] for t in MAP}; W = {t: np.linalg.inv(PR[t]["S"]) for t in MAP}
dbar = np.linalg.solve(sum(W.values()), sum(W[t] @ d[t] for t in MAP)); Q = float(sum((d[t] - dbar) @ W[t] @ (d[t] - dbar) for t in MAP)); df = 2 * (len(MAP) - 1)
P(f"\n=== I-H2 (기술): 지도 형질 11개 다변량 Cochran Q = {Q:.1f}, 자유도 {df}, p = {stats.chi2.sf(Q, df):.3g}")
def ivw(ts, k):
    w = np.array([1 / T.loc[t, f"{k}_se"] ** 2 for t in ts]); v = np.array([T.loc[t, k] for t in ts]); return (w * v).sum() / w.sum(), np.sqrt(1 / w.sum())
P("\n=== I-H3 (기술): 3f 사전 대비 (양측, α 0.05/3 = 0.0167 기준) ===")
for lab, A, Bs, k in [("AD 대 키", ["AD"], ["HEIGHT"], "c1"), ("뇌 형질 4개 대 키", ["EA", "SCZ", "BIP", "NEURO"], ["HEIGHT"], "c1"), ("IBD·RA 대 AD·단핵구 수", ["IBD", "RA"], ["AD", "MONO"], "c2")]:
    ma, sa = ivw(A, k); mb, sb = ivw(Bs, k); df_, se_ = ma - mb, np.hypot(sa, sb); p = 2 * stats.norm.sf(abs(df_ / se_))
    P(f"  {lab} ({k}): {ma:+.2f} ± {sa:.2f} 대 {mb:+.2f} ± {sb:.2f} → 차이 {df_:+.2f} ± {se_:.2f}, p {p:.3g} {'(α 0.0167 미만)' if p < .0167 else ''}")
P("\n=== I-H4 (기술): 검증 형질 ===")
for t in VAL: P(f"  {KO[t]:7s} 1위 {T.loc[t, 'top']} (예측 {PRED[t]}) {'일치' if T.loc[t, 'top'] == PRED[t] else '불일치'}")
mc, sc = ivw(["NEUT", "TC"], "c1"); md, sd_ = ivw(["BMD"], "c1"); dv, sv = mc - md, np.hypot(sc, sd_)
P(f"  검증 대비 (c1): CAGE형 예측 {mc:+.2f} ± {sc:.2f} 대 DNase형 예측 {md:+.2f} ± {sd_:.2f} → 차이 {dv:+.2f} ± {sv:.2f}, 단측 p {stats.norm.sf(dv / sv):.3g}")
(EXP/"results"/f"i03_{STAGE}_tests.txt").write_text("\n".join(log), encoding="utf-8")
