"""실험 3k k03: 근위·원위 층위 지도 검정 (PROTOCOL_exp3k_2026-10-15).
 형질 × 자기 계통의 PD 모형(단일 실행, 정렬 블록): 구간 k ∈ {D 원위, P 근위}마다 층위 τ*, c1_k = CAGE_k − DNase_k, c2_k = CAGE_k − H3K27ac_k.
 K-H1 (주, 원위): (a) CAGE형 7개 중 CAGE_D 1위 수 ≥ 5 그리고 키 DNase_D 1위, (b) 지도 형질 11개 (c1_D, c2_D) 다변량 Q (공동 jackknife, χ²) P < 0.05.
   둘 다 → 원위에서 유지 / 하나 → 부분 / 둘 다 아님 → 원위에서 유지 안 됨.
 K-H2 원위 대비 3개(양측, α 0.05/3)와 검증 대비(단측) / K-H3 근위의 같은 검정 / K-H4 원위 1위 대 2위(단측), F 근사, 주석 크기 — 모두 기술.
 대비의 묶음 평균은 이 분석의 SE로 IVW, 대비 SE는 가중 결합한 delete 값의 jackknife SE(3j와 같다).
 사용: python k03_tests.py dry | main
출력: results/k03_<단계>_layers.csv, k03_<단계>_tests.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); K = Path(f"{DATA}/exp3k")
RUN = K/"test"/"runs" if STAGE == "dry" else K/"runs"
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LAY = ["CAGE", "DNASE", "H3K27ac"]; SH = {"CAGE": "CAGE", "DNASE": "DNase", "H3K27ac": "H3K27ac"}
SD = pd.read_csv(K/"annot_sd_common.csv", index_col=0).sd_common
B0 = X/"GRCh38"/"baselineLD_v2.2"; M_TOT = sum(float(open(B0/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
OWN = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
       "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal", "BMD": "stromal", "NEUT": "myeloid", "TC": "liver"}
MAP = ["AD", "MONO", "EA", "SCZ", "BIP", "NEURO", "BMI", "IBD", "RA", "LDL", "HEIGHT"]; VAL = ["BMD", "NEUT", "TC"]
CTYPE = ["AD", "MONO", "EA", "SCZ", "BIP", "BMI", "LDL"]   # 주 지도 확정 CAGE형 7개
KO = {"AD": "AD", "MONO": "단핵구 수", "EA": "교육 연수", "SCZ": "조현병", "BIP": "양극성장애", "NEURO": "신경증", "BMI": "BMI", "IBD": "IBD", "RA": "RA", "LDL": "LDL", "HEIGHT": "키",
      "BMD": "골밀도", "NEUT": "호중구 수", "TC": "총콜레스테롤", "FAKE1": "가짜"}
TR = ["FAKE1"] if STAGE == "dry" else MAP + VAL
if STAGE == "dry": OWN["FAKE1"] = "myeloid"
strip = lambda n: re.sub(r"L2_\d+$", "", n)
def jcov(Z): B = Z.shape[0]; z = Z - Z.mean(0); return (B - 1) / B * z.T @ z
def jse(v): return float(np.sqrt(jcov(np.asarray(v).reshape(-1, 1))[0, 0]))
PR, rows = {}, []
for t in TR:
    g = OWN[t]; z = np.load(RUN/f"{t}__PD_{g}.npz"); h2 = float(z["h2_log"]); names = [strip(n) for n in z["names"]]
    est, dl = {}, {}
    for k, n in enumerate(names):
        if n in SD.index:
            f = SD[n] * M_TOT / h2; cls, m = n.split("_")[0], n.split("_")[-1]; est[(cls, m)] = z["coef"][k] * f; dl[(cls, m)] = z["part_delete"][:, k] * f
    assert len(est) == 6, (t, list(est))
    rec = dict(trait=t, lineage=g)
    for cls in ["D", "P"]:
        e = {m: est[(cls, m)] for m in LAY}; o = sorted(LAY, key=lambda m: -e[m]); d12 = dl[(cls, o[0])] - dl[(cls, o[1])]
        rec.update({f"{cls}_{m}_tau": e[m] for m in LAY}); rec.update({f"{cls}_{m}_se": jse(dl[(cls, m)]) for m in LAY})
        rec[f"{cls}_top"] = o[0]; rec[f"{cls}_top_vs_second_p1"] = stats.norm.sf((e[o[0]] - e[o[1]]) / jse(d12))
        c1, d1 = e["CAGE"] - e["DNASE"], dl[(cls, "CAGE")] - dl[(cls, "DNASE")]; c2, d2 = e["CAGE"] - e["H3K27ac"], dl[(cls, "CAGE")] - dl[(cls, "H3K27ac")]
        rec.update({f"{cls}_c1": c1, f"{cls}_c1_se": jse(d1), f"{cls}_c2": c2, f"{cls}_c2_se": jse(d2)}); PR[(t, cls)] = (np.array([c1, c2]), np.column_stack([d1, d2]))
    rows.append(rec)
L = pd.DataFrame(rows); L.to_csv(EXP/"results"/f"k03_{STAGE}_layers.csv", index=False)
for cls, lab in [("D", "원위 (> 1 kb)"), ("P", "근위 (≤ 1 kb)")]:
    P(f"=== {lab}: 층위 τ* ± SE (PD 모형 6주석 동시; * = 1위가 2위보다 단측 p < 0.05) ===")
    for r in L.itertuples():
        P(f"  {KO[r.trait]:7s} {r.lineage:11s} " + "  ".join(f"{SH[m]} {getattr(r, f'{cls}_{m}_tau'):+5.2f} ± {getattr(r, f'{cls}_{m}_se'):.2f}" for m in LAY)
          + f"  → 1위 {SH[getattr(r, f'{cls}_top')]}{'*' if getattr(r, f'{cls}_top_vs_second_p1') < .05 else ' '} | c1 {getattr(r, f'{cls}_c1'):+5.2f} ± {getattr(r, f'{cls}_c1_se'):.2f}")
if STAGE == "dry":
    (EXP/"results"/"k03_dry_tests.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit()
T = L.set_index("trait")
def Qtest(ts, cls):
    d = np.concatenate([PR[(t, cls)][0] for t in ts]); S = jcov(np.hstack([PR[(t, cls)][1] for t in ts])); Kn = len(ts)
    Si = np.linalg.inv(S); Xm = np.tile(np.eye(2), (Kn, 1)); mu = np.linalg.solve(Xm.T @ Si @ Xm, Xm.T @ Si @ d); r = d - Xm @ mu
    Q = float(r @ Si @ r); df = 2 * (Kn - 1); nu = 199; Fs = (nu - df + 1) / (nu * df) * Q; return Q, df, stats.chi2.sf(Q, df), stats.f.sf(Fs, df, nu - df + 1)
def comb(groups, cls, k):
    a, keys = [], []
    for sgn, ts in groups:
        w = np.array([1 / jse(PR[(t, cls)][1][:, k]) ** 2 for t in ts]); a += list(sgn * w / w.sum()); keys += ts
    est = float(sum(ai * PR[(t, cls)][0][k] for ai, t in zip(a, keys))); se = jse(sum(ai * PR[(t, cls)][1][:, k] for ai, t in zip(a, keys))); return est, se
CON = [("AD 대 키 (c1)", [(1, ["AD"]), (-1, ["HEIGHT"])], 0), ("뇌 형질 4개 대 키 (c1)", [(1, ["EA", "SCZ", "BIP", "NEURO"]), (-1, ["HEIGHT"])], 0),
       ("IBD·RA 대 AD·단핵구 수 (c2)", [(1, ["IBD", "RA"]), (-1, ["AD", "MONO"])], 1)]
VCON = ("검증 대비: 호중구 수·총콜레스테롤 대 골밀도 (c1)", [(1, ["NEUT", "TC"]), (-1, ["BMD"])], 0)
res = {}
for cls, lab in [("D", "원위"), ("P", "근위")]:
    n_c = sum(T.loc[t, f"{cls}_top"] == "CAGE" for t in CTYPE); h_ok = T.loc["HEIGHT", f"{cls}_top"] == "DNASE"
    Q, df, pQ, pF = Qtest(MAP, cls); res[cls] = (n_c, h_ok, pQ)
    P(f"\n=== {lab}: 1위 층위 — CAGE형 7개 중 CAGE 1위 {n_c}/7 (" + ", ".join(f"{KO[t]} {SH[T.loc[t, f'{cls}_top']]}" for t in CTYPE) + f"), 키 {SH[T.loc['HEIGHT', f'{cls}_top']]}")
    P(f"  미정 3개: " + ", ".join(f"{KO[t]} {SH[T.loc[t, f'{cls}_top']]}" for t in ["NEURO", "IBD", "RA"]) + " | 검증: " + ", ".join(f"{KO[t]} {SH[T.loc[t, f'{cls}_top']]}" for t in VAL))
    P(f"  1위 > 2위 (단측 p < 0.05): {int((T.loc[MAP, f'{cls}_top_vs_second_p1'] < .05).sum())}/11")
    P(f"  다변량 Q (지도 형질 11개, 공동 jackknife) = {Q:.1f}, 자유도 {df}, p = {pQ:.3g} (F 근사 p {pF:.3g})")
    for l2, gr, k in CON:
        e, se = comb(gr, cls, k); p2 = 2 * stats.norm.sf(abs(e / se)); P(f"  {l2}: {e:+.2f} ± {se:.2f}, p {p2:.3g}{' (α 0.0167 미만)' if p2 < .0167 else ''}")
    e, se = comb(VCON[1], cls, VCON[2]); P(f"  {VCON[0]}: {e:+.2f} ± {se:.2f}, 단측 p {stats.norm.sf(e / se):.3g}")
n_c, h_ok, pQ = res["D"]; a_ok = n_c >= 5 and h_ok; b_ok = pQ < .05
verdict = "원위에서 유지" if (a_ok and b_ok) else ("부분" if (a_ok or b_ok) else "원위에서 유지 안 됨")
P(f"\n=== K-H1 (주) 판정: (a) {'충족' if a_ok else '불충족'} (CAGE형 {n_c}/7, 키 {'DNase' if h_ok else '다른 층위'}), (b) {'충족' if b_ok else '불충족'} (Q p {pQ:.3g}) → {verdict}")
(EXP/"results"/f"k03_{STAGE}_tests.txt").write_text("\n".join(log), encoding="utf-8")
