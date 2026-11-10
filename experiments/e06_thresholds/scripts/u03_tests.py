"""실험 3g u03: 역치별 층위 프로필과 검정 (PROTOCOL_exp3g_2026-10-11). 3f v03의 통계량을 단일 실행 jackknife로 계산한다.
 역치 q: L0p1 (0.1 %, 주 판정), L0p5, L1p0 (u02), L5 (5 % 참조 = 3f T3all).
 G-H1 다변량 Cochran Q (본 형질 11개, d = (c1, c2), 자유도 20) / G-H2 1위 층위가 3f 1단계 시드 평균 1위와 같은 형질 수 (≥ 9/11 그리고 키 = DNASE)
 G-H3 쌍 비교 L-H2a–c (양측, α 0.05/3) / G-H4 V-H1 (호중구 수·총콜레스테롤 c1 가중 평균 − 골밀도 c1 > 0, 단측)
 판정(0.1 %): 재현 = G-H1 통과 + G-H2 충족 / 방향만 재현 = G-H2 충족, G-H1 불통과 / 재현 안 됨 = G-H2 불충족. 나머지 역치는 기술.
 사용: python u03_tests.py dry | main
출력: results/u03_<단계>_layers.csv, u03_<단계>_tests.txt, figures/exp3g_threshold.png (main)
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); Y = Path(f"{DATA}/exp3d"); F = Path(f"{DATA}/exp3f"); G = Path(f"{DATA}/exp3g")
R3F = EXP.parent/"e05_layer_map"/"results"; (EXP/"figures").mkdir(exist_ok=True)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LAY = ["CAGE", "DNASE", "H3K27ac"]; TH = ["L0p1", "L0p5", "L1p0", "L5"]; THL = {"L0p1": "0.1 %", "L0p5": "0.5 %", "L1p0": "1 %", "L5": "5 %"}
SDF = G/("annot_sd_common.csv" if STAGE == "main" or (G/"annot_sd_common.csv").exists() else "test_annot_sd_common.csv")
SD = pd.concat([pd.read_csv(SDF, index_col=0).sd_common, pd.read_csv(Y/"annot_sd_common.csv", index_col=0).sd_common])
B0 = X/"GRCh38"/"baselineLD_v2.2"; M_TOT = sum(float(open(B0/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
MAIN = ["AD", "MONO", "EA", "SCZ", "BIP", "NEURO", "BMI", "IBD", "RA", "LDL", "HEIGHT"]; VAL = ["BMD", "NEUT", "TC"]
LINE = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
        "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal", "BMD": "stromal", "NEUT": "myeloid", "TC": "liver", "FAKE": "myeloid"}
KO = {"AD": "AD", "MONO": "단핵구 수", "EA": "교육 연수", "SCZ": "조현병", "BIP": "양극성장애", "NEURO": "신경증", "BMI": "BMI", "IBD": "IBD", "RA": "RA", "LDL": "LDL", "HEIGHT": "키",
      "BMD": "골밀도", "NEUT": "호중구 수", "TC": "총콜레스테롤", "FAKE": "가짜"}
TR = MAIN + VAL if STAGE == "main" else ["FAKE"]
def path(t, q):
    if q == "L5": return (F/"sldsc" if STAGE == "main" else F/"test"/"sldsc")/f"{t}__T3all"
    return (G/"sldsc" if STAGE == "main" else G/"test"/"sldsc")/f"{t}__G3_{q}"
def lay_of(n): return [m for m in LAY if n.endswith("_" + m)][0]
def load(t, q):
    out = path(t, q); R = pd.read_csv(f"{out}.results", sep="\t"); R["name"] = R.Category.str.replace(r"L2_\d+$", "", regex=True)
    h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log", encoding="utf-8").read()).group(1))
    D = pd.read_csv(f"{out}.part_delete", sep=r"\s+", header=None).to_numpy(); res = {}
    for k, n in enumerate(R.name):
        if n.startswith(q + "_"):
            f = SD[n] * M_TOT / h2; res[lay_of(n)] = dict(tau=R.Coefficient[k] * f, se=R.Coefficient_std_error[k] * f, dl=D[:, k] * f,
                                                         enr=R.Enrichment[k], enr_se=R.Enrichment_std_error[k], psnp=R["Prop._SNPs"][k])
    assert len(res) == 3, (t, q, list(res)); return res
def jse(dl): B = len(dl); return np.sqrt((B - 1) / B * ((dl - dl.mean()) ** 2).sum())
def jcov(dls): B = dls.shape[0]; z = dls - dls.mean(0); return (B - 1) / B * z.T @ z
PR = {(t, q): load(t, q) for t in TR for q in TH}
def con(t, q, a, b): p = PR[(t, q)]; return p[a]["tau"] - p[b]["tau"], p[a]["dl"] - p[b]["dl"]
rows = []
for t in TR:
    for q in TH:
        p = PR[(t, q)]; order = sorted(LAY, key=lambda m: -p[m]["tau"]); c, dl = con(t, q, order[0], order[1])
        rec = dict(trait=t, lineage=LINE[t], thr=q, top=order[0], second=order[1], top_vs_second=c, top_vs_second_p1=stats.norm.sf(c / jse(dl)))
        for m in LAY: rec.update({f"{m}_tau": p[m]["tau"], f"{m}_se": p[m]["se"], f"{m}_enr": p[m]["enr"], f"{m}_enr_se": p[m]["enr_se"], f"{m}_psnp": p[m]["psnp"]})
        for k, (a, b) in {"c1": ("CAGE", "DNASE"), "c2": ("CAGE", "H3K27ac"), "c3": ("DNASE", "H3K27ac")}.items(): c, dl = con(t, q, a, b); rec[k] = c; rec[f"{k}_se"] = jse(dl)
        rows.append(rec)
L = pd.DataFrame(rows); L.to_csv(EXP/"results"/f"u03_{STAGE}_layers.csv", index=False); LI = L.set_index(["trait", "thr"])
P("=== 역치별 3-층위 프로필 (전체 트랙, 단일 실행; τ* ± SE; * = 1위가 2위보다 단측 p < 0.05) ===")
for t in TR:
    P(f"{KO[t]} ({LINE[t]})")
    for q in TH:
        r = LI.loc[(t, q)]
        P(f"   {THL[q]:6s} " + "  ".join(f"{m} {r[f'{m}_tau']:+6.2f} ± {r[f'{m}_se']:.2f}" for m in LAY) + f"  → 1위 {r.top}{'*' if r.top_vs_second_p1 < .05 else ' '} "
          f"| 농축 " + " / ".join(f"{r[f'{m}_enr']:.1f}" for m in LAY) + f" (주석 비율 {r['CAGE_psnp'] * 100:.2f} %)")
if STAGE == "main":
    C1 = pd.read_csv(R3F/"v03_main_layers.csv").set_index("trait"); PRED = {"BMD": "DNASE", "NEUT": "CAGE", "TC": "CAGE"}   # 3f 3단계 예측 (프로토콜)
    def ivw(ts, q, k):
        w = np.array([1 / LI.loc[(t, q), f"{k}_se"] ** 2 for t in ts]); v = np.array([LI.loc[(t, q), k] for t in ts]); return (w * v).sum() / w.sum(), np.sqrt(1 / w.sum())
    S = {}
    for q in TH:
        d = {t: np.array([con(t, q, "CAGE", "DNASE")[0], con(t, q, "CAGE", "H3K27ac")[0]]) for t in MAIN}
        W = {t: np.linalg.inv(jcov(np.column_stack([con(t, q, "CAGE", "DNASE")[1], con(t, q, "CAGE", "H3K27ac")[1]]))) for t in MAIN}
        db = np.linalg.solve(sum(W.values()), sum(W[t] @ d[t] for t in MAIN)); Q = float(sum((d[t] - db) @ W[t] @ (d[t] - db) for t in MAIN)); pQ = stats.chi2.sf(Q, 2 * (len(MAIN) - 1))
        agree = [t for t in MAIN if LI.loc[(t, q), "top"] == C1.loc[t, "top"]]; h_ok = LI.loc[("HEIGHT", q), "top"] == "DNASE"
        pairs = []
        for lab, A, Bg, k in [("L-H2a AD 대 키 (c1)", ["AD"], ["HEIGHT"], "c1"), ("L-H2b 뇌 형질 4개 대 키 (c1)", ["EA", "SCZ", "BIP", "NEURO"], ["HEIGHT"], "c1"),
                              ("L-H2c 림프계 대 골수계 (c2)", ["IBD", "RA"], ["AD", "MONO"], "c2")]:
            (ma, sa), (mb, sb) = ivw(A, q, k), ivw(Bg, q, k); dd = ma - mb; se = np.sqrt(sa ** 2 + sb ** 2); pairs.append((lab, dd, se, 2 * stats.norm.sf(abs(dd / se))))
        (ma, sa), (mb, sb) = ivw(["NEUT", "TC"], q, "c1"), ivw(["BMD"], q, "c1"); vd = ma - mb; vse = np.sqrt(sa ** 2 + sb ** 2); vp = stats.norm.sf(vd / vse)
        v3 = {}
        for v in VAL:
            X1 = PRED[v]; ps = [stats.norm.sf(c / jse(dl)) for c, dl in (con(v, q, X1, o) for o in LAY if o != X1)]; v3[v] = (LI.loc[(v, q), "top"], all(p < .05 for p in ps), max(ps))
        ser = np.median([LI.loc[(t, q), "c1_se"] / LI.loc[(t, "L5"), "c1_se"] for t in MAIN])
        nconf = int(sum(LI.loc[(t, q), "top_vs_second_p1"] < .05 for t in MAIN))
        S[q] = dict(Q=Q, pQ=pQ, agree=agree, h_ok=h_ok, pairs=pairs, vd=vd, vse=vse, vp=vp, v3=v3, ser=ser, nconf=nconf)
    P("\n=== 역치별 검정 요약 ===")
    for q in TH:
        s = S[q]
        P(f"\n[{THL[q]}]" + (" ← 주 판정" if q == "L0p1" else " (기술)") + (" — 5 % 참조 = 3f T3all" if q == "L5" else ""))
        P(f"  G-H1 다변량 Q = {s['Q']:.1f}, 자유도 20, p = {s['pQ']:.3g} → {'통과' if s['pQ'] < .05 else '불통과'}")
        diff = ", ".join(f"{KO[t]} {C1.loc[t, 'top']}→{LI.loc[(t, q), 'top']}" for t in MAIN if t not in s["agree"])
        P(f"  G-H2 1위 층위가 3f 1단계와 같은 형질 {len(s['agree'])}/11, 키 1위 {LI.loc[('HEIGHT', q), 'top']} → {'충족' if len(s['agree']) >= 9 and s['h_ok'] else '불충족'}"
          + (f" (다른 형질: {diff})" if diff else ""))
        P(f"  1위가 2위보다 단측 p < 0.05인 형질 {s['nconf']}/11; c1 SE 비(이 역치 / 5 %) 중앙값 {s['ser']:.2f}")
        for lab, dd, se, p2 in s["pairs"]: P(f"  G-H3 {lab}: 차이 {dd:+.2f} ± {se:.2f}, p {p2:.3g} → {'유의' if p2 < .05 / 3 else '유의하지 않음'}")
        P(f"  G-H4 V-H1 (호중구 수·총콜레스테롤 c1 − 골밀도 c1): {s['vd']:+.2f} ± {s['vse']:.2f}, 단측 p {s['vp']:.3g} → {'통과' if s['vp'] < .05 else '불통과'}")
        P("  검증 형질 1위 층위 (예측): " + ", ".join(f"{KO[v]} {s['v3'][v][0]} (예측 {PRED[v]}, 엄격 일치 {'예' if s['v3'][v][1] else '아니오'})" for v in VAL))
    s = S["L0p1"]; g2 = len(s["agree"]) >= 9 and s["h_ok"]
    verdict = "재현" if (s["pQ"] < .05 and g2) else ("방향만 재현" if g2 else "재현 안 됨")
    P(f"\n=== 판정 (q = 0.1 %): {verdict} ===")
    P({"재현": "층위 지도는 Drusinsky & Pollard가 믿을 만하다고 한 상위 0.1 % 예측만으로도 나타난다.",
       "방향만 재현": f"1위 층위는 같지만 상위 0.1 % 주석에서는 형질 간 차이를 검정할 검정력이 부족하다 (c1 SE {s['ser']:.1f}배).",
       "재현 안 됨": "층위 지도는 역치에 의존한다. 바뀐 형질과 방향을 보고하고 본문 주장을 줄인다."}[verdict])
    P("\n=== 1위 층위 농축(Prop. h² / Prop. SNP) 역치별 (기술) ===")
    for t in MAIN + VAL:
        m = C1.loc[t, "top"] if t in MAIN else PRED[t]
        P(f"  {KO[t]:9s} {m:8s} " + "  ".join(f"{THL[q]} {LI.loc[(t, q), f'{m}_enr']:6.1f} ± {LI.loc[(t, q), f'{m}_enr_se']:.1f}" for q in TH))
    # 그림
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    plt.rcParams["font.family"] = "Malgun Gothic"; plt.rcParams["axes.unicode_minus"] = False
    COL = {"CAGE": "#1f4e8c", "DNASE": "#d9822b", "H3K27ac": "#4c9a6a"}; TXT = {"CAGE": "white", "DNASE": "#111", "H3K27ac": "white"}
    TCOL = {"L0p1": "#08306b", "L0p5": "#2b6aa8", "L1p0": "#6aa1d1", "L5": "#b3cde3"}
    ORD = MAIN + VAL; fig, ax = plt.subplots(1, 2, figsize=(13, 6.8), gridspec_kw=dict(width_ratios=[1, 1.25]))
    ny = len(ORD) + 1; ypos = {t: ny - 1 - i - (1 if t in VAL else 0) for i, t in enumerate(ORD)}
    for j, q in enumerate(TH):
        for t in ORD:
            r = LI.loc[(t, q)]; sig = r.top_vs_second_p1 < .05
            ax[0].add_patch(Rectangle((j + .04, ypos[t] - .44), .92, .88, facecolor=COL[r.top], alpha=1 if sig else .45, edgecolor="white", lw=1.5))
            ax[0].text(j + .5, ypos[t], f"{r.top.replace('DNASE', 'DNase')}{' *' if sig else ''}", ha="center", va="center", fontsize=8, color=TXT[r.top] if sig else "#111")
    ax[0].set_xlim(0, 4); ax[0].set_ylim(-.6, ny - .4); ax[0].set_xticks(np.arange(4) + .5); ax[0].set_xticklabels([f"상위 {THL[q]}\nQ p {S[q]['pQ']:.2g}" for q in TH], fontsize=8.5)
    ax[0].set_yticks([ypos[t] for t in ORD]); ax[0].set_yticklabels([f"{KO[t]} · {LINE[t]}" for t in ORD], fontsize=9); ax[0].tick_params(length=0)
    ax[0].axhline(ny - 1 - len(MAIN) + .0, color="#999", lw=.8, ls=":"); ax[0].text(4.05, ypos["TC"], "검증", fontsize=8, color="#555", va="center")
    ax[0].set_title("A. 역치별 1위 층위 (진한 칸·* = 2위보다 단측 p < 0.05)", loc="left", fontsize=11)
    for s_ in ["top", "right", "left", "bottom"]: ax[0].spines[s_].set_visible(False)
    for j, q in enumerate(TH):
        yy = [ypos[t] + (1.5 - j) * .18 for t in ORD]
        ax[1].errorbar([LI.loc[(t, q), "c1"] for t in ORD], yy, xerr=[1.96 * LI.loc[(t, q), "c1_se"] for t in ORD], fmt="o", ms=4.5, lw=1.2, color=TCOL[q], capsize=0, label=f"상위 {THL[q]}")
    ax[1].axvline(0, color="#555", lw=.8); ax[1].set_yticks([ypos[t] for t in ORD]); ax[1].set_yticklabels([KO[t] for t in ORD], fontsize=9)
    ax[1].axhline(ny - 1 - len(MAIN), color="#999", lw=.8, ls=":"); ax[1].set_ylim(-.6, ny - .4)
    ax[1].set_xlabel("CAGE - DNase (τ*, ± 95 % CI), 3-층위 동시 모형, 전체 트랙"); ax[1].set_title("B. 전사 시작 대 접근성, 역치별", loc="left", fontsize=11)
    ax[1].legend(fontsize=8, frameon=False, loc="lower right"); ax[1].spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(EXP/"figures"/"exp3g_threshold.png", dpi=200)
(EXP/"results"/f"u03_{STAGE}_tests.txt").write_text("\n".join(log), encoding="utf-8")
