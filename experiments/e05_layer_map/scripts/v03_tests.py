"""실험 3f v03: 시드 평균 층위 프로필과 검정 (PROTOCOL_exp3f_2026-10-09).
 시드 평균: 추정값은 시드 20개 평균, jackknife 삭제값도 같은 블록끼리 시드 평균 → SE = √((B−1)/B Σ (δ_b − δ̄)²).  대비 c1 = CAGE − DNASE, c2 = CAGE − H3K27ac, c3 = DNASE − H3K27ac (τ*).
 main: L-H1 다변량 Cochran Q (d = (c1, c2), 자유도 2(K−1)), L-H2a–c 쌍 비교(양측, α 0.05/3), 단변량 Q, 층위형 분류(확정: 시드 평균 1위 + 15/20 이상 1위 + 2위와 대비 단측 p < 0.05), T3all 비교.
 val: 예측 층위 = 유사 형질의 1단계 분류(HEIGHT→BMD, MONO→NEUT, LDL→TC). V-H1 형질 간 대비(단측), V-H2 예측 계통 M1 τ* > 0, V-H3 예측 층위 대 나머지 두 층위 모두 단측 p < 0.05.
 사용: python v03_tests.py dry | main | val
출력: results/v03_<단계>_layers.csv, v03_<단계>_tests.txt, figures/exp3f_layer_map.png (main)
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); C3 = Path(f"{DATA}/exp3c"); Y = Path(f"{DATA}/exp3d"); F = Path(f"{DATA}/exp3f")
SL = F/"test"/"sldsc" if STAGE == "dry" else F/"sldsc"; (EXP/"figures").mkdir(exist_ok=True)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LAY = ["CAGE", "DNASE", "H3K27ac"]
SDF = F/("annot_sd_common.csv" if STAGE != "dry" or (F/"annot_sd_common.csv").exists() else "test_annot_sd_common.csv")
SD = pd.concat([pd.read_csv(SDF, index_col=0).sd_common, pd.read_csv(Y/"annot_sd_common.csv", index_col=0).sd_common, pd.read_csv(C3/"annot_sd_common.csv", index_col=0).sd_common])
B0 = X/"GRCh38"/"baselineLD_v2.2"; M_TOT = sum(float(open(B0/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
if STAGE == "main":
    TR = ["AD", "MONO", "EA", "SCZ", "BIP", "NEURO", "BMI", "IBD", "RA", "LDL", "HEIGHT"]
    LINE = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue", "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal"}
elif STAGE == "val": TR = ["BMD", "NEUT", "TC"]; LINE = {"BMD": "stromal", "NEUT": "myeloid", "TC": "liver"}
else: TR = ["FAKE"]; LINE = {"FAKE": "myeloid"}
SEEDS = [f"{s:02d}" for s in range(1, 3 if STAGE == "dry" else 21)]
KO = {"AD": "AD", "MONO": "단핵구 수", "EA": "교육 연수", "SCZ": "조현병", "BIP": "양극성장애", "NEURO": "신경증", "BMI": "BMI", "IBD": "IBD", "RA": "RA", "LDL": "LDL", "HEIGHT": "키",
      "BMD": "골밀도", "NEUT": "호중구 수", "TC": "총콜레스테롤", "FAKE": "가짜"}
def load(t, model):
    out = SL/f"{t}__{model}"; R = pd.read_csv(f"{out}.results", sep="\t"); R["name"] = R.Category.str.replace(r"L2_\d+$", "", regex=True)
    h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log", encoding="utf-8").read()).group(1))
    D = pd.read_csv(f"{out}.part_delete", sep=r"\s+", header=None).to_numpy(); res = {}
    for n in R.name:
        if n.startswith(("L5r", "L5_")):
            k = int(R.index[R.name == n][0]); f = SD[n] * M_TOT / h2; res[n] = (R.Coefficient[k] * f, D[:, k] * f, R.Coefficient_std_error[k] * f)
    return res, h2
def jse(dl): B = len(dl); return np.sqrt((B - 1) / B * ((dl - dl.mean()) ** 2).sum())
def jcov(dls): B = dls.shape[0]; z = dls - dls.mean(0); return (B - 1) / B * z.T @ z
def lay_of(n): return [m for m in LAY if n.endswith("_" + m)][0]
PROF = {}   # 형질 → {"est": {층위: 시드 평균 τ*}, "del": {층위: 시드 평균 삭제값}, "seed_est": [{층위: τ*}], "all": {층위: τ*}}
for t in TR:
    g = LINE[t]; se_, sd_ = [], []
    for s in SEEDS:
        r, _ = load(t, f"T3_{s}"); assert len(r) == 3, (t, s, list(r)); se_.append({lay_of(n): v[0] for n, v in r.items()}); sd_.append({lay_of(n): v[1] for n, v in r.items()})
    ra, _ = load(t, "T3all"); assert len(ra) == 3
    PROF[t] = dict(est={m: np.mean([e[m] for e in se_]) for m in LAY}, dl={m: np.mean([d[m] for d in sd_], axis=0) for m in LAY}, seed_est=se_, all={lay_of(n): v[0] for n, v in ra.items()},
                   all_se={lay_of(n): v[2] for n, v in ra.items()})
def contrast(t, a, b): p = PROF[t]; return p["est"][a] - p["est"][b], p["dl"][a] - p["dl"][b]
rows = []
for t in TR:
    p = PROF[t]; top = pd.Series([max(e, key=e.get) for e in p["seed_est"]]).value_counts()
    order = sorted(LAY, key=lambda m: -p["est"][m]); X1, X2 = order[0], order[1]; cx, dx = contrast(t, X1, X2); px = stats.norm.sf(cx / jse(dx))
    ntop = int(top.get(X1, 0)); conf = (ntop >= .75 * len(SEEDS)) and (px < .05)
    rec = dict(trait=t, lineage=LINE[t], top=X1, n_top=ntop, second=X2, top_vs_second=cx, top_vs_second_p1=px, cls=(f"확정 {X1}형" if conf else f"미정 (최빈 {top.index[0]} {int(top.iloc[0])}/{len(SEEDS)})"),
               top_full=max(p["all"], key=p["all"].get))
    for m in LAY: rec[f"{m}_tau"] = p["est"][m]; rec[f"{m}_se"] = jse(p["dl"][m]); rec[f"{m}_ntop"] = int(top.get(m, 0)); rec[f"{m}_full"] = p["all"][m]; rec[f"{m}_full_se"] = p["all_se"][m]
    for k, (a, b) in {"c1": ("CAGE", "DNASE"), "c2": ("CAGE", "H3K27ac"), "c3": ("DNASE", "H3K27ac")}.items(): c, d = contrast(t, a, b); rec[k] = c; rec[f"{k}_se"] = jse(d)
    rows.append(rec)
L = pd.DataFrame(rows); L.to_csv(EXP/"results"/f"v03_{STAGE}_layers.csv", index=False)
P(f"=== 시드 평균 3-층위 프로필 (T3, 시드 {len(SEEDS)}개; 괄호 = 시드 중 1위 횟수) ===")
P("형질        계통        CAGE                DNase               H3K27ac             분류                         전체 트랙 1위")
for r in L.itertuples():
    P(f"{KO[r.trait]:9s}  {r.lineage:10s}  " + "  ".join(f"{getattr(r, m + '_tau'):+5.2f} ± {getattr(r, m + '_se'):.2f} ({getattr(r, m + '_ntop'):2d})" for m in LAY) + f"   {r.cls:26s} {r.top_full}")
P("\n대비 (시드 평균 τ*, ± SE): c1 = CAGE − DNase, c2 = CAGE − H3K27ac, c3 = DNase − H3K27ac")
for r in L.itertuples(): P(f"  {KO[r.trait]:9s} c1 {r.c1:+5.2f} ± {r.c1_se:.2f}   c2 {r.c2:+5.2f} ± {r.c2_se:.2f}   c3 {r.c3:+5.2f} ± {r.c3_se:.2f}")
def ivw(ts, k):
    w = np.array([1 / L.set_index("trait").loc[t, f"{k}_se"] ** 2 for t in ts]); v = np.array([L.set_index("trait").loc[t, k] for t in ts]); return (w * v).sum() / w.sum(), np.sqrt(1 / w.sum())
if STAGE == "main":
    d = {t: np.array([contrast(t, "CAGE", "DNASE")[0], contrast(t, "CAGE", "H3K27ac")[0]]) for t in TR}
    S = {t: jcov(np.column_stack([contrast(t, "CAGE", "DNASE")[1], contrast(t, "CAGE", "H3K27ac")[1]])) for t in TR}
    W = {t: np.linalg.inv(S[t]) for t in TR}; dbar = np.linalg.solve(sum(W.values()), sum(W[t] @ d[t] for t in TR))
    Q = float(sum((d[t] - dbar) @ W[t] @ (d[t] - dbar) for t in TR)); df = 2 * (len(TR) - 1); pQ = stats.chi2.sf(Q, df)
    P(f"\n=== L-H1 (주): 다변량 Cochran Q = {Q:.1f}, 자유도 {df}, p = {pQ:.3g} → {'통과: 유전율을 운반하는 조절 층위는 형질마다 다르다' if pQ < .05 else '불통과: 트랙 균형을 맞추면 형질 간 층위 차이는 검출되지 않는다'}")
    P(f"  가중 평균 프로필 d̄ = (CAGE − DNase {dbar[0]:+.2f}, CAGE − H3K27ac {dbar[1]:+.2f})")
    P("단변량 Q (자유도 10, 기술):")
    for k in ["c1", "c2", "c3"]:
        w = 1 / L[f"{k}_se"] ** 2; m = (w * L[k]).sum() / w.sum(); q = float((w * (L[k] - m) ** 2).sum()); P(f"  {k}: Q {q:.1f}, p {stats.chi2.sf(q, len(TR) - 1):.3g} (가중 평균 {m:+.2f})")
    P("\n=== L-H2 쌍 비교 (양측, α 0.05/3 = 0.0167) ===")
    for lab, A, Bg, k in [("L-H2a AD 대 키 (c1 CAGE − DNase)", ["AD"], ["HEIGHT"], "c1"), ("L-H2b 뇌 형질 4개 대 키 (c1)", ["EA", "SCZ", "BIP", "NEURO"], ["HEIGHT"], "c1"),
                          ("L-H2c 림프계(IBD·RA) 대 골수계(AD·단핵구 수) (c2 CAGE − H3K27ac)", ["IBD", "RA"], ["AD", "MONO"], "c2")]:
        (ma, sa), (mb, sb) = ivw(A, k), ivw(Bg, k); dd = ma - mb; se = np.sqrt(sa ** 2 + sb ** 2); p2 = 2 * stats.norm.sf(abs(dd / se))
        P(f"  {lab}: {ma:+.2f} ± {sa:.2f} 대 {mb:+.2f} ± {sb:.2f} → 차이 {dd:+.2f} ± {se:.2f}, p {p2:.3g} → {'유의' if p2 < .05 / 3 else '유의하지 않음'}")
    P(f"\n층위형 분류: " + ", ".join(f"{KO[r.trait]} {r.cls}" for r in L.itertuples()))
    P(f"전체 트랙 판(T3all)과 1위 층위가 다른 형질: {', '.join(f'{KO[r.trait]} (전체 {r.top_full} → 균형 {r.top})' for r in L.itertuples() if r.top_full != r.top) or '없음'}")
    # 그림
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "Malgun Gothic"; plt.rcParams["axes.unicode_minus"] = False
    COL = {"CAGE": "#1f4e8c", "DNASE": "#d9822b", "H3K27ac": "#4c9a6a"}; LAB = {"CAGE": "CAGE (전사 시작)", "DNASE": "DNase (접근성)", "H3K27ac": "H3K27ac (활성 크로마틴)"}
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.6), gridspec_kw=dict(width_ratios=[1.25, 1])); y = np.arange(len(L))[::-1]
    for j, m in enumerate(LAY):
        ax[0].errorbar(L[f"{m}_tau"], y + (1 - j) * .22, xerr=1.96 * L[f"{m}_se"], fmt="o", color=COL[m], ms=5.5, lw=1.4, capsize=0, label=LAB[m])
    ax[0].axvline(0, color="#555", lw=.8); ax[0].set_yticks(y); ax[0].set_yticklabels([f"{KO[t]} · {g}" for t, g in zip(L.trait, L.lineage)])
    for yy, c in zip(y, L.cls): ax[0].text(1.01, yy, c.split(" (")[0], transform=ax[0].get_yaxis_transform(), va="center", fontsize=8, color="#333")
    ax[0].set_xlabel("시드 평균 τ* (± 95 % CI), 3-층위 동시 모형"); ax[0].set_title("A. 형질별 층위 프로필 (트랙 균형, 시드 20개)", loc="left", fontsize=11); ax[0].legend(fontsize=8, frameon=False, loc="lower right")
    o = L.sort_values("c1"); yo = np.arange(len(o))
    ax[1].errorbar(o.c1, yo, xerr=1.96 * o.c1_se, fmt="o", color="#333", ms=5.5, lw=1.4, capsize=0)
    ax[1].axvline(0, color="#555", lw=.8); ax[1].set_yticks(yo); ax[1].set_yticklabels([KO[t] for t in o.trait])
    ax[1].set_xlabel("CAGE - DNase (τ*, ± 95 % CI)"); ax[1].set_title(f"B. 전사 시작 대 접근성 (다변량 Q p = {pQ:.2g})", loc="left", fontsize=11)
    for a in ax: a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(EXP/"figures"/"exp3f_layer_map.png", dpi=200)
if STAGE == "val":
    C1 = pd.read_csv(EXP/"results"/"v03_main_layers.csv").set_index("trait"); ANA = {"BMD": "HEIGHT", "NEUT": "MONO", "TC": "LDL"}
    PRED = {v: (C1.loc[a, "top"] if C1.loc[a, "cls"].startswith("확정") else None) for v, a in ANA.items()}
    P("\n=== 예측 (프로토콜 규칙: 유사 형질의 1단계 분류) ===")
    for v, a in ANA.items(): P(f"  {KO[v]} ← {KO[a]} ({C1.loc[a, 'cls']}) → 예측 층위 {PRED[v] or '없음(기술로만 보고)'}")
    P("\n=== V-H2 예측 계통 (M1 τ* > 0, 단측 α 0.05) ===")
    LIN12 = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
    for v in TR:
        m1 = {}
        for h in LIN12:
            r, _ = load(v, f"M1_{h}"); (n, (e, dl, se)), = r.items(); m1[h] = (e, se, stats.norm.sf(e / se))
        rk = sorted(LIN12, key=lambda h: -m1[h][0]).index(LINE[v]) + 1; e, se, p1 = m1[LINE[v]]
        P(f"  {KO[v]} {LINE[v]}: τ* {e:+.2f} ± {se:.2f}, p {p1:.3g} → {'통과' if p1 < .05 else '불통과'} | 12 계통 중 {rk}위 (1위 {max(m1, key=lambda h: m1[h][0])} {max(x[0] for x in m1.values()):+.2f})")
    P("\n=== V-H3 예측 층위형 (예측 층위 대 나머지 두 층위, 모두 단측 p < 0.05) ===")
    for v in TR:
        X1 = PRED[v]
        if X1 is None: P(f"  {KO[v]}: 예측 없음 — 기술: 1위 {L.set_index('trait').loc[v, 'top']}"); continue
        ps = []
        for o2 in [m for m in LAY if m != X1]: c, dl = contrast(v, X1, o2); ps.append((o2, c, stats.norm.sf(c / jse(dl))))
        ok = all(p < .05 for _, _, p in ps)
        P(f"  {KO[v]} 예측 {X1}: " + ", ".join(f"{X1} − {o2} {c:+.2f} (p {p:.3g})" for o2, c, p in ps) + f" → {'일치' if ok else '불일치'} | 시드 중 {X1} 1위 {int(L.set_index('trait').loc[v, X1 + '_ntop'])}/20")
    types = sorted(set(x for x in PRED.values() if x))
    if len(types) >= 2:
        a, b = ("CAGE", "DNASE") if {"CAGE", "DNASE"} <= set(types) else (types[0], types[1])
        GA = [v for v in TR if PRED[v] == a]; GB = [v for v in TR if PRED[v] == b]
        def grp(ts):
            vals = [contrast(v, a, b) for v in ts]; w = np.array([1 / jse(dl) ** 2 for _, dl in vals]); return (w * np.array([c for c, _ in vals])).sum() / w.sum(), np.sqrt(1 / w.sum())
        (ma, sa), (mb, sb) = grp(GA), grp(GB); dd = ma - mb; se = np.sqrt(sa ** 2 + sb ** 2); p1 = stats.norm.sf(dd / se)
        P(f"\n=== V-H1 (주): ({a} − {b}) 대비, {a}형 예측 {[KO[v] for v in GA]} {ma:+.2f} ± {sa:.2f} 대 {b}형 예측 {[KO[v] for v in GB]} {mb:+.2f} ± {sb:.2f}"
          f" → 차이 {dd:+.2f} ± {se:.2f}, 단측 p {p1:.3g} → {'통과' if p1 < .05 else '불통과'} (α 0.05)")
    else: P("\n=== V-H1: 예측된 층위형이 두 가지 미만이라 판정하지 않는다(기술) ===")
(EXP/"results"/f"v03_{STAGE}_tests.txt").write_text("\n".join(log), encoding="utf-8")
