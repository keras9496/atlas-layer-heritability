"""원고 v2 그림 (영문). 기존 결과 CSV만 읽는다(새 계산 없음).
 Fig 1 설계 / Fig 2 계통 지도(3c M1) / Fig 3 층위 지도(3f 1단계) / Fig 4 사전 예측 검증(3f 3단계) / Fig 5 자기 층위 효과의 실측·활성도 조건화(3e w05, 3l l05)
 Supp Fig S1 계통 특이 지도(3c M2) / S2 단일 주석 계통 × 층위 지도(3d D1) / S3 역치별 1위 층위(3g) / S4 CAGE − DNase 포함 관계 역치(3m M2)·연속·구간(3o) / S5 예측 대 실측·eQTL 크기 비교(3e, v1 그림 5b–d)
출력: figures/fig1.png … fig5.png, figS1.png … figS5.png, numbers.json (원고에 쓴 수치 대조용)
"""
import json
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch
P = Path(__file__).resolve().parent; OUT = P/"figures"; OUT.mkdir(exist_ok=True)
R = Path(__file__).resolve().parents[1]/"experiments"
plt.rcParams.update({"font.family": "Arial", "font.size": 8.5, "axes.titlesize": 9.5, "axes.labelsize": 8.5, "axes.unicode_minus": False, "svg.fonttype": "none"})
COL = {"CAGE": "#1f4e8c", "DNASE": "#d9822b", "H3K27ac": "#4c9a6a"}
LAYLAB = {"CAGE": "CAGE (transcription initiation)", "DNASE": "DNase (accessibility)", "H3K27ac": "H3K27ac (active chromatin)"}
SHORT = {"CAGE": "CAGE", "DNASE": "DNase", "H3K27ac": "H3K27ac"}
LAY = ["CAGE", "DNASE", "H3K27ac"]
TN = {"AD": "Alzheimer's disease", "MONO": "Monocyte count", "EA": "Educational attainment", "SCZ": "Schizophrenia", "BIP": "Bipolar disorder", "NEURO": "Neuroticism",
      "BMI": "BMI", "IBD": "IBD", "RA": "Rheumatoid arthritis", "LDL": "LDL cholesterol", "HEIGHT": "Height", "BMD": "Bone mineral density", "NEUT": "Neutrophil count", "TC": "Total cholesterol"}
LN = {"CNS_tissue": "Brain tissue", "neural_cells": "Neural cells", "myeloid": "Myeloid", "lymphoid": "Lymphoid", "liver": "Liver", "GI": "Gastrointestinal",
      "cardiovascular": "Cardiovascular", "skeletal_muscle": "Skeletal muscle", "adipose": "Adipose", "stromal": "Stromal/mesenchymal",
      "reproductive_endocrine": "Reproductive/endocrine", "epithelial_barrier": "Epithelial/barrier"}
LIN12 = list(LN)
ORD = ["EA", "SCZ", "BIP", "NEURO", "BMI", "AD", "MONO", "IBD", "RA", "LDL", "HEIGHT"]   # 뇌 형질 → AD → 혈액·면역 → LDL → 키
INK, MUTED, GRID = "#222222", "#666666", "#d0d0d0"
NUM = {}
def clean(ax): ax.spines[["top", "right"]].set_visible(False); ax.tick_params(length=2.5, color=MUTED)
def panel(ax, s, x=-0.02, y=1.02): ax.text(x, y, s, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="right")
DIV = LinearSegmentedColormap.from_list("div", ["#2b6aa8", "#9ec1e0", "#efefef", "#f2b38f", "#c4501f"])
def fmtp(p):
    if p >= 0.001: return f"{p:.3f}".rstrip("0").rstrip(".") if p >= 0.01 else f"{p:.3f}"
    e = int(np.floor(np.log10(p))); m = p / 10 ** e; return f"{m:.0f} × 10$^{{{e}}}$" if round(m, 1) >= 9.95 or abs(m - round(m)) < .05 else f"{m:.1f} × 10$^{{{e}}}$"

# ------------------------------------------------------------------ 자료
S3C = pd.read_csv(R/"e02_lineage_map/results/z04_sldsc.csv")
L3F = pd.read_csv(R/"e05_layer_map/results/v03_main_layers.csv").set_index("trait")
V3F = pd.read_csv(R/"e05_layer_map/results/v03_val_layers.csv").set_index("trait")
OWN = pd.read_csv(R/"e04_measured_annotations/results/w05_own_layer.csv").set_index("trait")
SM = pd.read_csv(R/"e04_measured_annotations/results/w04_Sm.csv")
D1 = pd.read_csv(R/"e03_single_layer_map/results/y04_B.csv")

# ------------------------------------------------------------------ Fig 1 설계
FH = 2.66   # 그림 1 높이(인치); 상자 안 글 위치는 인치로 맞춘다
def box(ax, x, y, w, h, title, body, fc="#f4f6f9", ec="#9aa7b8"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.005,rounding_size=0.014", fc=fc, ec=ec, lw=1))
    ax.text(x + w / 2, y + h - .11 / FH, title, ha="center", va="top", fontsize=8.7, fontweight="bold", color=INK)
    ax.text(x + w / 2, y + h - .33 / FH, body, ha="center", va="top", fontsize=7.2, color=INK, linespacing=1.35)
def arrow(ax, x0, y0, x1, y1): ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="-|>", mutation_scale=9, color="#5a6b80", lw=1.1, shrinkA=0, shrinkB=0))
fig = plt.figure(figsize=(7.2, FH)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
L0, R0, GAP = .045, .985, .05; BW = (R0 - L0 - 3 * GAP) / 4; XS = [L0 + i * (BW + GAP) for i in range(4)]; XC = [x + BW / 2 for x in XS]   # 위아래 같은 4열 격자
HA, HB = 1.2 / FH, .86 / FH; YA = 1 - .05 / FH - HA; YB = .09 / FH; LINE = "#5a6b80"; PAD = .005; BUS = (YA - PAD + YB + HB + PAD) / 2   # 상자 높이 = 글 분량
ax.text(.008, YA + HA + PAD, "a", fontsize=11, fontweight="bold", va="top"); ax.text(.008, YB + HB + PAD, "b", fontsize=11, fontweight="bold", va="top")
TOP = [("AlphaGenome Atlas", "Allelic effect\n|log$_2$(ΣALT+1)−log$_2$(ΣREF+1)|\nfor every common SNP\n\n5.83 million SNPs (MAF ≥ 5%)\n1,312 tracks, 433 biosamples\n(cell lines excluded)"),
       ("Annotations", "Tracks grouped into\n12 lineages × 6 layers\n\nper SNP: maximum effect\nover the tracks of a cell\n→ top 5% of SNPs"),
       ("S-LDSC", "Conditional on baselineLD v2.2\n(97 annotations, including\nfine-mapped eQTLs,\nFANTOM5 enhancers, TSS)\n\nτ* = per-SD heritability effect"),
       ("Lineage × layer map", "11 traits\n+ 3 prediction traits\n\nwhich lineage,\nand which regulatory layer\nwithin that lineage")]
for i, (t, b) in enumerate(TOP): box(ax, XS[i], YA, BW, HA, t, b, **(dict(fc="#eef3f9", ec="#5a7fa8") if i == 3 else {}))
for i in range(3): arrow(ax, XS[i] + BW + PAD + .008, YA + HA / 2, XS[i + 1] - PAD - .008, YA + HA / 2)   # 상자 사이 간격의 가운데
QS = [("Which lineage?", "Lineage annotation\nper lineage (Fig. 2)"),
      ("Which layer?", "CAGE, DNase, H3K27ac\nfitted jointly; up to\nthree tracks per layer,\n20 resamples (Fig. 3)"),
      ("Does it predict?", "Layer of 3 new traits\npredicted from analogous\ntraits before analysis\n(Fig. 4)"),
      ("Beyond measurement?", "Measured FANTOM5 and\nENCODE annotations;\neQTLs; predicted\nactivity (Fig. 5)")]
ax.plot([XC[3], XC[3]], [YA - PAD, BUS], color=LINE, lw=1.1, solid_capstyle="butt"); ax.plot([XC[0], XC[3]], [BUS, BUS], color=LINE, lw=1.1, solid_capstyle="butt")   # 지도 → 네 질문
for i, (t, b) in enumerate(QS): arrow(ax, XC[i], BUS, XC[i], YB + HB + PAD); box(ax, XS[i], YB, BW, HB, t, b, fc="#ffffff", ec="#b9c3cf")
fig.savefig(OUT/"fig1.png", dpi=300); plt.close(fig)

# ------------------------------------------------------------------ Fig 2 계통 지도 (M1)
def lineage_map(model, fname, title):
    T = S3C[S3C.model == model].pivot(index="trait", columns="lineage", values="tau_star").reindex(index=ORD, columns=LIN12)
    Q = S3C[S3C.model == model].pivot(index="trait", columns="lineage", values="q_bh").reindex(index=ORD, columns=LIN12)
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    im = ax.imshow(T.values, cmap=DIV, norm=TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=4), aspect="auto")
    for i in range(T.shape[0]):
        for j in range(T.shape[1]):
            v, q = T.iat[i, j], Q.iat[i, j]; sig = q < .05
            ax.text(j, i, ("0.0" if abs(v) < .05 else f"{v:+.1f}") + ("*" if sig else ""), ha="center", va="center", fontsize=6.8 if sig else 6.2, color=INK if abs(v) < 3.2 else "white", fontweight="bold" if sig else "normal")
    ax.set_xticks(range(12)); ax.set_xticklabels([LN[l] for l in LIN12], rotation=40, ha="right"); ax.set_yticks(range(len(ORD))); ax.set_yticklabels([TN[t] for t in ORD])
    for y in [4.5, 5.5, 8.5, 9.5]: ax.axhline(y, color="white", lw=2.2)
    ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]
    cb = fig.colorbar(im, ax=ax, fraction=.03, pad=.015); cb.set_label("τ* (per-SD heritability effect)"); cb.outline.set_visible(False)
    fig.tight_layout(); fig.savefig(OUT/fname, dpi=300); plt.close(fig)
    return T, Q
T2, Q2 = lineage_map("M1", "fig2.png", "Heritability enrichment of predicted effects by lineage (baselineLD + one lineage annotation)")
lineage_map("M2", "figS1.png", "Lineage-specific signal (baselineLD + all-lineage annotation + one lineage annotation)")
NUM["fig2_sig_cells"] = int((Q2 < .05).values.sum())

# ------------------------------------------------------------------ Fig 3 층위 지도 (3f 1단계)
CLS = {"확정 CAGE형": "CAGE type", "확정 DNASE형": "DNase type"}
def cls_en(c): return CLS.get(c, "Unresolved (CAGE modal)" if "CAGE" in c else c)
# 더 큰 GWAS(3n, 재추출 평균)는 b·c에 빈 기호로 겹친다. (2026-10-29: 3o의 방향 판정 열은 뺐다 — 역치 민감도 내용을 본문에서 제거)
N3L = pd.read_csv(R/"e13_larger_gwas/results/n03_layers.csv").query("gwas == 'new'").set_index("trait")
R3 = pd.read_csv(R/"e17_larger_ad_gwas/results/r03_summary.csv").set_index("test")   # 더 큰 AD GWAS(3r): c1·c2만 쓴다
N3L.loc["AD", ["c1", "c1_se", "c2", "c2_se"]] = [R3.loc["R-H4 c1", "new"], R3.loc["R-H4 c1", "new_se"], R3.loc["R-H4 c2", "new"], R3.loc["R-H4 c2", "new_se"]]
NEWT = list(N3L.index)
fig = plt.figure(figsize=(7.2, 4.8)); a0 = fig.add_axes([.205, .1, .325, .755]); a1 = fig.add_axes([.565, .1, .185, .755], sharey=a0); a2 = fig.add_axes([.785, .1, .185, .755], sharey=a0)
y = np.arange(len(ORD))[::-1]
for j, m in enumerate(LAY):
    a0.errorbar([L3F.loc[t, f"{m}_tau"] for t in ORD], y + (1 - j) * .24, xerr=[1.96 * L3F.loc[t, f"{m}_se"] for t in ORD], fmt="o", ms=3.8, lw=1.1, color=COL[m], capsize=0, label=SHORT[m])
a0.axvline(0, color=MUTED, lw=.7); a0.set_yticks(y); a0.set_yticklabels([TN[t] for t in ORD]); a0.set_ylim(-.6, len(ORD) - .4)
for a in (a0, a1, a2):
    for yy in [5.5, 4.5, 1.5, .5]: a.axhline(yy, color=GRID, lw=.6, ls=":")
a0.set_xlabel("τ* (seed average, ± 95% CI)"); a0.legend(fontsize=7.5, frameon=False, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, handletextpad=.3, columnspacing=1.2)
clean(a0); panel(a0, "a", x=-.36, y=1.0)
yo = {t: yy + (.14 if t in NEWT else 0) for yy, t in zip(y, ORD)}                       # 원래 GWAS (채운 기호)
for a, c, lab, p in [(a1, "c1", "CAGE − DNase", "b"), (a2, "c2", "CAGE − H3K27ac", "c")]:
    a.errorbar([L3F.loc[t, c] for t in ORD], [yo[t] for t in ORD], xerr=[1.96 * L3F.loc[t, f"{c}_se"] for t in ORD], fmt="o", ms=4, lw=1.1, color=INK, capsize=0, label="Main GWAS")
    a.errorbar([N3L.loc[t, c] for t in NEWT], [yo[t] - .28 for t in NEWT], xerr=[1.96 * N3L.loc[t, f"{c}_se"] for t in NEWT], fmt="o", ms=4, lw=1.1, color=INK, mfc="white", mew=1.1, capsize=0, label="Larger GWAS")
    a.axvline(0, color=MUTED, lw=.7); a.tick_params(labelleft=False); a.set_xlabel(lab); clean(a); panel(a, p, x=-.08, y=1.0)
h, l = a1.get_legend_handles_labels(); fig.legend(h, l, fontsize=7, frameon=False, loc="center left", bbox_to_anchor=(.545, .895), ncol=2, handletextpad=.2, columnspacing=.8)
fig.text(.565, .972, "Main map: Q = 74.3 (20 df), P = 3.6 × 10$^{-8}$", fontsize=7, va="top", color=INK); fig.text(.565, .942, "Larger GWAS for six traits: Q = 64.8 (20 df), P = 1.3 × 10$^{-6}$", fontsize=7, va="top", color=INK)   # 3j·3n·3r 공동 jackknife, 재추출 평균
fig.savefig(OUT/"fig3.png", dpi=300); plt.close(fig)
NUM["fig3_new_c1"] = {t: round(float(N3L.loc[t, "c1"]), 2) for t in NEWT}

# ------------------------------------------------------------------ Fig 4 사전 예측 검증
PAIRS = [("HEIGHT", "BMD", "Height → BMD", "stromal", "DNase"), ("LDL", "TC", "LDL → total cholesterol", "liver", "CAGE"), ("MONO", "NEUT", "Monocytes → neutrophils", "myeloid", "CAGE")]
ABB = {"HEIGHT": "Height", "BMD": "BMD", "LDL": "LDL", "TC": "Total chol.", "MONO": "Monocytes", "NEUT": "Neutrophils"}
fig = plt.figure(figsize=(7.2, 4.6)); gs = fig.add_gridspec(2, 3, height_ratios=[1, .85], hspace=.75, wspace=.35)
for k, (a_, v_, title, lin, pred) in enumerate(PAIRS):
    ax = fig.add_subplot(gs[0, k])
    for i, (src, t) in enumerate([(L3F, a_), (V3F, v_)]):
        for j, m in enumerate(LAY):
            x = i + (j - 1) * .24; r = src.loc[t]
            ax.bar(x, r[f"{m}_tau"], width=.21, color=COL[m], alpha=1 if i == 1 else .5, label=SHORT[m] if (k == 0 and i == 1) else None)
            ax.errorbar(x, r[f"{m}_tau"], yerr=1.96 * r[f"{m}_se"], color=INK, lw=.8, capsize=0)
    ax.axhline(0, color=MUTED, lw=.7); ax.set_xticks([0, 1]); ax.set_xticklabels([f"{ABB[a_]}\n(analog)", f"{ABB[v_]}\n(new)"], fontsize=7.5)
    ok = V3F.loc[v_, "cls"].startswith("확정")
    ax.set_title(f"{title} ({lin})\npredicted {pred}; " + (f"observed {pred}" if ok else "CAGE vs H3K27ac unresolved"), loc="left", fontsize=7.5); clean(ax); panel(ax, "abc"[k], x=-.1, y=1.13)
    if k == 0: ax.set_ylabel("τ* (± 95% CI)")
fig.legend(loc="lower center", ncol=3, frameon=False, fontsize=7.5, bbox_to_anchor=(.5, .985))
ax = fig.add_subplot(gs[1, :]); rows = [("NEUT", "CAGE"), ("TC", "CAGE"), ("BMD", "DNASE")]
for i, (t, pr) in enumerate(rows): ax.errorbar(V3F.loc[t, "c1"], i, xerr=1.96 * V3F.loc[t, "c1_se"], fmt="o", ms=5, lw=1.2, color=COL[pr], capsize=0)
def ivw(ts):
    w = np.array([1 / V3F.loc[t, "c1_se"] ** 2 for t in ts]); v = np.array([V3F.loc[t, "c1"] for t in ts]); return (w * v).sum() / w.sum(), np.sqrt(1 / w.sum())
(mc, sc), (md, sd) = ivw(["NEUT", "TC"]), ivw(["BMD"])
ax.errorbar(mc, 3.2, xerr=1.96 * sc, fmt="D", ms=6, color=COL["CAGE"], lw=1.4, capsize=0); ax.errorbar(md, 3.9, xerr=1.96 * sd, fmt="D", ms=6, color=COL["DNASE"], lw=1.4, capsize=0)
ax.set_yticks([0, 1, 2, 3.2, 3.9]); ax.set_yticklabels([TN["NEUT"], TN["TC"], TN["BMD"], "predicted CAGE type (pooled)", "predicted DNase type"]); ax.invert_yaxis()
ax.axvline(0, color=MUTED, lw=.7); ax.set_xlabel("CAGE − DNase (τ*, ± 95% CI)"); clean(ax); panel(ax, "d", x=-.31, y=1.0)
dd = mc - md; dse = np.sqrt(sc ** 2 + sd ** 2)
from scipy import stats
p1 = stats.norm.sf(dd / dse); NUM["fig4_contrast"] = dict(diff=dd, se=dse, p1=p1)
J3 = pd.read_csv(R/"e09_joint_jackknife/results/j03_main_summary.csv"); vj = J3[J3.test.str.startswith("검증 대비")].iloc[0]   # 3j 공동 jackknife
NUM["fig4_contrast_joint"] = dict(diff=float(vj.est), se=float(vj.stat_c), p1=float(vj.p_c))
ax.text(.99, .05, f"difference {vj.est:+.2f} ± {vj.stat_c:.2f}\none-sided P = {fmtp(vj.p_c)}", transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5)
fig.savefig(OUT/"fig4.png", dpi=300, bbox_inches="tight"); plt.close(fig)

# ------------------------------------------------------------------ Fig 5 자기 층위 효과: 실측 조건화(3e w05)와 활성도 조건화(3l l05)
ORD5 = ORD
V5 = pd.read_csv(R/"e11_activity/results/l05_values.csv"); Q5 = pd.read_csv(R/"e11_activity/results/l04_main_summary.csv")
def v5(t, model, kind, m): r = V5[(V5.trait == t) & (V5.model == model) & (V5.kind == kind) & (V5.layer == m)]; assert len(r) == 1; return float(r.tau.iloc[0]), float(r.se.iloc[0])
yy = np.arange(len(ORD5))[::-1]
def dumb(ax, a, b, ea, eb, la, lb, ca, cb, xl):
    for y_, x0, x1 in zip(yy, b, a): ax.plot([x0, x1], [y_, y_], color=GRID, lw=1.6, zorder=1)
    ax.errorbar(b, yy, xerr=1.96 * np.array(eb), fmt="o", ms=4, color=cb, lw=.8, capsize=0, label=lb, zorder=2)
    ax.errorbar(a, yy, xerr=1.96 * np.array(ea), fmt="o", ms=4.5, color=ca, lw=.8, capsize=0, label=la, zorder=3)
    ax.axvline(0, color=MUTED, lw=.7); ax.set_xlabel(xl); clean(ax); ax.set_ylim(-.7, len(ORD5) - .3)
fig = plt.figure(figsize=(7.2, 4.4)); B0, HH = .2, .68
axs = [fig.add_axes([.255, B0, .2, HH])]; axs += [fig.add_axes([.255 + k * .245, B0, .2, HH], sharey=axs[0]) for k in (1, 2)]
ax = axs[0]
dumb(ax, OWN.loc[ORD5, "atlas_J1"], OWN.loc[ORD5, "atlas_D1"], OWN.loc[ORD5, "atlas_J1_se"], OWN.loc[ORD5, "atlas_D1_se"], "+ measured annotations", "predicted effect alone", INK, "#b5b5b5", "effect τ* (± 95% CI)")
ax.set_yticks(yy); ax.set_yticklabels([f"{TN[t]} ({SHORT[OWN.loc[t, 'own_layer']]})" for t in ORD5]); ax.set_title("vs measured annotations", loc="left", fontsize=8)
ax.legend(fontsize=6.5, frameon=False, loc="upper left", bbox_to_anchor=(-.05, -.13)); panel(ax, "a", x=-.05, y=1.03)
# b: 활성도 조건화 세 구성(3p p04·p03): E(효과만), EA(이진 공변량), EAC(연속 공변량), S(활성도 층화 효과 주석). 채운 기호 = 11형질 BH q < 0.05
PV = pd.read_csv(R/"e15_activity_constructions/results/p04_values.csv").query("gwas == 'original'"); PS = pd.read_csv(R/"e15_activity_constructions/results/p03_main_summary.csv")
def pv(t, model, m): r = PV[(PV.trait == t) & (PV.model == model) & (PV.layer == m) & PV.kind.isin(["E", "S"])]; assert len(r) == 1, (t, model, m); return float(r.tau.iloc[0]), float(r.se.iloc[0])
def pq(t, model): r = PS[(PS.test == f"P-H1 {model}") & (PS.trait == t)]; assert len(r) == 1; return float(r.q.iloc[0])
ax = axs[1]; ret = []
MODELS = [("E", .3, "o", "#b5b5b5", "effect layers only"), ("EA", .1, "o", INK, "+ binary activity"), ("EAc", -.1, "s", "#2f6f9f", "+ continuous activity"), ("S", -.3, "D", "#8c4a2f", "stratified on activity")]
for t, y_ in zip(ORD5, yy):
    m = OWN.loc[t, "own_layer"]; e0 = pv(t, "E", m); assert abs(e0[0] - v5(t, "E", "E", m)[0]) < 1e-9
    ax.plot([min(pv(t, k, m)[0] for k, *_ in MODELS), max(pv(t, k, m)[0] for k, *_ in MODELS)], [y_, y_], color=GRID, lw=1.2, zorder=1)
    for k, dy, mk, c, lab in MODELS:
        v = pv(t, k, m); filled = (k == "E") or pq(t, k) < .05
        ax.errorbar(v[0], y_ + dy, xerr=1.96 * v[1], fmt=mk, ms=3.6 if k == "E" else 4, color=c, mfc=c if filled else "white", mec=c, mew=.9, lw=.8, capsize=0, zorder=3, label=lab if t == ORD5[0] else None)
    ret.append(pv(t, "EA", m)[0] / e0[0])
ax.axvline(0, color=MUTED, lw=.7); ax.tick_params(labelleft=False); ax.set_xlabel("effect τ* (± 95% CI)"); clean(ax); ax.set_ylim(-.7, len(ORD5) - .3)
ax.set_title("vs predicted activity", loc="left", fontsize=8); ax.legend(fontsize=6.5, frameon=False, loc="upper left", bbox_to_anchor=(-.05, -.13), ncol=2, columnspacing=.8, handletextpad=.3); panel(ax, "b", x=-.05, y=1.03)
for k in ("EA", "EAc", "S"): NUM[f"fig5b_n_fdr05_{k}"] = int(sum(pq(t, k) < .05 for t in ORD5))
NUM["fig5b_median_retention_EA"] = float(np.median(ret))
# c: 활성 요소 안(활성도 상위 10 %) 활성도 맞춤 효과 상·하 반분의 차이 Δ(3q q03, 사후). 채운 기호 = 11형질 BH q < 0.05
QS = pd.read_csv(R/"e16_within_active/results/q03_main_summary.csv"); QS = QS[QS.test.str.startswith("Q-H1") & QS.trait.isin(ORD5)].set_index("trait"); assert len(QS) == 11
ax = axs[2]
for t, y_ in zip(ORD5, yy):
    r = QS.loc[t]; assert r.layer == OWN.loc[t, "own_layer"]
    ax.errorbar(r.delta, y_, xerr=1.96 * r.se, fmt="o", ms=4.5, color=INK, mfc=INK if r.q < .05 else "white", mec=INK, mew=.9, lw=.8, capsize=0, label="larger − smaller effect,\nactivity-matched" if t == ORD5[0] else None)
ax.axvline(0, color=MUTED, lw=.7); ax.tick_params(labelleft=False); ax.set_xlabel("Δτ* (± 95% CI)"); clean(ax); ax.set_ylim(-.7, len(ORD5) - .3)
ax.set_title("within active elements (post hoc)", loc="left", fontsize=8); ax.legend(fontsize=6.5, frameon=False, loc="upper left", bbox_to_anchor=(-.05, -.13)); panel(ax, "c", x=-.05, y=1.03)
NUM["fig5c_n_fdr05"] = int((QS.q < .05).sum())
fig.savefig(OUT/"fig5.png", dpi=300); plt.close(fig)

# ------------------------------------------------------------------ Supp Fig S5 예측 대 실측·eQTL 크기 비교 (v1 그림 5b–d)
fig = plt.figure(figsize=(7.2, 4.4))
axs = [fig.add_axes([.255, B0, .2, HH])]; axs.append(fig.add_axes([.5, B0, .2, HH], sharey=axs[0])); axs.append(fig.add_axes([.8, B0, .15, HH]))
ax = axs[0]
dumb(ax, OWN.loc[ORD5, "atlas_S"], OWN.loc[ORD5, "meas_S"], OWN.loc[ORD5, "atlas_S_se"], OWN.loc[ORD5, "meas_S_se"], "Atlas", "measured, same layer", INK, "#c08a5a", "τ* (± 95% CI)")
ax.set_yticks(yy); ax.set_yticklabels([f"{TN[t]} ({SHORT[OWN.loc[t, 'own_layer']]})" for t in ORD5]); ax.set_title("vs measured", loc="left", fontsize=8)
for t, y_ in zip(ORD5, yy):
    if OWN.loc[t, "atlas_minus_meas_p2"] < .05: ax.text(1.0, y_, "*", transform=ax.get_yaxis_transform(), va="center", fontsize=9)
ax.legend(fontsize=6.5, frameon=False, loc="upper left", bbox_to_anchor=(-.05, -.13)); panel(ax, "a", x=-.05, y=1.03)
ax = axs[1]
dumb(ax, OWN.loc[ORD5, "atlas_D1"], OWN.loc[ORD5, "eqtl_D1"], OWN.loc[ORD5, "atlas_D1_se"], OWN.loc[ORD5, "eqtl_D1_se"], "Atlas (own layer)", "fine-mapped eQTL", INK, "#9b9b9b", "τ* (± 95% CI)")
ax.tick_params(labelleft=False); ax.set_title("vs eQTL", loc="left", fontsize=8)
for t, y_ in zip(ORD5, yy):
    if OWN.loc[t, "atlas_minus_eqtl_p2"] < .05: ax.text(1.0, y_, "*", transform=ax.get_yaxis_transform(), va="center", fontsize=9)
ax.legend(fontsize=6.5, frameon=False, loc="upper left", bbox_to_anchor=(-.05, -.13)); panel(ax, "b", x=-.05, y=1.03)
ax = axs[2]
for j, m in enumerate(LAY):
    d = SM[SM["mod"] == m]; xs = np.full(len(d), j) + np.linspace(-.2, .2, len(d))
    ax.scatter(xs, d["diff"], s=12, color=COL[m], alpha=.85, zorder=3, edgecolor="white", lw=.4)
    ax.plot([j - .3, j + .3], [d["diff"].median()] * 2, color=INK, lw=1.6, zorder=4)
    NUM[f"figS5c_{m}_median_diff"] = float(d["diff"].median()); NUM[f"figS5c_{m}_n_sig_pos"] = int(((d["diff"] > 0) & (d["diff_p2"] < .05)).sum())
ax.axhline(0, color=MUTED, lw=.7); ax.set_xticks(range(3)); ax.set_xticklabels([SHORT[m] for m in LAY], fontsize=7, rotation=35, ha="right"); ax.set_ylabel("Atlas − measured τ*")
ax.set_title("by layer", loc="left", fontsize=8); clean(ax); panel(ax, "c", x=-.42, y=1.03)
fig.savefig(OUT/"figS5.png", dpi=300); plt.close(fig)

# ------------------------------------------------------------------ Supp Fig S4 CAGE − DNase: a 포함 관계 역치(3m M2), b 연속·겹치지 않는 구간(3o)
M2 = pd.read_csv(R/"e12_diagnostics/results/m02_threshold_contrasts.csv").set_index(["trait", "threshold"])
O3 = pd.read_csv(R/"e14_continuous_bins/results/o03_contrasts_long.csv"); O3 = O3[O3.gwas == "original"].set_index(["trait", "criterion"])
TH4 = ["0.1%", "0.5%", "1%", "5%"]; TCOL = {"0.1%": "#c6dbef", "0.5%": "#6baed6", "1%": "#2171b5", "5%": "#08306b"}
CR5 = ["CONT", "B1", "B2", "B3", "B4"]; CRN = {"CONT": "continuous", "B1": "bin top 0–0.1%", "B2": "bin 0.1–0.5%", "B3": "bin 0.5–1%", "B4": "bin 1–5%"}
CCOL = {"CONT": "#333333", "B1": "#7f2704", "B2": "#d94801", "B3": "#fd8d3c", "B4": "#fdbe85"}
ROWS4 = ORD + ["BMD", "NEUT", "TC"]; ny4 = len(ROWS4) + 1; yp4 = {t: ny4 - 1 - i - (1 if t in ["BMD", "NEUT", "TC"] else 0) for i, t in enumerate(ROWS4)}
fig, axs = plt.subplots(1, 2, figsize=(7.2, 6.0), sharey=True); fig.subplots_adjust(left=.2, right=.98, top=.94, bottom=.17, wspace=.08)
ax = axs[0]
for j, q in enumerate(TH4):
    for t in ROWS4:
        r = M2.loc[(t, q)]; ax.errorbar(r.c1, yp4[t] + (1.5 - j) * .18, xerr=[[r.c1 - r.c1_lo], [r.c1_hi - r.c1]], fmt="o", ms=3.3, lw=.9, color=TCOL[q], capsize=0, label=f"top {q}" if t == ROWS4[0] else None)
ax.set_title("nested thresholds", loc="left", fontsize=8.5)
ax = axs[1]
for j, c in enumerate(CR5):
    for t in ROWS4:
        r = O3.loc[(t, c)]; ax.errorbar(r.c1, yp4[t] + (2 - j) * .15, xerr=1.96 * r.c1_se, fmt="D" if c == "CONT" else "o", ms=3.2, lw=.9, color=CCOL[c], capsize=0, label=CRN[c] if t == ROWS4[0] else None)
ax.set_title("continuous and non-overlapping bins", loc="left", fontsize=8.5)
for k, ax in enumerate(axs):
    ax.axvline(0, color=MUTED, lw=.7); ax.axhline(ny4 - 1 - len(ORD), color="#999", lw=.8, ls=":"); ax.set_xlabel("CAGE − DNase (τ* difference, ± 95% CI)"); clean(ax); panel(ax, "ab"[k], x=-.02, y=1.0)
    ax.legend(fontsize=6.8, frameon=False, loc="upper left", bbox_to_anchor=(-.02, -.085), ncol=2 if k == 0 else 3, columnspacing=1.0, handletextpad=.3)
axs[0].set_yticks([yp4[t] for t in ROWS4]); axs[0].set_yticklabels([TN[t] for t in ROWS4]); axs[0].set_ylim(-.7, ny4 - .4)
fig.savefig(OUT/"figS4.png", dpi=300); plt.close(fig)

# ------------------------------------------------------------------ Supp Fig S2 단일 주석 계통 × 층위 지도 (D1)
MODS6 = ["CAGE", "DNASE", "H3K27ac", "H3K4me1", "H3K4me3", "RNA"]
D1["col"] = D1.unit + "|" + D1["mod"]
cols = [f"{l}|{m}" for l in LIN12 for m in MODS6]
T = D1.pivot(index="trait", columns="col", values="tau_star").reindex(index=ORD, columns=cols); Q = D1.pivot(index="trait", columns="col", values="q_bh").reindex(index=ORD, columns=cols)
fig, ax = plt.subplots(figsize=(7.2, 3.6))
im = ax.imshow(T.values, cmap=DIV, norm=TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=4.5), aspect="auto")
for i in range(T.shape[0]):
    for j in range(T.shape[1]):
        if Q.iat[i, j] < .05: ax.text(j, i, "•", ha="center", va="center", fontsize=6, color=INK)
ax.set_xticks([k * 6 + 2.5 for k in range(12)]); ax.set_xticklabels([LN[l] for l in LIN12], rotation=40, ha="right", fontsize=7)
for k in range(1, 12): ax.axvline(k * 6 - .5, color="white", lw=1.6)
ax.set_yticks(range(len(ORD))); ax.set_yticklabels([TN[t] for t in ORD], fontsize=7.5); ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]
LAB6 = {"CAGE": "CAGE", "DNASE": "DNase", "H3K27ac": "K27ac", "H3K4me1": "K4me1", "H3K4me3": "K4me3", "RNA": "RNA"}
top = ax.secondary_xaxis("top"); top.set_xticks(range(len(cols))); top.set_xticklabels([LAB6[c.split("|")[1]] for c in cols], rotation=90, fontsize=4.3, color=INK)
top.tick_params(length=0, pad=1); [sp.set_visible(False) for sp in top.spines.values()]
cb = fig.colorbar(im, ax=ax, fraction=.025, pad=.01); cb.set_label("τ*"); cb.outline.set_visible(False)
fig.tight_layout(); fig.savefig(OUT/"figS2.png", dpi=300); plt.close(fig)
NUM["figS2_sig_cells"] = int((Q < .05).values.sum())
# ------------------------------------------------------------------ Supp Fig S3 역치별 1위 층위 (실험 3g)
from matplotlib.patches import Rectangle
G3 = pd.read_csv(R/"e06_thresholds/results/u03_main_layers.csv").set_index(["trait", "thr"])
THS = ["L0p1", "L0p5", "L1p0", "L5"]; THN = {"L0p1": "0.1%", "L0p5": "0.5%", "L1p0": "1%", "L5": "5%"}; QP = {"L0p1": "4 × 10$^{-4}$", "L0p5": "0.005", "L1p0": "0.055", "L5": "2 × 10$^{-7}$"}   # 3j 공동 jackknife (표 S7)
VALS = ["BMD", "NEUT", "TC"]; ROWS = ORD + VALS; ny = len(ROWS) + 1
yp = {t: ny - 1 - i - (1 if t in VALS else 0) for i, t in enumerate(ROWS)}
fig, ax = plt.subplots(figsize=(5.2, 5.4))
for j, q in enumerate(THS):
    for t in ROWS:
        r = G3.loc[(t, q)]; sig = r.top_vs_second_p1 < .05
        ax.add_patch(Rectangle((j + .04, yp[t] - .42), .92, .84, facecolor=COL[r.top], alpha=1 if sig else .4, edgecolor="white", lw=1.2))
        ax.text(j + .5, yp[t], SHORT[r.top] + (" *" if sig else ""), ha="center", va="center", fontsize=7, color=("white" if sig and r.top != "DNASE" else INK))
ax.set_xlim(0, 4); ax.set_ylim(-.6, ny - .4); ax.set_xticks(np.arange(4) + .5); ax.set_xticklabels([f"top {THN[q]}\nQ P = {QP[q]}" for q in THS], fontsize=7.5)
ax.set_yticks([yp[t] for t in ROWS]); ax.set_yticklabels([TN[t] for t in ROWS], fontsize=7.5); ax.tick_params(length=0)
ax.axhline(ny - 1 - len(ORD), color="#999", lw=.8, ls=":"); [s.set_visible(False) for s in ax.spines.values()]
fig.tight_layout(); fig.savefig(OUT/"figS3.png", dpi=300); plt.close(fig)
(P/"numbers.json").write_text(json.dumps(NUM, indent=1, default=float), encoding="utf-8"); print(json.dumps(NUM, indent=1, default=float))
