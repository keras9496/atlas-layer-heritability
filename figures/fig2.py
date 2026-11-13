# -*- coding: utf-8 -*-
"""Figure 2 — 조직 묶음 × 조절 층위 지도와 신호의 유형 (벤치마크 14형질 + 106개 GWAS).
 a 대표 형질 다섯의 띠(조현병, AD, 단핵구 수, LDL, 키). 띠 안에서 열 = 조직 묶음 12개(그림 1b와 같은 순서), 줄 = 조절 층위 6개.
   칸 색 = τ*(주석 하나를 baselineLD v2.2와 함께 넣은 모형), 점 = FDR < 0.05(14형질 × 72주석 = 1,008칸에 걸친 BH). 왼쪽 끝의 세로 글자와 선 = 유형.
   자료: e22_benchmark14/results/w04_map72.csv, 유형은 w04_types.csv와 대조한다.
 b 유형 지도: 원 = 106개 GWAS의 분포(e21_panel_single_layer_map/results/v02_gwas106.csv). 이름 = 패널 a의 다섯 형질의 자리와, 106개 GWAS 가운데 조직·층위 모두 넓은 유형의 보기 넷.
 유형의 색은 층위 색(파랑·주황·초록)과 열지도 색(올리브·남색)에 없는 계열이다: focal 회색, layer-spanning 황토, tissue-spanning 자홍, 둘 다 넓은 유형 보라.
 그림 안에는 결과와 이름표만 둔다(읽는 법·기준은 그림 설명). 기존 결과 파일만 읽는다. 새 통계 없음.
실행: python figures/fig2.py      (저장소 뿌리에서; 모두 한 번에는 python figures/make_figures.py)
출력: figures/out/fig2.png, fig2.pdf, numbers_fig2.json
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches
from matplotlib.colors import TwoSlopeNorm, to_hex
from matplotlib.lines import Line2D
import common as c
from common import INK, MUTED, LIN12, LN, COL
from taumap import TAUMAP, wcag

plt.rcParams.update({"font.family": "Arial", "font.size": 7, "text.color": INK, "axes.edgecolor": INK,
                     "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
                     "pdf.fonttype": 42, "axes.unicode_minus": True})
OUT = Path(__file__).resolve().parent / "out"; OUT.mkdir(parents=True, exist_ok=True)
FS, FSS, FS_LET = 7.0, 6.5, 9
ORD = ["EA", "N_SCZ", "N_BIP", "NEURO", "N_BMI", "N_AD", "MONO", "NEUT", "IBD", "RA", "N_LDL", "TC", "N_HEIGHT", "BMD"]      # 벤치마크 14형질
TN = {"EA": "Educational attainment", "N_SCZ": "Schizophrenia", "N_BIP": "Bipolar disorder", "NEURO": "Neuroticism", "N_BMI": "BMI", "N_AD": "Alzheimer's disease", "MONO": "Monocyte count",
      "NEUT": "Neutrophil count", "IBD": "IBD", "RA": "Rheumatoid arthritis", "N_LDL": "LDL cholesterol", "TC": "Total cholesterol", "N_HEIGHT": "Height", "BMD": "Bone mineral density"}
EX = ["N_SCZ", "N_AD", "MONO", "N_LDL", "N_HEIGHT"]              # 대표 형질 (벤치마크 GWAS)
MODS6 = ["CAGE", "DNASE", "H3K27ac", "H3K4me1", "H3K4me3", "RNA"]
LAB6 = {"CAGE": "CAGE", "DNASE": "DNase", "H3K27ac": "H3K27ac", "H3K4me1": "H3K4me1", "H3K4me3": "H3K4me3", "RNA": "RNA-seq"}
T_CUT, L_CUT = 3, 3.0                                             # e21_panel_single_layer_map v02_types.py와 같은 기준
# 유형의 색 (머리말 참고): 층위 색(CAGE 파랑, DNase 주황, H3K27ac 초록)과 열지도 색(올리브, 남색)에 없는 계열
TYPE_COL = {"focal": "#555555", "layer-spanning": "#C4A000", "tissue-spanning": "#D1478F", "tissue- and layer-spanning": "#6A2C91"}
RESERVED = {"CAGE": COL["CAGE"], "DNase": COL["DNASE"], "H3K27ac": COL["H3K27ac"], "tau* > 0 pole": "#757b10", "tau* < 0 pole": "#797ec7"}
TYPE_EN = {"focal": "Focal", "layer-spanning": "Layer-spanning", "tissue-spanning": "Tissue-spanning", "tissue- and layer-spanning": "Tissue- and layer-spanning"}

# ------------------------------------------------------------------ 자료
W4 = c.R / "e22_benchmark14" / "results"
D1 = pd.read_csv(W4 / "w04_map72.csv")
assert len(D1) == 1008 and D1.groupby(["trait", "group", "layer"]).size().max() == 1
assert set(D1.group) == set(LIN12) and set(D1.layer) == set(MODS6) and set(D1.trait) == set(ORD)
N_SIG = int((D1.q < .05).sum()); assert N_SIG == 107, N_SIG                                                # w04_summary.txt: 유의 107칸
N_RNA = int(((D1.q < .05) & (D1.layer == "RNA")).sum()); assert N_RNA == 2
def grid(t, col):
    return D1[D1.trait == t].pivot(index="layer", columns="group", values=col).reindex(index=MODS6, columns=LIN12)
TAUS = {t: grid(t, "tau") for t in EX}; QS = {t: grid(t, "q") for t in EX}
assert all(g.shape == (6, 12) and not g.isna().any().any() for g in TAUS.values())
VMIN, VMAX = -2.5, 4.5
assert all(VMIN <= g.values.min() and g.values.max() <= VMAX for g in TAUS.values()), {t: (g.values.min(), g.values.max()) for t, g in TAUS.items()}
NORM = TwoSlopeNorm(vmin=VMIN, vcenter=0, vmax=VMAX)
def typ(K, T, L):
    if K == 0: return "none"
    return ("tissue- and layer-spanning" if L >= L_CUT else "tissue-spanning") if T >= T_CUT else ("layer-spanning" if L >= L_CUT else "focal")
IDX = {}
for t in ORD:
    d = D1[(D1.trait == t) & (D1.q < .05)]; K = len(d); T = d.group.nunique(); L = K / T if T else float("nan"); IDX[t] = dict(K=K, T=T, L=L, type=typ(K, T, L))
V3 = pd.read_csv(W4 / "w04_types.csv", index_col=0)                                                        # w04_bench.py의 표와 대조
assert all(IDX[t]["K"] == V3.K[t] and IDX[t]["T"] == V3["T"][t] and IDX[t]["type"] == V3.type[t] for t in ORD)
TYPE_OF = {t: IDX[t]["type"] for t in EX}
assert [TYPE_OF[t] for t in EX] == ["focal", "focal", "layer-spanning", "tissue-spanning", "tissue- and layer-spanning"]
assert (IDX["N_SCZ"]["T"], IDX["N_SCZ"]["L"]) == (IDX["N_AD"]["T"], IDX["N_AD"]["L"]) == (2, 2.5) and (IDX["N_HEIGHT"]["T"], IDX["N_HEIGHT"]["L"]) == (6, 3.0)
NTYPE14 = {ty: sum(IDX[t]["type"] == ty for t in ORD) for ty in list(TYPE_COL) + ["none"]}
assert NTYPE14 == {"focal": 7, "layer-spanning": 2, "tissue-spanning": 3, "tissue- and layer-spanning": 2, "none": 0}

# ------------------------------------------------------------------ 배치 (인치)
W = 7.2
HX0, HW = 2.08, 4.40                                              # 그림 1b의 히트맵 왼쪽 끝과 폭
RH, GAPB = 0.105, 0.11
BAND = 6 * RH
BH_, BGAP = 2.02, 0.26                                            # 패널 b의 축 높이, a의 열 이름과 b 사이
MB = 0.98                                                         # b 아래(축 이름)
YB = MB + BH_ + BGAP + 0.98                                       # a의 마지막 띠 아래 끝(돌린 열 이름 0.98)
MT = 0.22
H = YB + len(EX) * BAND + (len(EX) - 1) * GAPB + MT
fig = plt.figure(figsize=(W, H))
def axin(l, b, w, h): return fig.add_axes([l / W, b / H, w / W, h / H])
ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(0, W); ov.set_ylim(0, H); ov.axis("off"); ov.set_zorder(5)
n_dots = 0; DOTS = {}; Y0 = {}
for k, t in enumerate(EX):
    y0 = YB + (len(EX) - 1 - k) * (BAND + GAPB); Y0[t] = y0
    ax = axin(HX0, y0, HW, BAND)
    ax.imshow(TAUS[t].values, cmap=TAUMAP, norm=NORM, aspect="auto", extent=(-.5, 11.5, 5.5, -.5), interpolation="nearest")
    for x in np.arange(.5, 11.5): ax.axvline(x, color="white", lw=.7)
    for y in np.arange(.5, 5.5): ax.axhline(y, color="white", lw=.45)
    ii, jj = np.where(QS[t].values < .05); n_dots += len(ii)
    DOTS[t] = [f"{LIN12[j]}|{MODS6[i]}" for i, j in zip(ii, jj)]; assert len(ii) == IDX[t]["K"]
    assert all(wcag(INK, TAUMAP(NORM(TAUS[t].values[i, j]))) >= 3.0 for i, j in zip(ii, jj))
    ax.plot(jj, ii, ls="none", marker="o", ms=2.8, mfc=INK, mec="white", mew=.5, zorder=4)
    ax.set_xlim(-.5, 11.5); ax.set_ylim(5.5, -.5)
    ax.set_yticks(range(6)); ax.set_yticklabels([LAB6[m] for m in MODS6], fontsize=FSS, color=MUTED)
    ax.tick_params(length=0, pad=3); [sp.set_visible(False) for sp in ax.spines.values()]
    if k == len(EX) - 1:
        ax.set_xticks(range(12)); ax.set_xticklabels([LN[l] for l in LIN12], rotation=40, ha="right", rotation_mode="anchor", fontsize=FS, color=INK)
    else: ax.set_xticks([])
    ov.text(HX0 - .62, y0 + BAND / 2, TN[t].replace("Alzheimer's disease", "Alzheimer's\ndisease").replace("LDL cholesterol", "LDL\ncholesterol").replace("Monocyte count", "Monocyte\ncount"),
            fontsize=FS, color=INK, ha="right", va="center", linespacing=1.1)
# 유형 이름(세로)과 묶는 선
TYPE_LAB = []
for ty in ["focal", "layer-spanning", "tissue-spanning", "tissue- and layer-spanning"]:
    ts_ = [t for t in EX if TYPE_OF[t] == ty]; lo = Y0[ts_[-1]]; hi = Y0[ts_[0]] + BAND
    ov.plot([.40, .40], [lo, hi], color=TYPE_COL[ty], lw=2.2, solid_capstyle="butt")
    two = ty == "tissue- and layer-spanning"                         # 띠 하나보다 긴 이름은 두 줄로
    TYPE_LAB.append(ov.text(.23 if two else .29, (lo + hi) / 2, "Tissue- and\nlayer-spanning" if two else TYPE_EN[ty], rotation=90, ha="center", va="center", fontsize=FSS, color=INK, linespacing=1.1))
# 색 막대
CBX, CBW, CBH = 6.70, 0.085, 1.45; ytop = H - MT
cb = axin(CBX, ytop - CBH, CBW, CBH); v = np.linspace(VMIN, VMAX, 512)
cb.imshow(TAUMAP(NORM(v))[::-1, None, :3], aspect="auto", extent=(0, 1, VMIN, VMAX))
cb.set_xticks([]); cb.yaxis.tick_right(); cb.set_yticks(range(-2, 5)); cb.set_yticklabels([("0" if q == 0 else f"{q:+d}".replace("-", "−")) for q in range(-2, 5)], fontsize=FS, color=INK)
cb.tick_params(axis="y", length=2, width=.5, color=MUTED, pad=2); [sp.set_visible(False) for sp in cb.spines.values()]; cb.set_ylim(VMIN, VMAX)
ov.text(CBX + CBW / 2, ytop + .05, r"$\tau$*", fontsize=FS, color=INK, ha="center", va="bottom", math_fontfamily="stixsans")
ky = ytop - CBH - .22
ov.plot([CBX - .03], [ky], ls="none", marker="o", ms=2.8, mfc=INK, mec="white", mew=.5)
ov.text(CBX + .05, ky, "FDR\n< 0.05", fontsize=FSS, color=MUTED, ha="left", va="center", linespacing=1.15)
ov.text(.04, H - .04, "a", fontsize=FS_LET, fontweight="bold", ha="left", va="top", color=INK)

# ------------------------------------------------------------------ b 유형 지도 (전체 분포 위에 대표 형질을 표시한다)
#   가로 = 유의한 주석이 있는 조직 묶음의 수 T, 세로 = 조직 묶음 하나에 유의한 층위의 평균 수 L = K / T. 축은 0부터. 옅은 바탕색 = 유형의 영역(경계 T = 2.5, L = 3), 원의 색 = 유형.
#   원 = 106개 GWAS의 분포(e21_panel_single_layer_map v02_gwas106.csv; 7,632칸에 걸친 FDR): 원 하나 = 같은 자리의 GWAS들, 넓이 ∝ GWAS 수. 유의한 주석이 없는 GWAS는 원점의 빈 원.
#   이름 = 패널 a의 대표 형질 다섯(자리는 벤치마크 GWAS, 14형질 1,008칸에 걸친 FDR — 패널 a의 점과 같은 기준; 표지는 그리지 않는다)과 조직·층위 모두 넓은 유형의 보기 넷(106개 GWAS).
XC, YC = T_CUT - .5, L_CUT
BX0, BW = 0.62, 6.40
G6 = pd.read_csv(c.R / "e21_panel_single_layer_map/results/v02_gwas106.csv", index_col=0); assert len(G6) == 106
assert all(typ(int(r.K), int(r.T), r.L) == r.type for r in G6.itertuples())                              # v02_types.py의 유형 = 이 스크립트의 규칙
NTYPE = {ty: int((G6.type == ty).sum()) for ty in list(TYPE_COL) + ["none"]}; assert sum(NTYPE.values()) == 106
XMAX = int(max(G6["T"].max(), max(IDX[t]["T"] for t in EX))); XLIM, YLIM = (-0.42, XMAX + .7), (-0.30, 5.30)
bx = axin(BX0, MB, BW, BH_)
XL0 = .15                                                         # 왼쪽 두 영역의 왼쪽 끝(두 줄 이름 "Schizophrenia, Alzheimer's disease"가 영역 안에 들도록 0.15)
REG = {"focal": (XL0, XC, .5, YC), "layer-spanning": (XL0, XC, YC, YLIM[1]), "tissue-spanning": (XC, XLIM[1], .5, YC), "tissue- and layer-spanning": (XC, XLIM[1], YC, YLIM[1])}
for ty, (a0, a1, b0, b1) in REG.items():
    bx.add_patch(matplotlib.patches.Rectangle((a0, b0), a1 - a0, b1 - b0, facecolor=TYPE_COL[ty], alpha=.11, edgecolor="none", zorder=0))
bx.plot([XC, XC], [.5, YLIM[1]], color="white", lw=1.6, zorder=.5, solid_capstyle="butt"); bx.plot([XL0, XLIM[1]], [YC, YC], color="white", lw=1.6, zorder=.5, solid_capstyle="butt")
def area_c(n): return 11.0 * n                                    # 원의 넓이(scatter의 s, pt²) ∝ 그 자리의 GWAS 수
PC = G6[G6.K > 0].groupby(["T", "L"]).agg(n=("K", "size"), type=("type", "first")).reset_index(); assert PC.n.sum() == 106 - NTYPE["none"]
for r in PC.sort_values("n", ascending=False).itertuples():
    bx.scatter([r.T], [r.L], s=area_c(r.n), facecolor=TYPE_COL[r.type], edgecolor="white", linewidth=.5, zorder=3)
bx.plot([0], [0], ls="none", marker="o", ms=np.sqrt(area_c(NTYPE["none"])), mfc="white", mec=MUTED, mew=.9, zorder=3, clip_on=False)
# 이름(표지는 그리지 않고, 대략의 자리에 이름만 둔다)
#  (1) 패널 a의 대표 형질 다섯: 자리 = 패널 a에 보이는 유의한 주석에서 센 값(주 GWAS). 그 자리에 원이 있으면 원 옆에, 없으면(키) 그 자리에 이름을 둔다.
#  (2) 조직·층위 모두 넓은 유형(패널 a에 없는 유형)의 보기: 106개 GWAS 가운데 넷. 원 옆에 이름을 둔다. x = 8의 두 원(무릎 골관절염이 위, 발꿈치 골밀도가 아래)은 두 줄로.
def rpt(T_, L_):                                                  # 그 자리에 있는 원의 반지름(pt); 원이 없으면 0
    m = PC[(PC["T"] == T_) & np.isclose(PC.L, L_)]; return float(np.sqrt(area_c(int(m.n.iloc[0]))) / 2) if len(m) else 0.0
PURPLE_EX = {"T063": "Total protein", "T005": "COPD", "T004": "Knee osteoarthritis", "T070": "Heel BMD"}   # id → 그림의 이름 (보충 표의 이름: Total protein, Chronic obstructive pulmonary disease, Knee osteoarthritis, Heel bone mineral density)
assert all(G6.type[i] == "tissue- and layer-spanning" for i in PURPLE_EX)
assert G6["T"]["T004"] == G6["T"]["T070"] and G6.L["T004"] > G6.L["T070"]                                  # 두 줄 이름의 위·아래 = 원의 위·아래
GAPN = 2.6                                                        # 원의 가장자리와 이름 사이(pt)
def side(T_, L_, where, dx=0.0, dy=0.0):
    """이름의 자리: where = left | right | above | below | at. 반환 = annotate에 줄 값들."""
    r = rpt(T_, L_) + GAPN
    return {"left": dict(xytext=(-r + dx, dy), ha="right", va="center"), "right": dict(xytext=(r + dx, dy), ha="left", va="center"),
            "above": dict(xytext=(dx, r - 1 + dy), ha="center", va="bottom"), "below": dict(xytext=(dx, -r + 1 + dy), ha="center", va="top"), "at": dict(xytext=(dx, dy), ha="center", va="center")}[where]
NAMES = [("N_SCZ+N_AD", "Schizophrenia,\nAlzheimer's disease", IDX["N_SCZ"]["T"], IDX["N_SCZ"]["L"], "left", 0, 1.5),
         ("MONO", TN["MONO"], IDX["MONO"]["T"], IDX["MONO"]["L"], "left", 0, 0), ("N_LDL", TN["N_LDL"], IDX["N_LDL"]["T"], IDX["N_LDL"]["L"], "below", 0, 0),
         ("N_HEIGHT", TN["N_HEIGHT"], IDX["N_HEIGHT"]["T"], IDX["N_HEIGHT"]["L"], "above", 0, 0),
         ("T063", PURPLE_EX["T063"], G6["T"]["T063"], G6.L["T063"], "right", 0, 0), ("T005", PURPLE_EX["T005"], G6["T"]["T005"], G6.L["T005"], "right", 0, 0),
         ("T004+T070", PURPLE_EX["T004"] + "\n" + PURPLE_EX["T070"], G6["T"]["T004"], (G6.L["T004"] + G6.L["T070"]) / 2, "left", -1.0, 1.5)]
assert rpt(IDX["N_HEIGHT"]["T"], IDX["N_HEIGHT"]["L"]) == 0 and all(rpt(IDX[t]["T"], IDX[t]["L"]) > 0 for t in ("N_SCZ", "MONO", "N_LDL"))   # 키의 자리에는 106개 GWAS의 원이 없다
LAB = {}; OWN = {}
for key, s, x, y, where, dx, dy in NAMES:
    kw = side(x, y, where, dx, dy)
    if key == "T004+T070": kw["xytext"] = (-(rpt(G6["T"]["T004"], G6.L["T004"]) + GAPN) + dx, dy)
    LAB[key] = bx.annotate(s, (x, y), textcoords="offset points", fontsize=FS, color=INK, linespacing=1.12, zorder=6, **kw)
    OWN[key] = [(G6["T"][i], G6.L[i]) for i in key.split("+")] if key.startswith("T0") else [(x, y)]
bx.set_xlim(*XLIM); bx.set_ylim(*YLIM); bx.set_xticks(range(0, XMAX + 1)); bx.set_yticks(range(0, 6))
for s_ in ("top", "right"): bx.spines[s_].set_visible(False)
for s_ in ("left", "bottom"): bx.spines[s_].set_color(MUTED); bx.spines[s_].set_linewidth(.6)
bx.spines["bottom"].set_bounds(0, XMAX); bx.spines["left"].set_bounds(0, 5)
bx.tick_params(labelsize=FS, length=2.5, width=.6, color=MUTED, labelcolor=INK, pad=2)
bx.set_xlabel("Tissue groups with a significant annotation", fontsize=FS, color=INK, labelpad=3)
bx.set_ylabel("Significant layers\nper tissue group", fontsize=FS, color=INK, labelpad=4, linespacing=1.15)
# ---- 풀이 (그림 아래쪽 두 줄): 위 = 기호와 크기, 아래 = 색
KY1, KY2 = .13, .40
fig.canvas.draw(); rd0 = fig.canvas.get_renderer()
def tw(s): t_ = ov.text(0, 0, s, fontsize=FSS); w_ = t_.get_window_extent(rd0).width / fig.dpi; t_.remove(); return w_
sz = [(v, np.sqrt(area_c(v))) for v in (1, 5, 20)]; head = "Number of GWAS"
w_size = tw(head) + .14 + sum(d / 72 + .04 + tw(str(v)) for v, d in sz) + .13 * (len(sz) - 1)
kx = W / 2 - w_size / 2
ov.text(kx, KY2, head, fontsize=FSS, color=MUTED, ha="left", va="center"); kx += tw(head) + .14
for v, d in sz:                                                   # 빈 원(테두리만): 유형의 색과 헷갈리지 않게
    ov.plot([kx + d / 144], [KY2], ls="none", marker="o", ms=d, mfc="white", mec=INK, mew=.7, clip_on=False)
    ov.text(kx + d / 72 + .04, KY2, str(v), fontsize=FSS, color=MUTED, ha="left", va="center"); kx += d / 72 + .04 + tw(str(v)) + .13
# 색 풀이: 유형 이름 뒤에 그 유형의 GWAS 수(106개가 모두 이 지도에 있다는 것이 보이게)
items = [(f"{TYPE_EN[ty]} ({NTYPE[ty]})", dict(marker="o", ms=5.2, mfc=TYPE_COL[ty], mec="white", mew=.5)) for ty in ("focal", "layer-spanning", "tissue-spanning", "tissue- and layer-spanning")]
items.append((f"No significant annotation ({NTYPE['none']})", dict(marker="o", ms=4.2, mfc="white", mec=MUTED, mew=.9)))
GAPK = .22; wtot = sum(.11 + tw(s) for s, _ in items) + GAPK * (len(items) - 1); kx = W / 2 - wtot / 2; assert wtot < W - .30, wtot
for s, kw in items:
    ov.plot([kx], [KY1], ls="none", clip_on=False, **kw); ov.text(kx + .11, KY1, s, fontsize=FSS, color=INK, ha="left", va="center"); kx += .11 + tw(s) + GAPK
ov.text(.04, MB + BH_ + .17, "b", fontsize=FS_LET, fontweight="bold", ha="left", va="top", color=INK)


# ------------------------------------------------------------------ 점검 후 저장
fig.canvas.draw(); rd = fig.canvas.get_renderer(); fb = fig.bbox
ts = [t_ for t_ in fig.findobj(matplotlib.text.Text) if t_.get_visible() and t_.get_text().strip()]
ext = [(t_, t_.get_window_extent(rd)) for t_ in ts]
clipped = [t_.get_text() for t_, e in ext if e.x0 < fb.x0 or e.x1 > fb.x1 or e.y0 < fb.y0 or e.y1 > fb.y1]
assert not clipped, clipped
assert all(6.5 <= t_.get_fontsize() <= 7 or (t_.get_fontsize() == FS_LET and t_.get_fontweight() == "bold") for t_ in ts)
assert all(to_hex(t_.get_color()) in (INK, MUTED) for t_ in ts)
flat = [(t_, e) for t_, e in ext if t_.get_rotation() % 90 == 0]            # 세로 글자(유형 이름)도 상자가 축에 나란하다
ov_ = [(a.get_text(), b.get_text()) for i, (a, x) in enumerate(flat) for b, y in flat[i + 1:] if x.overlaps(y)]
assert not ov_, ov_
# b: 이름이 어느 원과도 겹치지 않고(자기 원 포함), 자기 원이 아닌 원보다 자기 원에 더 가깝다. 이름은 축 안쪽에 있다.
axb = bx.get_window_extent(rd); GAPS = {}
for key, an in LAB.items():
    e = an.get_window_extent(rd); assert axb.x0 <= e.x0 and e.x1 <= axb.x1 and axb.y0 <= e.y0 and e.y1 <= axb.y1, key
    def gap(T_, L_, n_):                                           # 이름 상자와 원 사이의 거리(pt; 음수 = 겹침)
        px, py = bx.transData.transform((T_, L_)); dx_ = max(e.x0 - px, 0, px - e.x1); dy_ = max(e.y0 - py, 0, py - e.y1)
        return (np.hypot(dx_, dy_) * 72 / fig.dpi) - np.sqrt(area_c(n_)) / 2
    g = {(r_.T, round(r_.L, 3)): gap(r_.T, r_.L, r_.n) for r_ in PC.itertuples()}
    own = [k_ for k_ in g if any(k_[0] == a_ and np.isclose(k_[1], b_, atol=1e-3) for a_, b_ in OWN[key])]
    oth = min(v_ for k_, v_ in g.items() if k_ not in own)
    assert min(g.values()) > .8, (key, min(g.items(), key=lambda kv: kv[1]))
    if own: assert max(g[k_] for k_ in own) < oth, (key, [g[k_] for k_ in own], oth)
    if key == "N_SCZ+N_AD": assert bx.transData.transform((XL0, 0))[0] <= e.x0, "이름이 영역의 왼쪽 끝 밖으로 나간다"
    GAPS[key] = dict(own=[round(float(g[k_]), 1) for k_ in own], nearest_other=round(float(oth), 1))
print("b: name-to-circle gaps (pt):", GAPS)
fig.savefig(OUT / "fig2.png", dpi=300); fig.savefig(OUT / "fig2.pdf"); plt.close(fig)

NUM = {"source": "e22_benchmark14/results/w04_map72.csv: tau, q (BH over the 1,008 cells of the 14 benchmark traits x 72 annotations; one-sided tau* > 0); types by the rule of e21_panel_single_layer_map v02_types.py (checked against w04_types.csv)",
       "a_traits": EX, "a_trait_names": {t: TN[t] for t in EX}, "a_types": TYPE_OF, "rows": MODS6, "columns": LIN12, "colour_scale": {"vmin": VMIN, "vcenter": 0, "vmax": VMAX},
       "n_cells": 1008, "n_significant_all_14_traits": N_SIG, "n_significant_RNA_seq": N_RNA, "a_n_dots": n_dots, "a_dots": DOTS,
       "a_tau_star": {t: {l: {m: round(float(TAUS[t].loc[m, l]), 3) for m in MODS6} for l in LIN12} for t in EX},
       "b_rule": {"T": "tissue groups with >= 1 significant annotation", "L": "significant annotations / T", "cuts": {"T": ">= 3", "L": ">= 3"}, "shaded_regions_left_edge": XL0},
       "b_distribution": {"source": "e21_panel_single_layer_map/results/v02_gwas106.csv (106 GWAS; BH across the 7,632 models)", "type_counts": NTYPE, "positions": int(len(PC)), "largest_circle": int(PC.n.max()), "circle_area": "proportional to the number of GWAS at the position (11 pt^2 each)", "open_circle_at_origin": NTYPE["none"]},
       "b_names_traits_of_a": {t: {k: (None if isinstance(v, float) and v != v else v) for k, v in IDX[t].items()} for t in EX},
       "b_names_traits_of_a_source": "benchmark GWAS (w04_map72.csv; BH across the 1,008 cells): positions counted from the significant annotations shown in panel a; names only, no marker; schizophrenia and Alzheimer's disease share one position; no circle of the 106 GWAS lies at the position of height",
       "b_names_tissue_and_layer_spanning_examples": {i: dict(name=s_, file=str(G6.file[i]), K=int(G6.K[i]), T=int(G6["T"][i]), L=float(G6.L[i])) for i, s_ in PURPLE_EX.items()},
       "benchmark14": {t: {k: (None if isinstance(v, float) and v != v else v) for k, v in IDX[t].items()} | {"name": TN[t]} for t in ORD}, "benchmark14_type_counts": NTYPE14,
       "type_colours": TYPE_COL, "name_to_circle_gaps_pt": GAPS}
(OUT / "numbers_fig2.json").write_text(json.dumps(NUM, indent=1, ensure_ascii=False), encoding="utf-8")
print(f"saved {OUT / 'fig2.png'} ({W} x {H:.2f} in); a dots {n_dots}; b: 106 GWAS {NTYPE}, names of the traits of a at " + ", ".join(f"{t} ({IDX[t]['T']}, {IDX[t]['L']:.2f})" for t in EX))
