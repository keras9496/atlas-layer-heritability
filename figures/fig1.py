# -*- coding: utf-8 -*-
"""Figure 1 — 주석 만들기와 조직 묶음 지도 (벤치마크 14형질).
 a 모식도: 대립유전자마다의 예측 신호 → 조직 묶음 × 층위 주석(예측 효과 상위 5%) → S-LDSC. 트랙 수와 biosample 수(1,312, 433)는 e02의 z01_lineages.txt와 대조한다.
 b 조직 묶음 지도: 벤치마크 14형질 × 조직 묶음 12개의 τ*, S-LDSC 모형 M1(baselineLD v2.2 + 조직 묶음 주석 하나).
   자료: e22_benchmark14/results/w04_tissue.csv (tau, q = 168칸에 걸친 BH, 단측 τ* > 0). 색 범위 −2.5 ~ +4.5. q < 0.05인 칸에만 숫자.
   줄 순서와 왼쪽 띠: 형질 범주의 조직 묶음 — brain tissue(교육 연수, 조현병, 양극성장애, 신경증, BMI), myeloid(AD, 단핵구 수, 호중구 수),
   lymphoid(IBD, 류마티스관절염), liver(LDL, 총콜레스테롤), stromal(키, BMD). 띠는 색 없이 묶음 이름만 단다.
 기존 결과 파일만 읽는다. 새 통계 없음.
실행: python figures/fig1.py      (저장소 뿌리에서; 모두 한 번에는 python figures/make_figures.py)
출력: figures/out/fig1.png, fig1.pdf, numbers_fig1.json
"""
import sys, re, json, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))      # -I 로 실행하면 스크립트 폴더가 sys.path에 없다
import numpy as np, pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, TwoSlopeNorm, to_rgb
from matplotlib.patches import Rectangle, FancyBboxPatch
import common as c
from common import COL, COL_OTHER, INK, MUTED, GRID, LIN12, LN

plt.rcParams.update({"font.family": "Arial", "font.size": 7, "text.color": INK, "axes.unicode_minus": True, "mathtext.fontset": "custom",
                     "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "pdf.fonttype": 42})
from matplotlib import _mathtext as _mt
def _xheight_from_glyph(self, fontname, fontsize, dpi):
    """Arial의 PCLT x-height를 약 1.9배 크게 잡아 위·아래첨자가 어긋나므로 'x' 글자 높이를 잰다."""
    return self.get_metrics(fontname, matplotlib.rcParams["mathtext.default"], "x", fontsize, dpi).iceberg
_mt.TruetypeFonts.get_xheight = _xheight_from_glyph
OUT = Path(__file__).resolve().parent / "out"; OUT.mkdir(parents=True, exist_ok=True)

# ================================================================== 자료 (b): 벤치마크 14형질
W4 = c.R / "e22_benchmark14" / "results"
ORD = ["EA", "N_SCZ", "N_BIP", "NEURO", "N_BMI", "N_AD", "MONO", "NEUT", "IBD", "RA", "N_LDL", "TC", "N_HEIGHT", "BMD"]
TN = {"EA": "Educational attainment", "N_SCZ": "Schizophrenia", "N_BIP": "Bipolar disorder", "NEURO": "Neuroticism", "N_BMI": "BMI", "N_AD": "Alzheimer's disease", "MONO": "Monocyte count",
      "NEUT": "Neutrophil count", "IBD": "IBD", "RA": "Rheumatoid arthritis", "N_LDL": "LDL cholesterol", "TC": "Total cholesterol", "N_HEIGHT": "Height", "BMD": "Bone mineral density"}
OWN = {"EA": "CNS_tissue", "N_SCZ": "CNS_tissue", "N_BIP": "CNS_tissue", "NEURO": "CNS_tissue", "N_BMI": "CNS_tissue", "N_AD": "myeloid", "MONO": "myeloid", "NEUT": "myeloid",
       "IBD": "lymphoid", "RA": "lymphoid", "N_LDL": "liver", "TC": "liver", "N_HEIGHT": "stromal", "BMD": "stromal"}
M1 = pd.read_csv(W4 / "w04_tissue.csv")
assert len(M1) == 168 and M1.groupby(["trait", "group"]).size().max() == 1 and set(M1.trait) == set(ORD) and set(M1.group) == set(LIN12)
T = M1.pivot(index="trait", columns="group", values="tau").reindex(index=ORD, columns=LIN12)
Q = M1.pivot(index="trait", columns="group", values="q").reindex(index=ORD, columns=LIN12)
assert T.shape == (14, 12) and not T.isna().any().any() and not Q.isna().any().any()
assert all(bool(M1[(M1.trait == t) & (M1.group == OWN[t])].own.iloc[0]) for t in ORD) and int(M1.own.sum()) == 14      # 범주의 조직 묶음 = 표의 own 열
SIG = Q.values < 0.05
assert SIG.sum() == 41, SIG.sum()                                   # w04_summary.txt: FDR < 0.05인 칸 41개 / 168
assert (T.values[SIG] > 0).all()                                    # 단측 검정(τ* > 0): 적힌 값은 모두 양수
VMIN, VMAX = -2.5, 4.5                                              # 색 범위(이 저장소의 열지도는 모두 같은 범위를 쓴다)
assert VMIN <= T.values.min() and T.values.max() <= VMAX, (T.values.min(), T.values.max())   # 색 눈금에서 잘리는 값 없음
def tau(t, l): return float(T.loc[t, l])
SUMMARY_VALUES = {("EA", "CNS_tissue"): 0.87, ("N_SCZ", "CNS_tissue"): 1.42, ("N_BIP", "CNS_tissue"): 1.76, ("NEURO", "CNS_tissue"): 1.42, ("N_BMI", "CNS_tissue"): 0.83,
                  ("N_AD", "myeloid"): 2.33, ("MONO", "myeloid"): 2.82, ("NEUT", "myeloid"): 2.04, ("IBD", "lymphoid"): 3.43, ("RA", "lymphoid"): 2.48,
                  ("N_LDL", "liver"): 2.65, ("TC", "liver"): 1.97, ("N_HEIGHT", "stromal"): 1.35, ("N_HEIGHT", "adipose"): 1.79, ("BMD", "stromal"): 1.85, ("BMD", "adipose"): 2.59}
for (t, l), v in SUMMARY_VALUES.items():                            # w04_summary.txt 1절에 적힌 값(소수 둘째 자리)
    assert abs(tau(t, l) - v) <= 0.005 + 1e-9, (t, l, tau(t, l), v)
RANK1 = {t: T.loc[t].idxmax() for t in ORD}; N_OWN_FIRST = sum(RANK1[t] == OWN[t] for t in ORD)
assert N_OWN_FIRST == 12 and RANK1["N_HEIGHT"] == RANK1["BMD"] == "adipose"                      # 범주 묶음이 1위인 형질 12개; 키와 BMD는 지방이 1위
N_ENR = {t: int(M1[(M1.trait == t)].enriched.sum()) for t in ORD}                                # 농축된 묶음(P < 0.05/12)의 수
# 그림 설명에 쓰는 수(모식도 안에는 적지 않는다)
N_SNP, N_TRACK, N_BIOS, N_BASE = 5_828_323, 1312, 433, 97
assert f"{N_SNP / 1e6:.2f}" == "5.83"
tab = (c.R / "e02_lineage_map/results/z01_lineages.txt").read_text(encoding="utf-8").splitlines()
rows = [l.split() for l in tab if l.split() and l.split()[0] in LIN12 and len(l.split()) == 8]
cnt = np.array([[int(x) for x in r[1:]] for r in rows])
assert len(rows) == 12 and (cnt[:, :6] > 0).all() and cnt[:, :6].sum() == N_TRACK and cnt[:, 6].sum() == N_BIOS   # 1,312 트랙, 433 biosample, 빈 칸 없음

# ================================================================== 열지도 색 (발산 쌍: 올리브 = τ* > 0, 남색 = τ* < 0; taumap.py와 같다)
POS_HEX, NEG_HEX, MID_HEX = "#757b10", "#797ec7", "#efefef"
def _s2l(x): return x / 12.92 if x <= 0.04045 else ((x + 0.055) / 1.055) ** 2.4
def _l2s(x): x = min(max(x, 0.0), 1.0); return 12.92 * x if x <= 0.0031308 else 1.055 * x ** (1 / 2.4) - 0.055
def oklab(hexc):
    r, g, b = (_s2l(v) for v in to_rgb(hexc))
    l, m, s = (np.cbrt(0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b), np.cbrt(0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b),
               np.cbrt(0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b))
    return np.array([0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s, 1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
                     0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s])
def oklab_to_srgb(L, a, b):
    l_, m_, s_ = L + 0.3963377774 * a + 0.2158037573 * b, L - 0.1055613458 * a - 0.0638541728 * b, L - 0.0894841775 * a - 1.2914855480 * b
    l, m, s = l_ ** 3, m_ ** 3, s_ ** 3
    lin = [4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s, -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
           -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s]
    return [_l2s(x) for x in lin], max(max(0 - x, x - 1, 0) for x in lin)
def arm(pole, n=128, gamma=1.0):
    """중립 회색 → 극 색: 명도는 직선, 채도는 t**gamma."""
    L0 = oklab(MID_HEX)[0]; Lp, ap, bp = oklab(pole); Cp, hh = math.hypot(ap, bp), math.atan2(bp, ap); out = []
    for t in np.linspace(0, 1, n):
        rgb, oog = oklab_to_srgb(L0 + (Lp - L0) * t, Cp * t ** gamma * math.cos(hh), Cp * t ** gamma * math.sin(hh))
        assert oog < 2e-3, (pole, t, oog); out.append(rgb)
    return out
TAUMAP = ListedColormap(arm(NEG_HEX, gamma=1.2)[::-1] + arm(POS_HEX, gamma=0.6), name="tau_olive_indigo")
NORM = TwoSlopeNorm(vmin=VMIN, vcenter=0, vmax=VMAX)
def wcag(c1, c2):
    def Y(cc): r, g, b = (_s2l(v) for v in to_rgb(cc)); return 0.2126 * r + 0.7152 * g + 0.0722 * b
    a, b = sorted([Y(c1), Y(c2)], reverse=True); return (a + 0.05) / (b + 0.05)
CELL_CONTRAST = [wcag(INK, TAUMAP(NORM(v))) for v in T.values[SIG]]
assert min(CELL_CONTRAST) >= 4.5, min(CELL_CONTRAST)               # 숫자를 적는 칸 위의 INK 대비

# ================================================================== 그림 틀 (인치)
W, HA, GAP, HB = 7.2, 1.36, 0.06, 4.41      # a는 인치 좌표 y ∈ [0, HA]에 그린다 (b: 14줄, 줄 높이 0.236 in)
H = HA + GAP + HB
fig = plt.figure(figsize=(W, H))
FS, FSS, FSL = 7.0, 6.5, 9                  # 본문·단계 제목, 작은 글자, 패널 글자 (pt)
TAU = r"$\tau$*"                            # τ는 STIX Sans(Arial의 τ는 작은 대문자 T처럼 보인다), 별표는 Arial
def letter(x_in, y_in, s): fig.text(x_in / W, y_in / H, s, fontsize=FSL, fontweight="bold", va="top", ha="left", color=INK)

A = fig.add_axes([0, (H - HA) / H, 1, HA / H]); A.set_xlim(0, W); A.set_ylim(0, HA); A.axis("off")   # 자료 좌표 1 = 1 in
def T_(x, y, s, size=FSS, color=INK, **kw): return A.text(x, y, s, fontsize=size, color=color, **kw)
def header(x, num, s): return T_(x, HA - 0.02, f"{num}  {s}", size=FS, fontweight="bold", va="top", ha="left")
def arrow(x0, y0, x1, y1, color=MUTED, lw=0.7, style="-|>", ms=6, rad=0.0):
    A.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle=style, color=color, lw=lw, mutation_scale=ms, shrinkA=0, shrinkB=0,
                                                       connectionstyle=f"arc3,rad={rad}"))
def cells(x0, ytop, ncol, nrow, w, h, face="#e6e6e6", gap=0.012):
    for i in range(nrow):
        for j in range(ncol):
            A.add_patch(Rectangle((x0 + j * w + gap / 2, ytop - (i + 1) * h + gap / 2), w - gap, h - gap, fc=face, ec="none"))
def vdots(x, y, dy=0.022, color="#9a9a9a"): A.plot([x] * 3, [y + dy, y, y - dy], ls="none", marker="o", ms=0.9, mew=0, color=color)
LAYERS6 = [("CAGE", COL["CAGE"]), ("DNase", COL["DNASE"]), ("H3K27ac", COL["H3K27ac"]),     # 보충 그림 S2와 같은 층위 순서
           ("H3K4me1", COL_OTHER["H3K4me1"]), ("H3K4me3", COL_OTHER["H3K4me3"]), ("RNA-seq", COL_OTHER["RNA"])]
J_MY, I_CAGE = LIN12.index("myeloid"), 0      # 예로 든 쌍: myeloid × CAGE
REFC, CELLC = "#a3a3a3", "#e6e6e6"
YMID = HA - 0.80                              # 세 단계의 가운데 높이(단계 사이 화살표)

# ------------------------------------------------------------------ a1: 대립유전자마다의 예측 트랙과 allelic effect
X0, X1 = 0.70, 2.06; XV = (X0 + X1) / 2; SPAN = 4000.0          # 트랙 영역 = 변이 가운데 4 kb
def gx(u): return XV + u * (X1 - X0) / SPAN
header(0.20, "1", "Predicted allelic effects")
DX, YREF, YALT = 0.08, HA - 0.33, HA - 0.465
A.add_patch(FancyBboxPatch((XV - 0.04, YALT - 0.062), 0.08, YREF - YALT + 0.124, boxstyle="round,pad=0,rounding_size=0.02", fc="#e3e3e3", ec="none"))
for y, seq, lab in [(YREF, "CTGAGTCCA", "REF"), (YALT, "CTGAATCCA", "ALT")]:
    for k, ch in enumerate(seq):
        T_(XV + (k - 4) * DX, y, ch, size=FS, ha="center", va="center", fontweight="bold" if k == 4 else "normal")
    T_(XV - 4.7 * DX, y, "…", size=FS, ha="right", va="center"); T_(XV + 4.7 * DX, y, "…", size=FS, ha="left", va="center")
    T_(XV - 6.6 * DX, y, lab, size=FSS, color=MUTED, ha="right", va="center")
KX = XV - 6.6 * DX - 0.42    # 두 대립유전자 예측 신호의 선 표지: REF 회색, ALT 층위 색
A.plot([KX, KX + 0.17], [YREF, YREF], color=REFC, lw=1.1, solid_capstyle="butt")
for q, m in enumerate(["CAGE", "DNASE", "H3K27ac"]):
    A.plot([KX + q * 0.057, KX + (q + 1) * 0.057], [YALT, YALT], color=COL[m], lw=1.1, solid_capstyle="butt")
TH, TG = 0.125, 0.025                                               # 트랙 높이, 간격
tops = [HA - 0.615, HA - 0.615 - (TH + TG), HA - 0.615 - 2 * (TH + TG)]; tops.append(tops[2] - TH - 0.07)   # DNase, CAGE +, CAGE −, H3K27ac
SEQL, SEQR, YSB = XV - 4.7 * DX - 0.10, XV + 4.7 * DX + 0.10, YALT - 0.08
A.plot([SEQL, XV - 0.004, None, SEQR, XV + 0.004], [YSB, tops[0] + 0.012, None, YSB, tops[0] + 0.012], color=GRID, lw=0.5)   # 서열 = 변이 자리를 확대한 것
u = np.linspace(-2000, 2000, 801)
def gs(mu, sd): return np.exp(-0.5 * ((u - mu) / sd) ** 2)
SIGNALS = [("DNase", "DNASE", 250, 0.05 + 0.50 * gs(0, 120) + 0.16 * gs(-760, 170) + 0.12 * gs(820, 200) + 0.06 * gs(-1500, 150),
            0.05 + 0.92 * gs(0, 120) + 0.16 * gs(-760, 170) + 0.12 * gs(820, 200) + 0.06 * gs(-1500, 150)),
           ("CAGE +", "CAGE", 250, 0.02 + 0.90 * gs(55, 28) + 0.25 * gs(125, 22), 0.02 + 0.42 * gs(55, 28) + 0.12 * gs(125, 22)),
           ("CAGE −", "CAGE", 250, 0.02 + 0.38 * gs(-95, 30) + 0.12 * gs(-165, 25), 0.02 + 0.62 * gs(-95, 30) + 0.18 * gs(-165, 25)),
           ("H3K27ac", "H3K27ac", 1000, 0.08 + 0.42 * gs(-480, 300) + 0.50 * gs(470, 320) + 0.10 * gs(-1500, 250),
            0.08 + 0.64 * gs(-480, 300) + 0.76 * gs(470, 320) + 0.10 * gs(-1500, 250))]
A.add_patch(Rectangle((gx(-250), tops[2] - TH - 0.004), gx(250) - gx(-250), tops[0] - tops[2] + TH + 0.012, fc="#000000", alpha=0.055, ec="none"))   # 501-bp 창
A.add_patch(Rectangle((gx(-1000), tops[3] - TH - 0.004), gx(1000) - gx(-1000), TH + 0.012, fc="#000000", alpha=0.055, ec="none"))                # 2,001-bp 창
for (lab, m, half, ref, alt), top in zip(SIGNALS, tops):
    y0, sc = top - TH, 0.92 * TH / max(ref.max(), alt.max()); x = gx(u); win = np.abs(u) <= half
    A.plot([gx(-2000), gx(2000)], [y0, y0], color=GRID, lw=0.5)
    A.fill_between(x[win], y0, y0 + sc * ref[win], color=REFC, alpha=0.30, lw=0)
    A.fill_between(x[win], y0, y0 + sc * alt[win], color=COL[m], alpha=0.22, lw=0)
    A.plot(x, y0 + sc * ref, color=REFC, lw=0.8); A.plot(x, y0 + sc * alt, color=COL[m], lw=0.9)
    T_(X0 - 0.06, y0 + TH * 0.42, lab, size=FS, ha="right", va="center")
for ya, yb_ in [(tops[2] - TH - 0.006, tops[0] + 0.012), (tops[3] - TH - 0.006, tops[3] + 0.010)]:   # 변이 자리
    A.plot([XV, XV], [ya, yb_], color=INK, lw=0.5, alpha=0.6)
A1_BOTTOM = tops[3] - TH

# ------------------------------------------------------------------ a2: 조직 묶음 × 층위 쌍마다 SNP의 최대 효과 → 상위 5%
BX = 2.38
header(BX, "2", "Tissue group × layer annotations")
# 작은 막대: 한 SNP가 예로 든 쌍의 트랙마다 갖는 효과. 쌍의 점수는 그 최댓값
MB0, MBW, MBG, MBH = BX + 0.08, 0.052, 0.022, 0.30; MBY = HA - 1.00
vals = [0.30, 0.55, 0.22, 1.0, 0.40, 0.62, 0.27]
for k, v in enumerate(vals):
    A.add_patch(Rectangle((MB0 + k * (MBW + MBG), MBY), MBW, MBH * v, fc=COL["CAGE"], alpha=1.0 if v == 1.0 else 0.32, ec="none"))
A.plot([MB0 - 0.03, MB0 + len(vals) * (MBW + MBG) + 0.01], [MBY, MBY], color=MUTED, lw=0.6)
kmax = vals.index(1.0); xmx = MB0 + kmax * (MBW + MBG) + MBW / 2
T_(xmx, MBY + MBH + 0.025, "max", ha="center", va="bottom")
# 12 조직 묶음 × 6 층위 격자, 예로 든 칸(myeloid × CAGE)
GX0, GW, GH, GTOP = 3.14, 0.095, 0.105, HA - 0.44
cells(GX0, GTOP, 12, 6, GW, GH, face=CELLC)
for i, (lab, col) in enumerate(LAYERS6):
    yc = GTOP - (i + 0.5) * GH
    A.add_patch(Rectangle((GX0 + 12 * GW + 0.04, yc - 0.03), 0.06, 0.06, fc=col, ec="none"))
    T_(GX0 + 12 * GW + 0.13, yc, lab, ha="left", va="center")
GBOT = GTOP - 6 * GH
cx, cy = GX0 + J_MY * GW, GTOP - (I_CAGE + 1) * GH
A.add_patch(Rectangle((cx + 0.006, cy + 0.006), GW - 0.012, GH - 0.012, fc=COL["CAGE"], ec=INK, lw=0.9))
xb0, xb1, yb = GX0 + 0.006, GX0 + 12 * GW - 0.006, GBOT - 0.05
A.plot([xb0, xb0, xb1, xb1], [yb + 0.02, yb, yb, yb + 0.02], color=MUTED, lw=0.6)
T_((xb0 + xb1) / 2, yb - 0.015, "12 tissue groups", color=MUTED, ha="center", va="top")
arrow(xmx + 0.10, MBY + MBH + 0.075, cx + GW * 0.40, cy + GH + 0.006, color=INK, lw=0.6, ms=6, rad=-0.32)      # 최댓값 → 그 쌍의 점수
arrow(X1 + 0.07, YMID, MB0 - 0.09, YMID)                                                                       # 1 → 2
# 점수 순으로 줄 세운 SNP와 상위 5%
LABR = GX0 + 12 * GW + 0.13 + 0.42                                   # 층위 이름의 오른쪽 끝(대략; 아래에서 잰 값으로 확인)
CX0, CX1, CY0, CY1 = 5.10, 5.90, HA - 0.98, HA - 0.60
rr = np.linspace(0, 1, 400); s = 0.06 + 0.94 * np.exp(-(1 - rr) / 0.13); xr = CX0 + rr * (CX1 - CX0); yr = CY0 + s * (CY1 - CY0)   # SNP를 효과 순으로: 오른쪽이 큰 효과
top5, rest = rr >= 0.95, rr <= 0.95 + 1e-9                             # 오른쪽 끝 5% = 상위 5%
A.fill_between(xr[rest], CY0, yr[rest], color=CELLC, lw=0)
A.fill_between(xr[top5], CY0, yr[top5], color=COL["CAGE"], lw=0)
A.plot(xr, yr, color=MUTED, lw=0.7); A.plot([CX0, CX1], [CY0, CY0], color=MUTED, lw=0.6)
T_(CX1 - 0.11, CY1 - 0.005, "top 5%", ha="right", va="top")
T_((CX0 + CX1) / 2 - 0.05, CY0 - 0.035, "SNPs ranked by effect", color=MUTED, ha="center", va="top")
arrow(LABR + 0.05, YMID, CX0 - 0.06, YMID)                                                                     # 쌍의 점수 → 순위

# ------------------------------------------------------------------ a3: S-LDSC (baselineLD v2.2 + 그 쌍의 주석) → τ*
DX0 = 6.16
header(DX0, "3", "S-LDSC")
rng = np.random.default_rng(7)
SX0, SX1 = DX0 + 0.02, DX0 + 0.56
T_(SX0, HA - 0.34, "baselineLD v2.2", color=MUTED, ha="left", va="top")
ys = [HA - 0.56 - k * 0.06 for k in range(4)]; yn = CY0 + 0.018       # 기준선 주석 네 줄(… 줄임), 맨 아래 = 그 쌍의 주석
for y in ys:
    A.plot([SX0, SX1], [y, y], color=GRID, lw=0.5)
    for b0 in rng.uniform(0, SX1 - SX0 - 0.07, 3):
        A.add_patch(Rectangle((SX0 + b0, y - 0.013), rng.uniform(0.03, 0.07), 0.026, fc="#bdbdbd", ec="none"))
vdots((SX0 + SX1) / 2, (ys[-1] + yn) / 2 + 0.003)
A.plot([SX0, SX1], [yn, yn], color=GRID, lw=0.5)
for f in (0.15, 0.675):
    A.add_patch(Rectangle((SX0 + f * (SX1 - SX0), yn - 0.018), 0.035, 0.036, fc=COL["CAGE"], ec="none"))
arrow(CX1 + 0.05, yn, SX0 - 0.06, yn)                                                                          # 상위 5% → 이진 주석 (2 → 3)
YM3 = (ys[0] + yn) / 2
A.plot([SX1 + 0.05, SX1 + 0.07, SX1 + 0.07, SX1 + 0.05], [ys[0] + 0.02, ys[0] + 0.02, yn - 0.02, yn - 0.02], color=MUTED, lw=0.6)
arrow(SX1 + 0.07, YM3, SX1 + 0.22, YM3)
T_(SX1 + 0.26, YM3, TAU, size=FS, ha="left", va="center", math_fontfamily="stixsans")

# ================================================================== b: 조직 묶음 지도
letter(0.02, H - 0.02, "a"); letter(0.02, HB + 0.02, "b")
HX0, HW, HY0, HH = 2.08, 4.40, 0.98, 2.60 * 14 / 11      # 열지도 상자 (인치)
ax = fig.add_axes([HX0 / W, HY0 / H, HW / W, HH / H])
nr, ncl = T.shape
ax.imshow(T.values, cmap=TAUMAP, norm=NORM, aspect="auto", extent=(-0.5, ncl - 0.5, nr - 0.5, -0.5), interpolation="nearest")
for x in np.arange(0.5, ncl - 0.5): ax.axvline(x, color="white", lw=0.7)
for y in np.arange(0.5, nr - 0.5): ax.axhline(y, color="white", lw=0.7)
GROUP_BREAKS = [4.5, 7.5, 9.5, 11.5]
for y in GROUP_BREAKS: ax.axhline(y, color="white", lw=2.6)
SHOWN = []
for i in range(nr):
    for j in range(ncl):
        if SIG[i, j]:
            lab = f"{T.iat[i, j]:+.1f}"; SHOWN.append((ORD[i], LIN12[j], lab))
            ax.text(j, i, lab, ha="center", va="center", fontsize=FS, fontweight="bold", color=INK)
ax.set_xlim(-0.5, ncl - 0.5); ax.set_ylim(nr - 0.5, -0.5)
ax.set_xticks(range(ncl)); ax.set_xticklabels([LN[l] for l in LIN12], rotation=40, ha="right", rotation_mode="anchor", fontsize=FS, color=INK)
ax.set_yticks(range(nr)); ax.set_yticklabels([TN[t] for t in ORD], fontsize=FS, color=INK)
ax.tick_params(length=0, pad=3); [sp.set_visible(False) for sp in ax.spines.values()]
# 형질 범주의 조직 묶음: 띠(색 없이)와 묶음 이름 (형질 이름 왼쪽)
SX, SW = 0.70, 0.07
hs = fig.add_axes([SX / W, HY0 / H, SW / W, HH / H]); hs.set_xlim(0, 1); hs.set_ylim(nr - 0.5, -0.5); hs.axis("off")
cgroups = []
for i, t in enumerate(ORD):
    if not cgroups or cgroups[-1][0] != OWN[t]: cgroups.append([OWN[t], i, i])
    else: cgroups[-1][2] = i
assert [g[0] for g in cgroups] == ["CNS_tissue", "myeloid", "lymphoid", "liver", "stromal"] and [y - .5 for _, y, _ in cgroups[1:]] == GROUP_BREAKS
STRIP = "#8a8a8a"
for g, i0_, i1_ in cgroups:
    hs.add_patch(Rectangle((0, i0_ - 0.5 + 0.06), 1, i1_ - i0_ + 1 - 0.12, fc=STRIP, ec="none"))
    fig.text((SX - 0.06) / W, (HY0 + HH * (1 - ((i0_ + i1_) / 2 + 0.5) / nr)) / H, LN[g].replace("/", "/\n"), fontsize=FS, color=INK, ha="right", va="center", linespacing=1.1)
# 색 막대(τ*에 선형)와 숫자 표기 기준
CBX, CBW, CBH = 6.70, 0.085, 1.55
cb_ax = fig.add_axes([CBX / W, (HY0 + HH - CBH) / H, CBW / W, CBH / H])
cbv = np.linspace(VMIN, VMAX, 512)
cb_ax.imshow(TAUMAP(NORM(cbv))[::-1, None, :3], aspect="auto", extent=(0, 1, VMIN, VMAX))
CB_TICKS = list(range(-2, 5)); assert VMIN < CB_TICKS[0] and CB_TICKS[-1] < VMAX
cb_ax.set_xticks([]); cb_ax.yaxis.tick_right(); cb_ax.set_yticks(CB_TICKS); cb_ax.set_ylim(VMIN, VMAX)
cb_ax.set_yticklabels([f"{v:+d}".replace("-", "−") if v else "0" for v in CB_TICKS], fontsize=FS, color=INK)
cb_ax.tick_params(axis="y", length=2, width=0.5, color=MUTED, pad=2); [sp.set_visible(False) for sp in cb_ax.spines.values()]
fig.text((CBX + CBW / 2) / W, (HY0 + HH + 0.05) / H, TAU, fontsize=FS, color=INK, ha="center", va="bottom", math_fontfamily="stixsans")
# (열지도 숫자의 기준 'FDR < 0.05'는 그림 설명에 쓴다: 그림 안에는 결과만 둔다)

# ================================================================== 겹침·잘림 점검 후 저장
def _box(t, rd):
    """글자의 상자 네 꼭짓점(화면 좌표). 돌린 글자는 돌린 사각형으로 잰다(곧은 상자로 재면 나란한 기울인 눈금 글자가 겹친 것으로 나온다)."""
    rot = t.get_rotation() % 360
    if rot % 90 == 0:
        e = matplotlib.text.Text.get_window_extent(t, rd); return np.array([[e.x0, e.y0], [e.x1, e.y0], [e.x1, e.y1], [e.x0, e.y1]])
    assert t.get_rotation_mode() == "anchor" or (t.get_ha() == "center" and t.get_va() == "center"), ("돌린 글자의 상자를 잴 수 없다", t.get_text())
    t.set_rotation(0); e = matplotlib.text.Text.get_window_extent(t, rd); t.set_rotation(rot)
    a = np.asarray(t.get_transform().transform(t.get_unitless_position()), float); th = math.radians(rot)
    M = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
    return (np.array([[e.x0, e.y0], [e.x1, e.y0], [e.x1, e.y1], [e.x0, e.y1]]) - a) @ M.T + a
def _hit(P, Q_):
    for poly in (P, Q_):
        for i in range(len(poly)):
            e = poly[(i + 1) % len(poly)] - poly[i]; nrm = np.array([-e[1], e[0]]); a, b = P @ nrm, Q_ @ nrm
            if a.max() <= b.min() + 1e-6 or b.max() <= a.min() + 1e-6: return False
    return True
def text_issues(fig):
    fig.canvas.draw(); rd = fig.canvas.get_renderer(); fb = fig.bbox
    ts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    for t in ts:
        s = t.get_text(); assert "τ" not in s and "\t" not in s, repr(s)
        if r"\tau" in s: assert t.get_math_fontfamily() == "stixsans", s
        assert "lineage" not in s.lower() and "axis" not in s.lower() and "axes" not in s.lower(), s
        assert 6.5 <= t.get_fontsize() <= 7 or (t.get_fontsize() == FSL and t.get_fontweight() == "bold"), (s, t.get_fontsize())
        assert matplotlib.colors.to_hex(t.get_color()) in (INK, MUTED), (s, t.get_color())                 # 글자 색은 INK/MUTED만
    bb = [(t, _box(t, rd)) for t in ts]
    ov = [(a.get_text(), b.get_text()) for i, (a, x) in enumerate(bb) for b, y in bb[i + 1:] if _hit(x, y)]
    cl = [t.get_text() for t, x in bb if x[:, 0].min() < fb.x0 - 0.5 or x[:, 0].max() > fb.x1 - 2 or x[:, 1].min() < fb.y0 + 2 or x[:, 1].max() > fb.y1 + 0.5]
    return ov, cl, bb
ov, cl, bb = text_issues(fig)
_lab = [x for t, x in bb if t.axes is A and t.get_text() in {l for l, _ in LAYERS6}]
assert abs(max(x[:, 0].max() for x in _lab) / fig.dpi - LABR) < 0.03, (max(x[:, 0].max() for x in _lab) / fig.dpi, LABR)   # 층위 이름의 오른쪽 끝
print("text overlaps:", ov); print("clipped:", cl)
assert not ov and not cl
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"fig1.{ext}", dpi=300)
plt.close(fig)
NUM = {
    "source": "e22_benchmark14/results/w04_tissue.csv: tau, q (BH over the 168 cells, one-sided tau* > 0); model M1 = baselineLD v2.2 + one tissue-group annotation",
    "b_rows": [dict(key=t, name=TN[t], category_tissue_group=OWN[t], top_ranked_tissue_group=RANK1[t], n_enriched_groups=N_ENR[t]) for t in ORD],
    "b_row_blocks": [dict(tissue_group=g, name=LN[g], traits=ORD[a:b + 1]) for g, a, b in cgroups],
    "b_tau_star": {t: {l: round(tau(t, l), 3) for l in LIN12} for t in ORD},
    "b_q_bh": {t: {l: float(f"{Q.loc[t, l]:.4g}") for l in LIN12} for t in ORD},
    "b_n_cells": int(T.size), "b_n_significant": int(SIG.sum()), "b_numbers_shown": [f"{t}|{l}|{s_}" for t, l, s_ in SHOWN],
    "b_n_traits_category_group_ranked_first": int(N_OWN_FIRST), "b_traits_category_group_not_first": {t: dict(top=RANK1[t], rank_of_category_group=int((T.loc[t] > tau(t, OWN[t])).sum()) + 1) for t in ORD if RANK1[t] != OWN[t]},
    "b_colour_scale": {"vmin": VMIN, "vcenter": 0, "vmax": VMAX, "colour_bar_ticks": CB_TICKS, "pos_pole": POS_HEX, "neg_pole": NEG_HEX, "mid": MID_HEX,
                       "min_INK_contrast_on_shown_cells": round(min(CELL_CONTRAST), 2)},
    "b_tau_range": [round(float(T.values.min()), 3), round(float(T.values.max()), 3)],
    "b_strip": "neutral blocks (no colour coding) labelled with the tissue group of each trait's category",
    "a_schematic_steps": ["Predicted allelic effects", "Tissue group × layer annotations", "S-LDSC"],
    "a_step2_mini_curve": "SNPs ranked by effect, increasing to the right; the right-hand 5% (the top 5%) is shaded",
    "a_legend_constants_not_drawn": {"common_SNPs": N_SNP, "common_SNPs_label": f"{N_SNP / 1e6:.2f} million", "tracks": N_TRACK, "biosamples": N_BIOS,
                                     "tissue_groups": 12, "layers": 6, "top_fraction": "5%", "window_bp_DNase_CAGE": 501, "window_bp_histone": 2001,
                                     "baselineLD_v2.2_annotations": N_BASE, "track_counts_tissue_group_x_layer": {r[0]: [int(x) for x in r[1:7]] for r in rows}},
    "figure_size_in": [W, H],
}
(OUT / "numbers_fig1.json").write_text(json.dumps(NUM, indent=1, default=float, ensure_ascii=False), encoding="utf-8")
print("fig1 written:", OUT / "fig1.png", f"({W} × {H:.2f} in) | significant cells:", int(SIG.sum()), "| min INK contrast on shown cells:", round(min(CELL_CONTRAST), 2))
