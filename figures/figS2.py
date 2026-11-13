# -*- coding: utf-8 -*-
"""Supplementary Figure S2 — 106개 GWAS의 조직 묶음 지도.
 열지도: 행 = GWAS 106개, 열 = 조직 묶음 12개(그림 1b와 같은 순서), 색 = τ*(그림 1b와 같은 색 지도와 색 범위 −2.5 ~ +4.5),
   점 = FDR < 0.05인 칸(전체 1,272칸의 Benjamini–Hochberg, 단측 τ* > 0). 모형은 그림 1b와 같다: baselineLD v2.2 + 조직 묶음 주석 하나.
 행은 형질 범주(골격, 혈액·면역, 지질, 뇌·행동, 그 밖)로 묶고 가는 옆 띠와 범주 이름으로 표시한다(기술용).
   범주 안에서는 τ*가 가장 큰 조직 묶음(열 순서) → 그 τ*가 큰 순서로 놓는다(u02_analysis_lineage.csv의 top_lineage, top_tau; 정렬일 뿐 새 통계가 아니다).
   행이 많아 두 덩어리로 나란히 놓는다: 왼쪽 = 골격·혈액·면역·지질·뇌·행동(50행), 오른쪽 = 그 밖(56행).
 행 이름: figures/panel_names.json(GWAS 파일 이름 → 형질 이름). 같은 이름이 둘이면(류마티스 관절염, 전립선암) 연구 이름을 괄호로 붙인다. 칸보다 긴 이름은 낱말 경계에서 자르고 "…"를 붙인다.
 자료: e20_gwas_panel/results/u02_lineage.csv, u02_analysis_lineage.csv; 유의한 칸의 수는 u02_lineage.txt와 대조한다(다르면 멈춘다). 새 통계 없음.
실행: python figures/figS2.py      (저장소 뿌리에서; 모두 한 번에는 python figures/make_figures.py)
출력: figures/out/figS2.png, figS2.pdf, numbers_figS2.json
"""
import sys, json, math, re
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))      # -I 로 실행하면 스크립트 폴더가 sys.path에 없다
import numpy as np, pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, to_hex
from matplotlib.patches import Rectangle, Polygon
from matplotlib.lines import Line2D
from common import R as RESROOT, CATCOL, INK, MUTED, LIN12, LN
from taumap import TAUMAP, wcag

plt.rcParams.update({"font.family": "Arial", "font.size": 7, "text.color": INK, "axes.unicode_minus": True, "mathtext.fontset": "custom",
                     "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "pdf.fonttype": 42})
from matplotlib import _mathtext as _mt
def _xheight_from_glyph(self, fontname, fontsize, dpi):
    """Arial의 PCLT x-height를 약 1.9배 크게 잡아 위·아래첨자가 어긋나므로 'x' 글자 높이를 잰다."""
    return self.get_metrics(fontname, matplotlib.rcParams["mathtext.default"], "x", fontsize, dpi).iceberg
_mt.TruetypeFonts.get_xheight = _xheight_from_glyph
HERE = Path(__file__).resolve().parent
OUT = HERE / "out"; OUT.mkdir(parents=True, exist_ok=True)
FS, FSS, FS_ROW = 7.0, 6.5, 5.5                                # 본문, 작은 글자, 행 이름(행이 106개라 5.5 pt)
TAU = r"$\tau$*"

# ================================================================== 자료
U_DIR = RESROOT / "e20_gwas_panel" / "results"
U2 = pd.read_csv(U_DIR / "u02_lineage.csv", encoding="utf-8-sig")
UA = pd.read_csv(U_DIR / "u02_analysis_lineage.csv", encoding="utf-8-sig").set_index("id")
assert len(U2) == 1272 and (U2.fit == "ok").all() and not U2.duplicated(["id", "lineage"]).any() and set(U2.lineage) == set(LIN12)
T = U2.pivot(index="id", columns="lineage", values="tau").reindex(columns=LIN12)
Q = U2.pivot(index="id", columns="lineage", values="q").reindex(columns=LIN12)
assert T.shape == (106, 12) and not T.isna().any().any() and not Q.isna().any().any() and set(T.index) == set(UA.index)
META = U2.drop_duplicates("id").set_index("id")[["file", "group", "status"]]
assert (META.loc[UA.index, "file"] == UA.file).all() and (META.loc[UA.index, "group"] == UA.group).all()
assert (T.idxmax(axis=1).reindex(UA.index) == UA.top_lineage).all()                       # top_lineage = τ*가 가장 큰 조직 묶음
SIG = Q < .05
N_CELL, N_SIG, N_TR = int(SIG.size), int(SIG.values.sum()), int(SIG.any(axis=1).sum())
assert (T.values[SIG.values] > 0).all()                                                    # 단측 검정(τ* > 0): 점이 있는 칸은 모두 양수
# 개수 대조: u02_lineage.py가 적은 요약(u02_lineage.txt)의 값과 같아야 한다(다르면 멈춘다)
_m = re.search(r"FDR < 0\.05인 칸: (\d+) / (\d+)", (U_DIR / "u02_lineage.txt").read_text(encoding="utf-8")); assert _m, "u02_lineage.txt의 개수 줄을 찾지 못했다"
REP_SIG = int(_m[1])
assert (N_CELL, N_SIG) == (1272, REP_SIG) and int(_m[2]) == 1272, f"개수가 u02_lineage.txt와 다르다: CSV {N_SIG}칸, 요약 {REP_SIG}칸"
GSIZE = META.group.value_counts().to_dict(); assert GSIZE == {"U": 56, "BR": 25, "BI": 17, "S": 5, "LL": 3}, GSIZE
N_SAME14 = int((META.status == "discovery").sum()); assert set(META.status) == {"new", "discovery"} and N_SAME14 == 13

# ------------------------------------------------------------------ 이름: figures/panel_names.json (GWAS 파일 이름 → 형질 이름, 연구 이름의 규칙, 형질 범주의 이름)
_names = json.loads((HERE / "panel_names.json").read_text(encoding="utf-8"))
NAME, STUDY_FIXED, STUDY_SOURCE, GROUPNAME = _names["NAME"], _names["STUDY_FIXED"], _names["STUDY_SOURCE"], _names["GROUP"]
GORDER = ["S", "BI", "LL", "BR", "U"]; assert list(GROUPNAME) == GORDER and set(NAME) == set(META.file) and GROUPNAME["U"] == "Other"
def study_of(file):
    """연구 이름: 출처 접두(UK Biobank 등)가 있으면 그것, PASS 파일이면 저자·연도 접미."""
    src, _, rest = file.partition(".")
    if src in STUDY_SOURCE: return STUDY_SOURCE[src]
    assert src == "PASS" and "." in rest, file
    tok = rest.rsplit(".", 1)[1]
    if tok in STUDY_FIXED: return STUDY_FIXED[tok]
    m_ = re.fullmatch(r"([A-Z][a-z]+)(\d{4})", tok); assert m_, file
    return f"{m_.group(1)} {m_.group(2)}"
_dup = {n for n in NAME.values() if list(NAME.values()).count(n) > 1}; assert _dup == {"Rheumatoid arthritis", "Prostate cancer"}, _dup
def full_label(file):
    return f"{NAME[file]} ({study_of(file)})" if NAME[file] in _dup else NAME[file]
FULL = {i: full_label(META.loc[i, "file"]) for i in T.index}
assert len(set(FULL.values())) == 106

# ------------------------------------------------------------------ 행 순서: 형질 범주 → τ*가 가장 큰 조직 묶음(열 순서) → 그 τ*가 큰 순서
ROWS = {g: sorted(META.index[META.group == g], key=lambda i: (LIN12.index(UA.loc[i, "top_lineage"]), -UA.loc[i, "top_tau"], i)) for g in GORDER}
BLOCKS = [["S", "BI", "LL", "BR"], ["U"]]                         # 왼쪽 덩어리, 오른쪽 덩어리
assert sum(len(ROWS[g]) for g in GORDER) == 106 and [sum(len(ROWS[g]) for g in b) for b in BLOCKS] == [50, 56]

# ================================================================== 색 (그림 1b와 같은 색 지도·색 범위)
VMIN, VMAX = -2.5, 4.5
NORM = TwoSlopeNorm(vmin=VMIN, vcenter=0, vmax=VMAX)
SAT_LO = [(i, l, round(float(T.loc[i, l]), 2)) for i in T.index for l in LIN12 if T.loc[i, l] < VMIN]
SAT_HI = [(i, l, round(float(T.loc[i, l]), 2)) for i in T.index for l in LIN12 if T.loc[i, l] > VMAX]
assert len(SAT_LO) == 3 and len(SAT_HI) == 1 and not SIG.values[T.values < VMIN].any(), (SAT_LO, SAT_HI)   # 색 범위를 넘는 칸: 아래 3, 위 1
RING = "white"
DOT_CONTRAST = min(wcag(INK, TAUMAP(NORM(v))) for v in T.values[SIG.values]); assert DOT_CONTRAST >= 3.0, DOT_CONTRAST   # 점(INK, 흰 테) 대 칸: 3:1 이상
GCOL = {"S": CATCOL["Skeletal"], "BI": CATCOL["Blood & immune"], "LL": CATCOL["Lipid"], "BR": CATCOL["Brain"], "U": "#8f8f8f"}   # 옆 띠: 형질 범주의 색, 그 밖(Other)은 회색

# ================================================================== 배치 (인치)
W = 7.2
RP, GG = 0.118, 0.07                                              # 행 간격, 묶음 사이 간격
GLW, SW, G1, G2 = 0.13, 0.05, 0.05, 0.04                          # 묶음 이름(세로 글자) 칸, 옆 띠, 띠–행 이름 사이, 행 이름–열지도 사이
MLEFT, MRIGHT, BGAP = 0.05, 0.06, 0.16                            # 왼쪽·오른쪽 여백, 두 덩어리 사이
BW = (W - MLEFT - MRIGHT - BGAP) / 2                              # 덩어리 하나의 폭
HW = 1.34; CWI = HW / 12                                          # 열지도 폭, 열 하나의 폭
LW = BW - GLW - SW - G1 - G2 - HW                                 # 행 이름 칸
BH = [sum(len(ROWS[g]) for g in b) * RP + (len(b) - 1) * GG for b in BLOCKS]
fig = plt.figure(figsize=(W, 8), dpi=300); RD = fig.canvas.get_renderer()                  # 글자 폭은 저장 해상도(300 dpi)에서 잰다
def tsize(s, fs, **kw):
    t = fig.text(0, 0, s, fontsize=fs, **kw); e = t.get_window_extent(RD); t.remove(); return e.width / fig.dpi, e.height / fig.dpi
TOPLAB = max(tsize(LN[l], FSS)[0] for l in LIN12)                 # 위쪽 열 이름(세로 글자)의 높이
MB, MT = 0.10, 0.05
H = MB + max(BH) + 0.06 + TOPLAB + MT
fig.set_size_inches(W, H)
YTOP = MB + max(BH)                                               # 두 덩어리의 위쪽 끝(같은 높이에서 시작)
def axin(l, b, w, h): return fig.add_axes([l / W, b / H, w / W, h / H])
def fit_label(s):
    """행 이름이 칸(LW)보다 길면 낱말 경계에서 자르고 '…'를 붙인다."""
    if tsize(s, FS_ROW)[0] <= LW: return s
    words = s.split(" ")
    for k in range(len(words) - 1, 0, -1):
        cut = " ".join(words[:k]).rstrip(",;") + " …"
        if tsize(cut, FS_ROW)[0] <= LW: return cut
    raise AssertionError(("행 이름을 칸에 넣을 수 없다", s))
SHOWN = {i: fit_label(FULL[i]) for i in T.index}
TRUNC = {FULL[i]: SHOWN[i] for i in T.index if SHOWN[i] != FULL[i]}

AXES, n_dots = {}, 0
for k, block in enumerate(BLOCKS):
    x0 = MLEFT + k * (BW + BGAP); xs = x0 + GLW; xl = xs + SW + G1; xh = xl + LW + G2     # 띠, 행 이름, 열지도의 왼쪽 끝
    y = YTOP
    for j, g in enumerate(block):
        ids = ROWS[g]; n = len(ids); h = n * RP; y -= h
        ax = axin(xh, y, HW, h); AXES[g] = ax
        ax.imshow(T.loc[ids].values, cmap=TAUMAP, norm=NORM, aspect="auto", extent=(-0.5, 11.5, n - 0.5, -0.5), interpolation="nearest")
        for x in np.arange(0.5, 11.5): ax.axvline(x, color="white", lw=0.45)
        for yy in np.arange(0.5, n - 0.5): ax.axhline(yy, color="white", lw=0.45)
        ii, jj = np.where(SIG.loc[ids].values); n_dots += len(ii)
        ax.plot(jj, ii, ls="none", marker="o", ms=2.3, mfc=INK, mec=RING, mew=.35, zorder=4)
        ax.set_xlim(-0.5, 11.5); ax.set_ylim(n - 0.5, -0.5)
        ax.set_yticks(range(n)); ax.set_yticklabels([SHOWN[i] for i in ids], fontsize=FS_ROW, color=INK)
        ax.tick_params(length=0, pad=G2 * 72); [sp.set_visible(False) for sp in ax.spines.values()]
        if j == 0:                                               # 열 이름은 덩어리마다 맨 위에(세로 글자)
            ax.xaxis.tick_top(); ax.set_xticks(range(12)); ax.set_xticklabels([LN[l] for l in LIN12], rotation=90, ha="center", va="bottom", fontsize=FSS, color=INK)
            ax.tick_params(axis="x", pad=2)
        else: ax.set_xticks([])
        fig.add_artist(Rectangle((xs / W, (y + .012) / H), SW / W, (h - .024) / H, fc=GCOL[g], ec="none", transform=fig.transFigure))   # 옆 띠
        fig.text((xs - .035) / W, (y + h / 2) / H, GROUPNAME[g], rotation=90, ha="right", va="center", fontsize=FSS, color=INK)        # 묶음 이름
        assert tsize(GROUPNAME[g], FSS)[0] <= h, (g, "묶음 이름이 띠보다 길다")
        y -= GG
assert n_dots == N_SIG                                           # FDR < 0.05인 칸마다 점 하나

# 색 막대(가로)와 점 풀이: 왼쪽 덩어리 아래(오른쪽 덩어리보다 짧아 남는 자리)
XH1 = MLEFT + GLW + SW + G1 + LW + G2
CBY = YTOP - BH[0] - 0.20; CBH = 0.07; assert BH[1] - BH[0] > 0.45
cb = axin(XH1, CBY, HW, CBH)
cb.imshow(TAUMAP(NORM(np.linspace(VMIN, VMAX, 512)))[None, :, :3], aspect="auto", extent=(VMIN, VMAX, 0, 1))
CB_TICKS = list(range(-2, 5))
cb.set_yticks([]); cb.set_xticks(CB_TICKS); cb.set_xlim(VMIN, VMAX)
cb.set_xticklabels([f"{v:+d}".replace("-", "−") if v else "0" for v in CB_TICKS], fontsize=FSS, color=INK)
cb.tick_params(axis="x", length=2, width=0.5, color=MUTED, pad=1.5); [sp.set_visible(False) for sp in cb.spines.values()]
d = .035 * (VMAX - VMIN)                                          # 색 범위를 넘는 값이 있다는 표시(양쪽 끝의 세모)
cb.add_patch(Polygon([(VMIN, 0), (VMIN, 1), (VMIN - d, .5)], closed=True, fc=TAUMAP(0.0), ec="none", clip_on=False))
cb.add_patch(Polygon([(VMAX, 0), (VMAX, 1), (VMAX + d, .5)], closed=True, fc=TAUMAP(1.0), ec="none", clip_on=False))
fig.text((XH1 - .10) / W, (CBY + CBH / 2) / H, TAU, fontsize=FS, color=INK, ha="right", va="center", math_fontfamily="stixsans")
KX = XH1 - 1.00                                                   # 점 풀이
fig.add_artist(Line2D([KX / W], [(CBY + CBH / 2) / H], ls="none", marker="o", ms=2.3, mfc=INK, mec=RING, mew=.35, transform=fig.transFigure))
fig.text((KX + .06) / W, (CBY + CBH / 2) / H, "FDR < 0.05", fontsize=FSS, color=MUTED, ha="left", va="center")

# ================================================================== 겹침·잘림 점검 후 저장
def _box(t, rd):
    rot = t.get_rotation() % 360; assert rot % 90 == 0, t.get_text()
    e = matplotlib.text.Text.get_window_extent(t, rd); return e
def text_issues(fig):
    fig.canvas.draw(); rd = fig.canvas.get_renderer(); fb = fig.bbox
    ts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    rowlab = {id(lb) for ax in AXES.values() for lb in ax.get_yticklabels()}
    for t in ts:
        s = t.get_text(); assert "τ" not in s and "\t" not in s, repr(s)
        if r"\tau" in s: assert t.get_math_fontfamily() == "stixsans", s
        assert not any(w in s.lower() for w in ("lineage", "axis", "confirmed")), s
        assert to_hex(t.get_color()) in (INK, MUTED), (s, t.get_color())
        assert (t.get_fontsize() == FS_ROW) if id(t) in rowlab else (6.5 <= t.get_fontsize() <= 7), (s, t.get_fontsize())   # 행 이름 5.5 pt, 나머지 6.5–7 pt
    bb = [(t, _box(t, rd)) for t in ts]
    ov = [(a.get_text(), b.get_text()) for i, (a, x) in enumerate(bb) for b, y in bb[i + 1:] if x.overlaps(y)]
    cl = [t.get_text() for t, x in bb if x.x0 < fb.x0 + 1 or x.x1 > fb.x1 - 1 or x.y0 < fb.y0 + 1 or x.y1 > fb.y1 - 1]
    wide = [t.get_text() for t, x in bb if id(t) in rowlab and x.width / fig.dpi > LW + 1e-3]                              # 행 이름이 칸 안에 드는지
    return ov, cl, wide
ov, cl, wide = text_issues(fig)
print("text overlaps:", ov[:10]); print("clipped:", cl); print("row labels wider than their column:", wide); print("truncated labels:", TRUNC)
assert not ov and not cl and not wide
fig.savefig(OUT / "figS2.png", dpi=300); fig.savefig(OUT / "figS2.pdf")
plt.close(fig)
NUM = {
    "source": {"tau_fdr": "e20_gwas_panel/results/u02_lineage.csv", "top_tissue_group": "e20_gwas_panel/results/u02_analysis_lineage.csv",
               "counts_checked_against": "e20_gwas_panel/results/u02_lineage.txt", "names": "figures/panel_names.json"},
    "model": "baselineLD v2.2 + one tissue-group annotation (e20_gwas_panel u01_fit.py, M1_<tissue group>; annotations of e02_lineage_map), as in Fig. 1b",
    "counts": {"gwas": 106, "tissue_groups": 12, "cells": N_CELL, "cells_fdr_lt_0.05": N_SIG, "gwas_with_at_least_one": N_TR,
               "reported_in_u02_lineage_txt": {"cells_fdr_lt_0.05": REP_SIG}, "same_trait_as_one_of_the_14": N_SAME14,
               "group_sizes": {GROUPNAME[g]: len(ROWS[g]) for g in GORDER},
               "cells_fdr_lt_0.05_by_tissue_group": {LN[l]: int(SIG[l].sum()) for l in LIN12}},
    "colour_scale": {"vmin": VMIN, "vcenter": 0, "vmax": VMAX, "same_as": "Fig. 1b", "colour_bar_ticks": CB_TICKS,
                     "below_vmin": [dict(gwas=FULL[i], tissue_group=LN[l], tau_star=v) for i, l, v in SAT_LO],
                     "above_vmax": [dict(gwas=FULL[i], tissue_group=LN[l], tau_star=v) for i, l, v in SAT_HI],
                     "tau_range": [round(float(T.values.min()), 3), round(float(T.values.max()), 3)], "min_INK_dot_vs_cell": round(DOT_CONTRAST, 2)},
    "columns": [LN[l] for l in LIN12],
    "row_order_rule": "trait category (S, BI, LL, BR, U); within a group by the tissue group with the largest tau* (column order), then by decreasing tau* in it",
    "blocks": [[GROUPNAME[g] for g in b] for b in BLOCKS],
    "truncated_labels": TRUNC,
    "rows": [dict(id=i, file=META.loc[i, "file"], label=SHOWN[i], full_name=FULL[i], study=study_of(META.loc[i, "file"]), group=GROUPNAME[g], status=META.loc[i, "status"],
                  top_tissue_group=LN[UA.loc[i, "top_lineage"]], tau_star={LN[l]: round(float(T.loc[i, l]), 3) for l in LIN12},
                  fdr_lt_0_05=[LN[l] for l in LIN12 if SIG.loc[i, l]]) for g in GORDER for i in ROWS[g]],
    "figure_size_in": [W, H],
}
(OUT / "numbers_figS2.json").write_text(json.dumps(NUM, indent=1, default=float, ensure_ascii=False), encoding="utf-8")
print(f"cells {N_CELL}, FDR < 0.05: {N_SIG} (u02_lineage.txt {REP_SIG}); GWAS with at least one: {N_TR}; below {VMIN}: {SAT_LO}; above {VMAX}: {SAT_HI}")
print(f"label column {LW:.2f} in, cell {CWI:.3f} × {RP:.3f} in; saved {OUT / 'figS2.png'} ({W} × {H:.2f} in)")
