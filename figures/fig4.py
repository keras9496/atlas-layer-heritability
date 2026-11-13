# -*- coding: utf-8 -*-
"""Figure 4 — 106개 GWAS의 농축된 GWAS–조직 묶음 쌍의 CAGE–DNase 균형.
 균형 = atan2(DNase τ*, CAGE τ*) (도): 0° = CAGE만(전사 시작), 90° = DNase만(접근성), 45° = 같음.
 a 균형이 정의되는 쌍 285개(전체 트랙 주석). 줄 = 조직 묶음 주석 12개(위에서 아래로 b의 가중 모형 평균이 커지는 순서), 가로 = 균형, 점 하나 = 쌍 하나, 점 크기 = 균형의 s.e. 세 단계.
   한 줄 안의 세로 어긋남은 점이 겹치지 않게 놓은 것일 뿐 값이 아니다.
 b 같은 12줄: 가중 모형의 평균(원) ± 95% CI(선; jackknife s.e.), 쌍들의 중앙값(짧은 세로선), 오른쪽에 쌍 수. (w10_means.csv, set = twelve)
 c 맞춘 주석이 있는 다섯 주석: 전체 트랙, 층위당 세 단위, 세포형 일치 주석에서의 가중 모형 평균 ± 95% CI. (set = five; all, res, mat)
 d 같은 다섯 주석: TSS에서 1 kb 밖(원위)과 안(근위)의 SNP. (set = five; dist, prox)
 자료: e22_benchmark14/results/w05_pairs.csv, w08_pairs.csv, w10_means.csv, w10_contrasts.csv. 기존 결과 파일만 읽는다. 새 통계 없음.
실행: python figures/fig4.py      (저장소 뿌리에서; 모두 한 번에는 python figures/make_figures.py)
출력: figures/out/fig4.png, fig4.pdf, numbers_fig4.json
"""
import sys, json, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))      # -I 로 실행하면 스크립트 폴더가 sys.path에 없다
import numpy as np, pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from common import R as RESROOT, COL, LN, INK, MUTED, GRID

plt.rcParams.update({"font.family": "Arial", "font.size": 7, "text.color": INK,
                     "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
                     "pdf.fonttype": 42, "axes.unicode_minus": True})
from matplotlib import _mathtext as _mt
def _xheight_from_glyph(self, fontname, fontsize, dpi):
    return self.get_metrics(fontname, matplotlib.rcParams["mathtext.default"], "x", fontsize, dpi).iceberg
_mt.TruetypeFonts.get_xheight = _xheight_from_glyph
FS_NOTE, FS_LET, Z = 6.5, 9, 1.96
PT = 1 / 72
BAL = "CAGE–DNase balance"                                   # 그림 안에서 균형을 부르는 말(정의는 그림 설명에)
NEUTRAL, BAND = "#555555", "#f3f3f1"

# ------------------------------------------------------------------ 자료
RES = RESROOT / "e22_benchmark14" / "results"
OUT = Path(__file__).resolve().parent / "out"; OUT.mkdir(parents=True, exist_ok=True)
USED = {}
def src(name): p = RES / name; assert p.exists(), p; USED[name] = "e22_benchmark14/results/" + name; return p
PR = pd.read_csv(src("w05_pairs.csv")); PD_ = pd.read_csv(src("w08_pairs.csv")); MEANS = pd.read_csv(src("w10_means.csv")); CON = pd.read_csv(src("w10_contrasts.csv"))
assert len(PR) == 288 and PR.id.nunique() == 77 and len(PD_) == 117
A = PR[PR.all_valid.astype(str).str.lower() == "true"].copy(); assert len(A) == 285
M12 = MEANS[MEANS.set == "twelve"].set_index("group"); assert len(M12) == 12 and set(M12.index) == set(LN)
M5 = MEANS[MEANS.set == "five"].set_index(["construction", "group"])
FIVE = ["CNS_tissue", "liver", "lymphoid", "myeloid", "stromal"]
ORDER = list(M12.model_mean.sort_values().index)                                       # 줄 순서: 가중 모형 평균이 작은 것(전사 시작 쪽)부터
ROWS5 = [g for g in ORDER if g in FIVE]                                                # c, d의 줄 순서도 같은 순서
for g in ORDER:                                                                        # w10의 쌍 수와 중앙값이 쌍 표와 같은지
    v = A[A.group == g].all_phiD; assert len(v) == M12.n_pairs[g] and abs(v.median() - M12["median"][g]) < 1e-9, g
for lab, col, df in (("all", "all", PR), ("res", "res", PR), ("mat", "mat", PR), ("dist", "dist", PD_), ("prox", "prox", PD_)):
    for g in FIVE:
        v = df[(df.group == g) & (df[col + "_valid"].astype(str).str.lower() == "true")][col + "_phiD"]; assert len(v) == M5.n_pairs[(lab, g)] and abs(v.median() - M5["median"][(lab, g)]) < 1e-9, (lab, g)
SE_CUT = (10.0, 20.0); DIAM = (3.9, 3.1, 2.3)                                          # 점 크기 단계: s.e. < 10°, 10–20°, > 20° (지름 pt)
def sclass(se): return 0 if se < SE_CUT[0] else (1 if se <= SE_CUT[1] else 2)
A["cls"] = A.all_phiD_se.map(sclass)

# ------------------------------------------------------------------ 그림 틀 (인치)
W = 7.2
LABW = 1.26                                                                             # 주석 이름 칸의 오른쪽 끝
AX0, AX1 = LABW + .10, 4.62                                                             # a
BX0, BX1 = 4.90, 6.80; XN = 7.08                                                        # b, 쌍 수의 오른쪽 끝
CX0, CX1 = LABW + .10, 4.02                                                             # c
DX0, DX1 = 4.30, 6.96                                                                   # d
XA = (-50.0, 125.0); XB = (-25.0, 97.0)                                                 # 가로 범위: a(쌍), b·c·d(평균; 셋이 같다)
PITCH_C = .30
SCALE = (AX1 - AX0) * 72 / (XA[1] - XA[0])                                              # a: pt per degree
def swarm_row(items, hb):
    """한 줄 안에서 점이 겹치지 않는 가장 작은 세로 어긋남(0, ±0.5, ±1, … pt). items = [(x pt, 지름 pt)] — 놓는 순서대로(s.e.가 작은 점부터 가운데 줄에). 자리가 없으면 None."""
    cand = [0.0] + [s * k * .5 for k in range(1, 80) for s in (1, -1)]; placed = []; out = []
    for x, d in items:
        for off in cand:
            if abs(off) + d / 2 > hb: continue
            if all((x - px) ** 2 + (off - py) ** 2 >= ((d + pd_) / 2 + .45) ** 2 for px, py, pd_ in placed): break
        else: return None
        placed.append((x, off, d)); out.append(off)
    return out
A = A.sort_values(["group", "cls", "all_phiD", "id"]).copy(); A["off"] = 0.0; HB = {}
for g in ORDER:                                                                         # 줄마다 점이 겹치지 않는 가장 작은 반폭(pt) → 줄 높이(줄마다 다르다; a와 b가 같은 줄을 쓴다)
    sub = A[A.group == g]; items = [((r.all_phiD - XA[0]) * SCALE, DIAM[r.cls]) for r in sub.itertuples()]; hb = 2.5
    while (offs := swarm_row(items, hb)) is None: hb += .25
    A.loc[sub.index, "off"] = offs; HB[g] = hb
ROWH = {g: max(2 * HB[g] + 2.4, 11.5) for g in ORDER}                                   # 줄 높이(pt)
EDGE_A = np.concatenate([[0.0], np.cumsum([ROWH[g] for g in ORDER])]); YC_A = {g: (EDGE_A[i] + EDGE_A[i + 1]) / 2 for i, g in enumerate(ORDER)}   # 줄 경계와 가운데(pt; 위 = 0)
EDGE_C = np.arange(6) * PITCH_C * 72.0
print("a: row heights (pt):", {g: round(ROWH[g], 1) for g in ORDER})
HEAD = .20                                                                              # 줄 위의 머리글(CAGE, equal, DNase) 자리
KEY = .26                                                                               # 축 이름 아래 기호 풀이 자리
AXLAB = .34                                                                             # 눈금과 축 이름 자리
MT, MB, GAPV = .10, .03, .14
HA = EDGE_A[-1] / 72; HC = 5 * PITCH_C
H = MT + HEAD + HA + AXLAB + KEY + GAPV + HEAD + HC + AXLAB + KEY + MB
Y_A0 = H - MT - HEAD - HA                                                               # a·b 축의 아래 끝
Y_C0 = MB + KEY + AXLAB                                                                 # c·d 축의 아래 끝
fig = plt.figure(figsize=(W, H)); RD = fig.canvas.get_renderer()
ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(0, W); ov.set_ylim(0, H); ov.axis("off"); ov.set_zorder(5)
def T(x, y, s, **kw):
    kw.setdefault("fontsize", FS_NOTE); kw.setdefault("color", INK); kw.setdefault("ha", "left"); kw.setdefault("va", "center")
    if "$" in s: kw.setdefault("math_fontfamily", "stixsans")
    return ov.text(x, y, s, **kw)
def mk(x, y, **kw): ov.plot([x], [y], ls="none", clip_on=False, **kw)
def letter(x, y, s): ov.text(x, y, s, fontsize=FS_LET, fontweight="bold", ha="left", va="top", color=INK)
def panel(x0, x1, y0, h, xlim, edges, xlabel):
    """줄 그림 하나(세로 좌표 = 위에서부터의 pt): 하나 걸러 옅은 띠, 0°(CAGE 색)·45°(점선)·90°(DNase 색) 세로선, 위쪽 머리글, 아래 눈금."""
    ax = fig.add_axes([x0 / W, y0 / H, (x1 - x0) / W, h / H])
    for i in range(0, len(edges) - 1, 2): ax.axhspan(edges[i], edges[i + 1], color=BAND, lw=0, zorder=0)
    ax.axvline(0, color=COL["CAGE"], lw=1.1, zorder=1); ax.axvline(90, color=COL["DNASE"], lw=1.1, zorder=1); ax.axvline(45, color=GRID, lw=.7, ls=(0, (3, 2)), zorder=1)
    ax.set_xlim(*xlim); ax.set_ylim(edges[-1], 0); ax.set_yticks([]); ax.set_xticks([0, 30, 60, 90]); ax.set_xticklabels(["0°", "30°", "60°", "90°"])
    for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(MUTED); ax.spines["bottom"].set_linewidth(.6)
    ax.tick_params(labelsize=FS_NOTE, length=2.5, width=.6, color=MUTED, labelcolor=INK, pad=2)
    ax.set_xlabel(xlabel, fontsize=FS_NOTE, color=INK, labelpad=2)
    xin = lambda d: x0 + (d - xlim[0]) / (xlim[1] - xlim[0]) * (x1 - x0); yh = y0 + h + .10
    T(xin(0), yh, "CAGE", ha="center", fontweight="bold"); T(xin(90), yh, "DNase", ha="center", fontweight="bold"); T(xin(45), yh, "equal", ha="center", color=MUTED)
    return ax, xin
def rownames(names, y0, h, centres):
    for s, yc in zip(names, centres): T(LABW, y0 + h - yc / 72, s, ha="right")
    T(LABW, y0 + h + .10, "Tissue-group annotation", ha="right", color=MUTED)
FILL = lambda d, c=NEUTRAL: dict(marker="o", ms=d, mfc=c, mec="white", mew=.5)
NUM = {"definition": "CAGE–DNase balance = atan2(DNase tau*, CAGE tau*) in degrees from the joint model of the CAGE, DNase and H3K27ac annotations of a tissue group; 0 = CAGE only (transcription initiation), "
                     "90 = DNase only (accessibility), 45 = equal; defined when CAGE tau* + DNase tau* > 0",
       "unit": "enriched GWAS–tissue-group pair: single-group annotation tau* > 0 with one-sided P < 0.05/12; 288 pairs from 77 of 106 GWAS; balance defined for 285 (all-track annotations)",
       "colours": dict(points=NEUTRAL, line_0deg_CAGE=COL["CAGE"], line_90deg_DNase=COL["DNASE"])}

# ------------------------------------------------------------------ a 쌍마다의 균형 (줄 = 조직 묶음 주석)
axa, xin_a = panel(AX0, AX1, Y_A0, HA, XA, EDGE_A, BAL + " of each enriched GWAS–tissue-group pair")
print("a: largest vertical offset per row (pt):", A.groupby("group").off.agg(lambda v: v.abs().max()).to_dict())
for r in A.sort_values("cls", ascending=False).itertuples():                           # s.e.가 작은(큰) 점을 위에 그린다
    assert XA[0] + 2 < r.all_phiD < XA[1] - 2, (r.id, r.all_phiD)
    axa.plot([r.all_phiD], [YC_A[r.group] - r.off], ls="none", zorder=3 + .1 * (2 - r.cls), **FILL(DIAM[r.cls]))   # 세로 좌표는 pt(위 = 0); 어긋남 +는 위쪽
rownames([LN[g] for g in ORDER], Y_A0, HA, [YC_A[g] for g in ORDER])
NUM["a_row_order_top_to_bottom"] = ORDER
NUM["a_counts"] = dict(pairs=len(A), gwas=int(A.id.nunique()), by_annotation={g: int((A.group == g).sum()) for g in ORDER},
                       by_size_class={("<10", "10-20", ">20")[k]: int(v) for k, v in A.cls.value_counts().sort_index().items()},
                       below_45=int((A.all_phiD < 45).sum()), above_45=int((A.all_phiD > 45).sum()), below_0=int((A.all_phiD < 0).sum()), above_90=int((A.all_phiD > 90).sum()),
                       min=float(A.all_phiD.min()), max=float(A.all_phiD.max()),
                       cage_gt_dnase_p05=int((PR.all_c1 / PR.all_c1_se > 1.645).sum()), dnase_gt_cage_p05=int((PR.all_c1 / PR.all_c1_se < -1.645).sum()))
NUM["a_pairs"] = [dict(id=r.id, trait=r.trait, category=r.category, annotation=r.group, balance=float(r.all_phiD), se=float(r.all_phiD_se), size_class=("<10", "10-20", ">20")[r.cls]) for r in A.itertuples()]

# ------------------------------------------------------------------ b 주석마다: 가중 모형의 평균 ± 95% CI, 중앙값, 쌍 수
axb, xin_b = panel(BX0, BX1, Y_A0, HA, XB, EDGE_A, "Mean " + BAL + " (95% CI)")
NUM["b"] = {}
for i, g in enumerate(ORDER):
    m, se, med, n = float(M12.model_mean[g]), float(M12.se_jackknife[g]), float(M12["median"][g]), int(M12.n_pairs[g]); lo, hi = m - Z * se, m + Z * se; assert XB[0] < lo and hi < XB[1] and XB[0] < med < XB[1], g
    yc = YC_A[g]; axb.plot([lo, hi], [yc, yc], color=NEUTRAL, lw=.9, solid_capstyle="butt", zorder=2); axb.plot([m], [yc], ls="none", zorder=3, **FILL(4.4))
    axb.plot([med, med], [yc + 1.7, yc + 5.3], color=INK, lw=1.0, solid_capstyle="butt", zorder=4)   # 중앙값: 평균의 선 아래쪽에 짧은 세로선(선을 가로지르면 신뢰구간 끝처럼 보인다)
    T(XN, Y_A0 + HA - yc / 72, f"{n}", ha="right", color=MUTED)
    NUM["b"][g] = dict(name=LN[g], n_pairs=n, model_mean=m, se_jackknife=se, ci95=[lo, hi], median=med)
T(XN, Y_A0 + HA + .10, "$n$", ha="right", color=MUTED)

# ------------------------------------------------------------------ c 다섯 주석 × 세 구성
axc, xin_c = panel(CX0, CX1, Y_C0, HC, XB, EDGE_C, "Mean " + BAL + " (95% CI)")
CSTY = [("all", "all tracks", dict(marker="o", ms=4.4, mfc=NEUTRAL, mec="white", mew=.5), -.27),
        ("res", "three units per layer", dict(marker="s", ms=3.7 - .85, mfc="white", mec=NEUTRAL, mew=.85), 0.0),
        ("mat", "cell-type-matched", dict(marker="D", ms=3.9, mfc=INK, mec="white", mew=.5), .27)]
def rows5(ax, stys, key):
    out = {}
    for i, g in enumerate(ROWS5):
        for lab, name, sty, dy in stys:
            m, se = float(M5.model_mean[(lab, g)]), float(M5.se_jackknife[(lab, g)]); lo, hi = m - Z * se, m + Z * se; assert XB[0] < lo and hi < XB[1], (lab, g, lo, hi)
            yy = (i + .5 + dy) * PITCH_C * 72; ax.plot([lo, hi], [yy, yy], color=sty.get("mec") if sty["mfc"] == "white" else sty["mfc"], lw=.9, solid_capstyle="butt", zorder=2); ax.plot([m], [yy], ls="none", zorder=3, **sty)
            out.setdefault(g, {})[lab] = dict(n_pairs=int(M5.n_pairs[(lab, g)]), model_mean=m, se_jackknife=se, ci95=[lo, hi], median=float(M5["median"][(lab, g)]))
    return out
NUM["c"] = rows5(axc, CSTY, "c"); NUM["c_row_order_top_to_bottom"] = ROWS5
rownames([LN[g] for g in ROWS5], Y_C0, HC, [(i + .5) * PITCH_C * 72 for i in range(5)])

# ------------------------------------------------------------------ d 다섯 주석: TSS에서 1 kb 밖과 안
axd, xin_d = panel(DX0, DX1, Y_C0, HC, XB, EDGE_C, "Mean " + BAL + " (95% CI)")
DSTY = [("dist", "more than 1 kb from a TSS", dict(marker="o", ms=4.4, mfc=NEUTRAL, mec="white", mew=.5), -.17),
        ("prox", "within 1 kb of a TSS", dict(marker="o", ms=4.4 - .85, mfc="white", mec=NEUTRAL, mew=.85), .17)]
NUM["d"] = rows5(axd, DSTY, "d")
def con(c, k):
    r = CON[(CON.construction == c) & (CON.contrast == k)]; assert len(r) == 1, (c, k); r = r.iloc[0]; return dict(estimate=float(r.estimate), se_jackknife=float(r.se_jackknife), p=float(r.p))
NUM["stromal_minus_mean_of_four"] = {c: con(c, "stromal − mean of four") for c in ("all", "res", "mat", "dist", "prox", "dist − prox")}

# ------------------------------------------------------------------ 기호 풀이 (패널 아래)
def key_row(x, y, items):
    """표지와 글을 가로로 늘어놓는다. items = [(그리는 함수(x, y), 글)]."""
    for draw, s in items:
        draw(x + .05, y); t = T(x + .15, y, s); x += .15 + t.get_window_extent(RD).width / fig.dpi + .16
    return x
YK_A = Y_A0 - AXLAB - .16
x = AX0; t = T(x, YK_A, "Standard error of the balance", color=MUTED); x += t.get_window_extent(RD).width / fig.dpi + .12
key_row(x, YK_A, [((lambda xx, yy, d=d: mk(xx, yy, **FILL(d))), s) for d, s in zip(DIAM, [f"< {SE_CUT[0]:.0f}°", f"{SE_CUT[0]:.0f}–{SE_CUT[1]:.0f}°", f"> {SE_CUT[1]:.0f}°"])])
def k_mean(xx, yy): ov.plot([xx - .07, xx + .07], [yy, yy], color=NEUTRAL, lw=.9, solid_capstyle="butt"); mk(xx, yy, **FILL(4.4))
def k_med(xx, yy): ov.plot([xx, xx], [yy - .025, yy + .025], color=INK, lw=1.0, solid_capstyle="butt")
key_row(BX0 + .04, YK_A, [(k_mean, "mean"), (k_med, "median of the pairs")])
YK_C = Y_C0 - AXLAB - .16
key_row(CX0, YK_C, [((lambda xx, yy, sty=sty: mk(xx, yy, **sty)), name) for _, name, sty, _ in CSTY])
key_row(DX0, YK_C, [((lambda xx, yy, sty=sty: mk(xx, yy, **sty)), name) for _, name, sty, _ in DSTY])

# ------------------------------------------------------------------ 패널 글자
YL_A = H - .03; YL_C = Y_C0 + HC + HEAD + .06
letter(.04, YL_A, "a"); letter(BX0 - .22, YL_A, "b"); letter(.04, YL_C, "c"); letter(DX0 - .22, YL_C, "d")

# ------------------------------------------------------------------ 겹침·잘림 점검 후 저장
def text_issues(fig):
    fig.canvas.draw(); rd = fig.canvas.get_renderer(); fb = fig.bbox
    ts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    for t in ts:
        s = t.get_text(); assert "τ" not in s and "φ" not in s and "\t" not in s, repr(s)
    bb = [(t, matplotlib.text.Text.get_window_extent(t, rd)) for t in ts]
    ovl = [(a.get_text(), b.get_text()) for i, (a, x) in enumerate(bb) for b, y in bb[i + 1:] if x.overlaps(y)]
    cl = [t.get_text() for t, x in bb if x.x0 < fb.x0 - 1 or x.x1 > fb.x1 + 1 or x.y0 < fb.y0 - 1 or x.y1 > fb.y1 + 1]
    return ovl, cl
ovl, cl = text_issues(fig); print("text overlaps:", ovl); print("clipped:", cl)
fig.savefig(OUT / "fig4.png", dpi=300); fig.savefig(OUT / "fig4.pdf"); plt.close(fig)
NUM["inputs"] = USED; NUM["figure_size_in"] = [W, H]
(OUT / "numbers_fig4.json").write_text(json.dumps(NUM, indent=1, default=float, ensure_ascii=False), encoding="utf-8")
print("a:", NUM["a_counts"]); print("b:", {LN[g]: (round(v["model_mean"], 1), round(v["median"], 1), v["n_pairs"]) for g, v in NUM["b"].items()})
print("c:", {LN[g]: {k: round(v["model_mean"], 1) for k, v in d.items()} for g, d in NUM["c"].items()}); print("d:", {LN[g]: {k: round(v["model_mean"], 1) for k, v in d.items()} for g, d in NUM["d"].items()})
print("stromal − mean of four:", {k: (round(v["estimate"], 1), round(v["se_jackknife"], 1)) for k, v in NUM["stromal_minus_mean_of_four"].items()})
assert not ovl and not cl
print(f"saved {OUT / 'fig4.png'} ({W} × {H:.2f} in)")
