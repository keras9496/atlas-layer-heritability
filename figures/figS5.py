# -*- coding: utf-8 -*-
"""Supplementary Figure S5 — 예측 효과 크기 구간마다의 CAGE − DNase 대비 (벤치마크 14형질).
 형질마다 조직 묶음의 CAGE·DNase·H3K27ac 예측 효과를 겹치지 않는 네 구간(상위 0.1%, 0.1–0.5%, 0.5–1%, 1–5%)으로 나눈 주석 12개를 함께 적합한 모형에서, 구간마다 CAGE τ* − DNase τ*.
 줄 = 14형질(그림 5와 같은 순서; 옅은 띠는 조직 묶음을 하나 걸러 나눈다), 칸 = 구간, 점 = 대비, 가로 구간 = 95% CI, 세로선 = 0. 네 칸의 가로 범위가 같다.
 점의 색 = 층위를 함께 적합한 모형(층위당 세 단위)에서 τ*가 가장 큰 층위(키·BMD는 DNase, 나머지는 CAGE; 그림 5와 같다). 그 모형의 대비는 w04_layers.csv의 c1.
 자료: e22_benchmark14/results/w04_bins.csv, w04_layers.csv. 새 통계 없음.
실행: python figures/figS5.py      (저장소 뿌리에서; 모두 한 번에는 python figures/make_figures.py)
출력: figures/out/figS5.png, figS5.pdf, numbers_figS5.json
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from common import R, COL, TN, INK, MUTED

plt.rcParams.update({"font.family": "Arial", "font.size": 7, "text.color": INK, "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
                     "pdf.fonttype": 42, "axes.unicode_minus": True})
from matplotlib import _mathtext as _mt
def _xheight_from_glyph(self, fontname, fontsize, dpi):
    return self.get_metrics(fontname, matplotlib.rcParams["mathtext.default"], "x", fontsize, dpi).iceberg
_mt.TruetypeFonts.get_xheight = _xheight_from_glyph
FS, FSS, FS_LET, Z = 7.0, 6.5, 9, 1.96
OUT = Path(__file__).resolve().parent / "out"; OUT.mkdir(parents=True, exist_ok=True)
W3 = R / "e22_benchmark14" / "results"
B14 = ["EA", "N_SCZ", "N_BIP", "NEURO", "N_BMI", "N_AD", "MONO", "NEUT", "IBD", "RA", "N_LDL", "TC", "N_HEIGHT", "BMD"]; GSIZE = [5, 3, 2, 2, 2]
LAYER = {t: ("DNASE" if t in ("N_HEIGHT", "BMD") else "CAGE") for t in B14}; NAME = {t: TN[t[2:] if t.startswith("N_") else t] for t in B14}
BINS = ["0–0.1%", "0.1–0.5%", "0.5–1%", "1–5%"]; HEADS = ["top 0.1%", "0.1–0.5%", "0.5–1%", "1–5%"]
D = pd.read_csv(W3 / "w04_bins.csv").set_index(["trait", "bin"]); assert len(D) == 56 and set(D.index.get_level_values(1)) == set(BINS)
LY = pd.read_csv(W3 / "w04_layers.csv").set_index("trait").loc[B14]; assert all((LY.c1[t] < 0) == (LAYER[t] == "DNASE") for t in B14)       # 재추출 모형의 대비의 부호와 색이 같다

W = 7.2; XL0, MR, GAP = 1.30, .08, .16; PW = (W - XL0 - MR - 3 * GAP) / 4
PITCH = .215; n = len(B14); PH = n * PITCH; MB, MT = .46, .36; H = MB + PH + MT
XLIM = (-5.8, 10.3); TICKS = [-4, 0, 4, 8]
fig = plt.figure(figsize=(W, H)); ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(0, W); ov.set_ylim(0, H); ov.axis("off"); ov.set_zorder(5)
edges = np.concatenate([[0], np.cumsum(GSIZE)]); NUM = {"row_order": B14, "names": NAME, "layer": LAYER, "x_range": list(XLIM), "bins": {}, "same_sign_as_three_unit_model": {}}
for k, (b, head) in enumerate(zip(BINS, HEADS)):
    x0 = XL0 + k * (PW + GAP); ax = fig.add_axes([x0 / W, MB / H, PW / W, PH / H])
    for j in range(0, len(GSIZE), 2): ax.axhspan(edges[j] - .5, edges[j + 1] - .5, color="#f3f3f1", lw=0, zorder=0)
    ax.axvline(0, color=MUTED, lw=.6, zorder=1); same = 0; NUM["bins"][b] = {}
    for i, t in enumerate(B14):
        c1, se = float(D.c1[(t, b)]), float(D.c1_se[(t, b)]); col = COL[LAYER[t]]; assert XLIM[0] < c1 - Z * se and c1 + Z * se < XLIM[1], (t, b, c1, se)
        ax.plot([c1 - Z * se, c1 + Z * se], [i, i], color=col, lw=.9, solid_capstyle="butt", zorder=2); ax.plot([c1], [i], marker="o", ms=3.8, mfc=col, mec="white", mew=.5, ls="none", zorder=3)
        same += int((c1 > 0) == (LY.c1[t] > 0)); NUM["bins"][b][t] = dict(contrast=c1, se=se)
    NUM["same_sign_as_three_unit_model"][b] = same
    ax.set_xlim(*XLIM); ax.set_ylim(n - .5, -.5); ax.set_xticks(TICKS); ax.set_yticks([])
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(MUTED); ax.spines["bottom"].set_linewidth(.6)
    ax.tick_params(labelsize=FS, length=2.5, width=.6, color=MUTED, labelcolor=INK, pad=2)
    ax.set_xlabel("CAGE − DNase (95% CI)", fontsize=FS, color=INK, labelpad=2)
    if k == 0:
        for i, t in enumerate(B14): ax.text(-.035, i, NAME[t], transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=FS, color=INK)
    ov.text(x0 + PW / 2, MB + PH + .09, head, fontsize=FS, color=INK, ha="center", va="center")
ov.text(XL0 + (W - XL0 - MR) / 2, MB + PH + .26, "SNPs ranked by predicted effect", fontsize=FSS, color=MUTED, ha="center", va="center")
yk = MB - .30; ov.text(.06, yk + .13, r"Layer with the largest $\tau$*", fontsize=FSS, color=MUTED, ha="left", va="center", math_fontfamily="stixsans"); x = .06; rd = fig.canvas.get_renderer()
for lay, s in (("CAGE", "CAGE"), ("DNASE", "DNase")):
    ov.plot([x + .05], [yk], marker="o", ms=3.8, mfc=COL[lay], mec="white", mew=.5, ls="none", clip_on=False); t_ = ov.text(x + .13, yk, s, fontsize=FSS, color=INK, ha="left", va="center")
    fig.canvas.draw(); x += .13 + t_.get_window_extent(rd).width / fig.dpi + .14
fig.canvas.draw(); fb = fig.bbox
ts = [t_ for t_ in fig.findobj(matplotlib.text.Text) if t_.get_visible() and t_.get_text().strip()]; ext = [(t_, matplotlib.text.Text.get_window_extent(t_, rd)) for t_ in ts]
cl = [t_.get_text() for t_, e in ext if e.x0 < fb.x0 or e.x1 > fb.x1 or e.y0 < fb.y0 or e.y1 > fb.y1]; ovl = [(a.get_text(), b.get_text()) for i, (a, x_) in enumerate(ext) for b, y_ in ext[i + 1:] if x_.overlaps(y_)]
print("clipped:", cl, "| overlaps:", ovl); assert not cl and not ovl
fig.savefig(OUT / "figS5.png", dpi=300); fig.savefig(OUT / "figS5.pdf"); plt.close(fig)
NUM["three_unit_model_contrast"] = {t: dict(contrast=float(LY.c1[t]), se=float(LY.c1_se[t])) for t in B14}
(OUT / "numbers_figS5.json").write_text(json.dumps(NUM, indent=1, ensure_ascii=False), encoding="utf-8")
print(f"saved {OUT / 'figS5.png'} ({W} x {H:.2f} in); same sign as the three-unit model: {NUM['same_sign_as_three_unit_model']}")
