# -*- coding: utf-8 -*-
"""Supplementary Figure S3 — 조직 묶음 × 조절 층위 지도: 그림 2a에 없는 벤치마크 아홉 형질.
 그림 2a와 같은 띠: 띠 안에서 열 = 조직 묶음 12개, 줄 = 조절 층위 6개. 칸 색 = τ*(주석 하나를 baselineLD v2.2와 함께 넣은 모형),
 점 = FDR < 0.05(14형질 × 72주석 = 1,008칸에 걸친 BH). 띠 왼쪽의 세로선 = 유형의 색(그림 2와 같은 색), 색 풀이는 아래쪽.
 형질 순서는 그림 1b와 같다(그림 2a의 다섯을 뺀 아홉). 자료: e22_benchmark14/results/w04_map72.csv, w04_types.csv. 새 통계 없음.
실행: python figures/figS3.py      (저장소 뿌리에서; 모두 한 번에는 python figures/make_figures.py)
출력: figures/out/figS3.png, figS3.pdf, numbers_figS3.json
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, to_hex
import common as c
from common import INK, MUTED, LIN12, LN
from taumap import TAUMAP, wcag

plt.rcParams.update({"font.family": "Arial", "font.size": 7, "text.color": INK, "axes.edgecolor": INK, "mathtext.fontset": "custom", "mathtext.rm": "Arial",
                     "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "pdf.fonttype": 42, "axes.unicode_minus": True})
OUT = Path(__file__).resolve().parent / "out"; OUT.mkdir(parents=True, exist_ok=True)
FS, FSS = 7.0, 6.5
W4 = c.R / "e22_benchmark14" / "results"
ORD = ["EA", "N_SCZ", "N_BIP", "NEURO", "N_BMI", "N_AD", "MONO", "NEUT", "IBD", "RA", "N_LDL", "TC", "N_HEIGHT", "BMD"]
MAIN5 = ["N_SCZ", "N_AD", "MONO", "N_LDL", "N_HEIGHT"]; EX = [t for t in ORD if t not in MAIN5]; assert len(EX) == 9
NAME = {"EA": "Educational\nattainment", "N_BIP": "Bipolar\ndisorder", "NEURO": "Neuroticism", "N_BMI": "BMI", "NEUT": "Neutrophil\ncount", "IBD": "IBD", "RA": "Rheumatoid\narthritis",
        "TC": "Total\ncholesterol", "BMD": "Bone mineral\ndensity"}
MODS6 = ["CAGE", "DNASE", "H3K27ac", "H3K4me1", "H3K4me3", "RNA"]
LAB6 = {"CAGE": "CAGE", "DNASE": "DNase", "H3K27ac": "H3K27ac", "H3K4me1": "H3K4me1", "H3K4me3": "H3K4me3", "RNA": "RNA-seq"}
TYPE_COL = {"focal": "#555555", "layer-spanning": "#C4A000", "tissue-spanning": "#D1478F", "tissue- and layer-spanning": "#6A2C91"}      # fig2.py와 같다
TYPE_EN = {"focal": "Focal", "layer-spanning": "Layer-spanning", "tissue-spanning": "Tissue-spanning", "tissue- and layer-spanning": "Tissue- and layer-spanning"}
assert f"TYPE_COL = {json.dumps(TYPE_COL)}" .replace("'", '"') in (Path(__file__).resolve().parent / "fig2.py").read_text(encoding="utf-8").replace("'", '"')   # 그림 2와 같은 색
D1 = pd.read_csv(W4 / "w04_map72.csv"); assert len(D1) == 1008 and int((D1.q < .05).sum()) == 107
V3 = pd.read_csv(W4 / "w04_types.csv", index_col=0)
def grid(t, col): return D1[D1.trait == t].pivot(index="layer", columns="group", values=col).reindex(index=MODS6, columns=LIN12)
TAUS = {t: grid(t, "tau") for t in EX}; QS = {t: grid(t, "q") for t in EX}
VMIN, VMAX = -2.5, 4.5; NORM = TwoSlopeNorm(vmin=VMIN, vcenter=0, vmax=VMAX)
CLIPPED = [(t, l, m, round(float(TAUS[t].loc[m, l]), 2)) for t in EX for l in LIN12 for m in MODS6 if not VMIN <= TAUS[t].loc[m, l] <= VMAX]   # 색 눈금 밖의 값(색만 끝 색으로)

W = 7.2; HX0, HW = 2.08, 4.40; RH, GAPB = .105, .11; BAND = 6 * RH; KEYH = .34; YB = KEYH + .98; MT = .22
H = YB + len(EX) * BAND + (len(EX) - 1) * GAPB + MT
fig = plt.figure(figsize=(W, H))
def axin(l, b, w, h): return fig.add_axes([l / W, b / H, w / W, h / H])
ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(0, W); ov.set_ylim(0, H); ov.axis("off"); ov.set_zorder(5)
DOTS = {}; n_dots = 0
for k, t in enumerate(EX):
    y0 = YB + (len(EX) - 1 - k) * (BAND + GAPB)
    ax = axin(HX0, y0, HW, BAND)
    ax.imshow(TAUS[t].values, cmap=TAUMAP, norm=NORM, aspect="auto", extent=(-.5, 11.5, 5.5, -.5), interpolation="nearest")
    for x in np.arange(.5, 11.5): ax.axvline(x, color="white", lw=.7)
    for y in np.arange(.5, 5.5): ax.axhline(y, color="white", lw=.45)
    ii, jj = np.where(QS[t].values < .05); n_dots += len(ii); DOTS[t] = [f"{LIN12[j]}|{MODS6[i]}" for i, j in zip(ii, jj)]
    assert len(ii) == V3.K[t] and len(set(jj)) == V3["T"][t], t                                       # 유의한 주석 수와 조직 묶음 수 = w04_types.csv
    assert all(wcag(INK, TAUMAP(NORM(TAUS[t].values[i, j]))) >= 3.0 for i, j in zip(ii, jj))
    ax.plot(jj, ii, ls="none", marker="o", ms=2.8, mfc=INK, mec="white", mew=.5, zorder=4)
    ax.set_xlim(-.5, 11.5); ax.set_ylim(5.5, -.5); ax.set_yticks(range(6)); ax.set_yticklabels([LAB6[m] for m in MODS6], fontsize=FSS, color=MUTED)
    ax.tick_params(length=0, pad=3); [sp.set_visible(False) for sp in ax.spines.values()]
    if k == len(EX) - 1:
        ax.set_xticks(range(12)); ax.set_xticklabels([LN[l] for l in LIN12], rotation=40, ha="right", rotation_mode="anchor", fontsize=FS, color=INK)
    else: ax.set_xticks([])
    ov.text(HX0 - .62, y0 + BAND / 2, NAME[t], fontsize=FS, color=INK, ha="right", va="center", linespacing=1.1)
    ov.plot([.40, .40], [y0, y0 + BAND], color=TYPE_COL[V3.type[t]], lw=2.2, solid_capstyle="butt")
CBX, CBW, CBH = 6.70, .085, 1.45; ytop = H - MT
cb = axin(CBX, ytop - CBH, CBW, CBH); v = np.linspace(VMIN, VMAX, 512)
cb.imshow(TAUMAP(NORM(v))[::-1, None, :3], aspect="auto", extent=(0, 1, VMIN, VMAX))
cb.set_xticks([]); cb.yaxis.tick_right(); cb.set_yticks(range(-2, 5)); cb.set_yticklabels([("0" if q == 0 else f"{q:+d}".replace("-", "−")) for q in range(-2, 5)], fontsize=FS, color=INK)
cb.tick_params(axis="y", length=2, width=.5, color=MUTED, pad=2); [sp.set_visible(False) for sp in cb.spines.values()]; cb.set_ylim(VMIN, VMAX)
ov.text(CBX + CBW / 2, ytop + .05, r"$\tau$*", fontsize=FS, color=INK, ha="center", va="bottom", math_fontfamily="stixsans")
ky = ytop - CBH - .22
ov.plot([CBX - .03], [ky], ls="none", marker="o", ms=2.8, mfc=INK, mec="white", mew=.5); ov.text(CBX + .05, ky, "FDR\n< 0.05", fontsize=FSS, color=MUTED, ha="left", va="center", linespacing=1.15)
# 유형의 색 풀이 (아래쪽 한 줄): 이 그림에 있는 유형만
fig.canvas.draw(); rd0 = fig.canvas.get_renderer()
def tw(s): t_ = ov.text(0, 0, s, fontsize=FSS); w_ = t_.get_window_extent(rd0).width / fig.dpi; t_.remove(); return w_
present = [ty for ty in TYPE_COL if any(V3.type[t] == ty for t in EX)]
items = [TYPE_EN[ty] for ty in present]; GAPK = .24; wtot = sum(.20 + tw(s) for s in items) + GAPK * (len(items) - 1); kx = W / 2 - wtot / 2
for ty, s in zip(present, items):
    ov.plot([kx, kx + .13], [.15, .15], color=TYPE_COL[ty], lw=2.2, solid_capstyle="butt"); ov.text(kx + .20, .15, s, fontsize=FSS, color=INK, ha="left", va="center"); kx += .20 + tw(s) + GAPK

fig.canvas.draw(); rd = fig.canvas.get_renderer(); fb = fig.bbox
ts = [t_ for t_ in fig.findobj(matplotlib.text.Text) if t_.get_visible() and t_.get_text().strip()]
ext = [(t_, t_.get_window_extent(rd)) for t_ in ts]
assert not [t_.get_text() for t_, e in ext if e.x0 < fb.x0 or e.x1 > fb.x1 or e.y0 < fb.y0 or e.y1 > fb.y1]
assert all(6.5 <= t_.get_fontsize() <= 7 and to_hex(t_.get_color()) in (INK, MUTED) for t_ in ts)
flat = [(t_, e) for t_, e in ext if t_.get_rotation() % 90 == 0]; ov_ = [(a.get_text(), b.get_text()) for i, (a, x) in enumerate(flat) for b, y in flat[i + 1:] if x.overlaps(y)]; assert not ov_, ov_
fig.savefig(OUT / "figS3.png", dpi=300); fig.savefig(OUT / "figS3.pdf"); plt.close(fig)
NUM = {"source": "e22_benchmark14/results/w04_map72.csv: tau, q (BH over the 1,008 cells of the 14 benchmark traits x 72 annotations; one-sided tau* > 0); types from w04_types.csv",
       "traits": EX, "trait_names": {t: NAME[t].replace("\n", " ") for t in EX}, "types": {t: V3.type[t] for t in EX}, "K_T_L": {t: dict(K=int(V3.K[t]), T=int(V3["T"][t]), L=float(V3.L[t])) for t in EX},
       "rows": MODS6, "columns": LIN12, "colour_scale": {"vmin": VMIN, "vcenter": 0, "vmax": VMAX}, "n_dots": n_dots, "dots": DOTS, "values_outside_colour_scale": CLIPPED,
       "tau_star": {t: {l: {m: round(float(TAUS[t].loc[m, l]), 3) for m in MODS6} for l in LIN12} for t in EX}, "type_colours": TYPE_COL, "figure_size_in": [W, H]}
(OUT / "numbers_figS3.json").write_text(json.dumps(NUM, indent=1, ensure_ascii=False), encoding="utf-8")
print(f"saved {OUT / 'figS3.png'} ({W} x {H:.2f} in); dots {n_dots}; types {NUM['types']}; outside colour scale: {CLIPPED}")
