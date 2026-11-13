# -*- coding: utf-8 -*-
"""Supplementary Figure S1 — 조직 묶음에 고유한 농축 (벤치마크 14형질 × 조직 묶음 12개, 모형 M2).
 모형 M2 = baselineLD v2.2 + 모든 조직 묶음의 큰 효과 SNP 주석(L5_all) + 조직 묶음 주석 하나. 칸 색 = 조직 묶음 주석의 τ*(tau_m2).
 숫자 = q_m2 < 0.05인 칸(168칸에 걸친 BH, 단측 τ* > 0). 줄 순서·왼쪽 띠·색 범위(−2.5 ~ +4.5)는 그림 1b와 같다.
 자료: e22_benchmark14/results/w04_tissue.csv. 새 통계 없음.
실행: python figures/figS1.py      (저장소 뿌리에서; 모두 한 번에는 python figures/make_figures.py)
출력: figures/out/figS1.png, figS1.pdf, numbers_figS1.json
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, to_hex
from matplotlib.patches import Rectangle
import common as c
from common import INK, MUTED, LIN12, LN
from taumap import TAUMAP, wcag

plt.rcParams.update({"font.family": "Arial", "font.size": 7, "text.color": INK, "axes.unicode_minus": True, "mathtext.fontset": "custom",
                     "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "pdf.fonttype": 42})
OUT = Path(__file__).resolve().parent / "out"; OUT.mkdir(parents=True, exist_ok=True)
FS = 7.0
W4 = c.R / "e22_benchmark14" / "results"
ORD = ["EA", "N_SCZ", "N_BIP", "NEURO", "N_BMI", "N_AD", "MONO", "NEUT", "IBD", "RA", "N_LDL", "TC", "N_HEIGHT", "BMD"]
TN = {"EA": "Educational attainment", "N_SCZ": "Schizophrenia", "N_BIP": "Bipolar disorder", "NEURO": "Neuroticism", "N_BMI": "BMI", "N_AD": "Alzheimer's disease", "MONO": "Monocyte count",
      "NEUT": "Neutrophil count", "IBD": "IBD", "RA": "Rheumatoid arthritis", "N_LDL": "LDL cholesterol", "TC": "Total cholesterol", "N_HEIGHT": "Height", "BMD": "Bone mineral density"}
OWN = {"EA": "CNS_tissue", "N_SCZ": "CNS_tissue", "N_BIP": "CNS_tissue", "NEURO": "CNS_tissue", "N_BMI": "CNS_tissue", "N_AD": "myeloid", "MONO": "myeloid", "NEUT": "myeloid",
       "IBD": "lymphoid", "RA": "lymphoid", "N_LDL": "liver", "TC": "liver", "N_HEIGHT": "stromal", "BMD": "stromal"}
M = pd.read_csv(W4 / "w04_tissue.csv"); assert len(M) == 168 and set(M.trait) == set(ORD) and set(M.group) == set(LIN12)
T = M.pivot(index="trait", columns="group", values="tau_m2").reindex(index=ORD, columns=LIN12)
Q = M.pivot(index="trait", columns="group", values="q_m2").reindex(index=ORD, columns=LIN12)
T1 = M.pivot(index="trait", columns="group", values="tau").reindex(index=ORD, columns=LIN12); Q1 = M.pivot(index="trait", columns="group", values="q").reindex(index=ORD, columns=LIN12)
assert not T.isna().any().any() and not Q.isna().any().any()
SIG = Q.values < .05; assert SIG.sum() == 29, SIG.sum()                # w04_summary.txt: 전체 묶음 보정 모형(M2)에서 29개
assert (T.values[SIG] > 0).all()
for (t, g), v in {("EA", "CNS_tissue"): 1.15, ("N_SCZ", "CNS_tissue"): 1.82, ("N_AD", "myeloid"): 3.19, ("IBD", "lymphoid"): 3.95, ("N_LDL", "liver"): 2.65, ("N_HEIGHT", "stromal"): 2.02, ("BMD", "stromal"): 2.31}.items():
    assert abs(float(T.loc[t, g]) - v) <= .005 + 1e-9, (t, g, float(T.loc[t, g]), v)      # w04_summary.txt 1절의 M2 값
VMIN, VMAX = -2.5, 4.5; NORM = TwoSlopeNorm(vmin=VMIN, vcenter=0, vmax=VMAX)
CLIPPED = [(t, g, round(float(T.loc[t, g]), 2)) for t in ORD for g in LIN12 if not VMIN <= T.loc[t, g] <= VMAX]   # 색 눈금 밖의 값(색만 끝 색으로 그린다)
MINC = min(wcag(INK, TAUMAP(NORM(v))) for v in T.values[SIG]); assert MINC >= 4.0, MINC      # 숫자를 적는 칸 위의 INK 대비(가장 어두운 칸: IBD lymphoid)
OWN_SIG = [t for t in ORD if Q.loc[t, OWN[t]] < .05]; BOTH = int((SIG & (Q1.values < .05)).sum()); ONLY_M2 = [(ORD[i], LIN12[j]) for i, j in zip(*np.where(SIG & ~(Q1.values < .05)))]

W = 7.2; HX0, HW, HY0 = 2.08, 4.40, 0.98; nr, ncl = T.shape; HH = 2.60 * nr / 11; H = HY0 + HH + .30
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([HX0 / W, HY0 / H, HW / W, HH / H])
ax.imshow(T.values, cmap=TAUMAP, norm=NORM, aspect="auto", extent=(-.5, ncl - .5, nr - .5, -.5), interpolation="nearest")
for x in np.arange(.5, ncl - .5): ax.axvline(x, color="white", lw=.7)
for y in np.arange(.5, nr - .5): ax.axhline(y, color="white", lw=.7)
BREAKS = [4.5, 7.5, 9.5, 11.5]
for y in BREAKS: ax.axhline(y, color="white", lw=2.6)
SHOWN = []
for i in range(nr):
    for j in range(ncl):
        if SIG[i, j]:
            lab = f"{T.iat[i, j]:+.1f}"; SHOWN.append(f"{ORD[i]}|{LIN12[j]}|{lab}"); ax.text(j, i, lab, ha="center", va="center", fontsize=FS, fontweight="bold", color=INK)
ax.set_xlim(-.5, ncl - .5); ax.set_ylim(nr - .5, -.5)
ax.set_xticks(range(ncl)); ax.set_xticklabels([LN[l] for l in LIN12], rotation=40, ha="right", rotation_mode="anchor", fontsize=FS, color=INK)
ax.set_yticks(range(nr)); ax.set_yticklabels([TN[t] for t in ORD], fontsize=FS, color=INK)
ax.tick_params(length=0, pad=3); [sp.set_visible(False) for sp in ax.spines.values()]
SX, SW = .70, .07
hs = fig.add_axes([SX / W, HY0 / H, SW / W, HH / H]); hs.set_xlim(0, 1); hs.set_ylim(nr - .5, -.5); hs.axis("off")
cg = []
for i, t in enumerate(ORD):
    if not cg or cg[-1][0] != OWN[t]: cg.append([OWN[t], i, i])
    else: cg[-1][2] = i
for g, i0, i1 in cg:
    hs.add_patch(Rectangle((0, i0 - .5 + .06), 1, i1 - i0 + 1 - .12, fc="#8a8a8a", ec="none"))
    fig.text((SX - .06) / W, (HY0 + HH * (1 - ((i0 + i1) / 2 + .5) / nr)) / H, LN[g].replace("/", "/\n"), fontsize=FS, color=INK, ha="right", va="center", linespacing=1.1)
CBX, CBW, CBH = 6.70, .085, 1.55
cb = fig.add_axes([CBX / W, (HY0 + HH - CBH) / H, CBW / W, CBH / H]); v = np.linspace(VMIN, VMAX, 512)
cb.imshow(TAUMAP(NORM(v))[::-1, None, :3], aspect="auto", extent=(0, 1, VMIN, VMAX))
cb.set_xticks([]); cb.yaxis.tick_right(); cb.set_yticks(range(-2, 5)); cb.set_ylim(VMIN, VMAX)
cb.set_yticklabels([f"{q:+d}".replace("-", "−") if q else "0" for q in range(-2, 5)], fontsize=FS, color=INK)
cb.tick_params(axis="y", length=2, width=.5, color=MUTED, pad=2); [sp.set_visible(False) for sp in cb.spines.values()]
fig.text((CBX + CBW / 2) / W, (HY0 + HH + .05) / H, r"$\tau$*", fontsize=FS, color=INK, ha="center", va="bottom", math_fontfamily="stixsans")

fig.canvas.draw(); rd = fig.canvas.get_renderer(); fb = fig.bbox
ts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
ext = [(t, t.get_window_extent(rd)) for t in ts]
assert not [t.get_text() for t, e in ext if e.x0 < fb.x0 or e.x1 > fb.x1 or e.y0 < fb.y0 or e.y1 > fb.y1]
assert all(to_hex(t.get_color()) in (INK, MUTED) and 6.5 <= t.get_fontsize() <= 7 for t in ts)
flat = [(t, e) for t, e in ext if t.get_rotation() % 90 == 0]; assert not [(a.get_text(), b.get_text()) for i, (a, x) in enumerate(flat) for b, y in flat[i + 1:] if x.overlaps(y)]
fig.savefig(OUT / "figS1.png", dpi=300); fig.savefig(OUT / "figS1.pdf"); plt.close(fig)
NUM = {"source": "e22_benchmark14/results/w04_tissue.csv: tau_m2, q_m2 (model M2 = baselineLD v2.2 + annotation of SNPs with large predicted effects in any tissue group + one tissue-group annotation; BH over the 168 cells, one-sided tau* > 0)",
       "rows": [dict(key=t, name=TN[t], category_tissue_group=OWN[t]) for t in ORD], "columns": LIN12, "tau_star_m2": {t: {l: round(float(T.loc[t, l]), 3) for l in LIN12} for t in ORD},
       "q_m2": {t: {l: float(f"{Q.loc[t, l]:.4g}") for l in LIN12} for t in ORD}, "n_significant": int(SIG.sum()), "numbers_shown": SHOWN,
       "n_significant_in_single_group_model_M1": int((Q1.values < .05).sum()), "n_significant_in_both_models": BOTH, "significant_in_M2_only": [f"{t}|{g}" for t, g in ONLY_M2],
       "traits_with_category_group_significant_in_M2": [TN[t] for t in OWN_SIG], "n_traits_with_category_group_significant_in_M2": len(OWN_SIG),
       "colour_scale": {"vmin": VMIN, "vcenter": 0, "vmax": VMAX, "min_INK_contrast_on_shown_cells": round(float(MINC), 2)}, "tau_range": [round(float(T.values.min()), 3), round(float(T.values.max()), 3)], "values_outside_colour_scale": CLIPPED, "figure_size_in": [W, H]}
(OUT / "numbers_figS1.json").write_text(json.dumps(NUM, indent=1, ensure_ascii=False), encoding="utf-8")
print(f"saved {OUT / 'figS1.png'} ({W} x {H:.2f} in); significant {int(SIG.sum())}; category group significant in {len(OWN_SIG)}/14; outside colour scale: {CLIPPED}")
