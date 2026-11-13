# -*- coding: utf-8 -*-
"""Supplementary Figure S4 — stromal 주석에 농축된 쌍 32개의 CAGE–DNase 균형 (106개 GWAS 가운데 stromal 주석이 농축된 GWAS).
 줄 = GWAS 하나(층위당 세 단위 주석의 균형의 s.e.가 작은 것부터; 세 단위 균형이 정의되지 않는 쌍은 맨 아래).
 진한 원과 선 = 층위당 세 단위 주석(재추출 20회 평균)의 균형 ± 95% CI. 옅은 표지 = 같은 쌍의 전체 트랙 주석(빈 원)과 세포형 일치 주석(빈 마름모)의 균형(구간은 그리지 않는다).
 95% CI가 가로 범위를 넘으면 끝에서 자르고 작은 세모로 표시한다. 균형이 정의되지 않는(CAGE τ* + DNase τ* ≤ 0) 구성은 그리지 않는다.
 오른쪽 칸 = 형질 범주(기술용). 세로선: 0°(CAGE만), 45°(같음), 90°(DNase만) — 그림 4와 같다.
 자료: e22_benchmark14/results/w05_pairs.csv. 새 통계 없음.
실행: python figures/figS4.py      (저장소 뿌리에서; 모두 한 번에는 python figures/make_figures.py)
출력: figures/out/figS4.png, figS4.pdf, numbers_figS4.json
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from common import R as RESROOT, COL, INK, MUTED, GRID

plt.rcParams.update({"font.family": "Arial", "font.size": 7, "text.color": INK, "pdf.fonttype": 42, "axes.unicode_minus": True})
FS, FSS, Z = 7.0, 6.5, 1.96
BAL = "CAGE–DNase balance"; NEUTRAL, LIGHT, BAND = "#555555", "#8a8a8a", "#f3f3f1"
OUT = Path(__file__).resolve().parent / "out"; OUT.mkdir(parents=True, exist_ok=True)
PR = pd.read_csv(RESROOT / "e22_benchmark14" / "results" / "w05_pairs.csv"); S = PR[PR.group == "stromal"].copy(); assert len(S) == 32
for c in ("all", "res", "mat"): S[c + "_ok"] = S[c + "_valid"].astype(str).str.lower() == "true"
assert S.all_ok.all() and S.res_ok.sum() == 31 and S.mat_ok.sum() == 30
NAME = {"T055": "Forced vital capacity", "T070": "Heel bone mineral density", "T086": "Balding", "T009": "Glaucoma", "T060": "Height", "T005": "Chronic obstructive pulmonary disease", "T031": "Inguinal hernia",
        "T064": "Waist–hip ratio adjusted for BMI", "T015": "Total hip replacement", "T028": "Varicose veins", "T072": "Birth weight, maternal effect", "T057": "Diastolic blood pressure",
        "T106": "Structural connectivity, global (brain imaging)", "T026": "Primary open-angle glaucoma", "T036": "Abdominal aortic aneurysm", "T083": "Platelet distribution width", "T082": "Reticulocyte count",
        "T004": "Knee osteoarthritis", "T069": "Birth weight, fetal effect", "T008": "Atrial fibrillation", "T076": "Aspartate aminotransferase", "T074": "IGF-1", "T016": "Diverticulosis and diverticulitis",
        "T014": "Gastrointestinal disease", "T063": "Total protein", "T032": "Endometriosis", "T011": "Coronary artery disease", "T092": "Total cholesterol", "T089": "Phosphate", "T043": "Acute appendicitis",
        "T087": "Alkaline phosphatase", "T066": "Creatinine"}
assert set(S.id) == set(NAME), set(S.id) ^ set(NAME)
S = pd.concat([S[S.res_ok].sort_values(["res_phiD_se", "id"]), S[~S.res_ok].sort_values("id")]); n = len(S)

W = 7.2; XNAME = 2.50; AX0, AX1 = XNAME + .10, 6.02; XCAT = AX1 + .12
PITCH = .158; PH = n * PITCH; MB, MT = .62, .30; H = MB + PH + MT; XL = (-60.0, 150.0)
fig = plt.figure(figsize=(W, H)); ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(0, W); ov.set_ylim(0, H); ov.axis("off"); ov.set_zorder(5)
ax = fig.add_axes([AX0 / W, MB / H, (AX1 - AX0) / W, PH / H])
for i in range(0, n, 2): ax.axhspan(i - .5, i + .5, color=BAND, lw=0, zorder=0)
ax.axvline(0, color=COL["CAGE"], lw=1.1, zorder=1); ax.axvline(90, color=COL["DNASE"], lw=1.1, zorder=1); ax.axvline(45, color=GRID, lw=.7, ls=(0, (3, 2)), zorder=1)
M_RES = dict(marker="o", ms=4.2, mfc=NEUTRAL, mec="white", mew=.5, ls="none"); M_ALL = dict(marker="o", ms=3.4 - .7, mfc="white", mec=LIGHT, mew=.7, ls="none"); M_MAT = dict(marker="D", ms=3.1 - .7, mfc="white", mec=LIGHT, mew=.7, ls="none")
NUM = {"rows_top_to_bottom": [], "x_range": list(XL), "truncated_intervals": []}
for i, r in enumerate(S.itertuples()):
    if r.all_ok and XL[0] < r.all_phiD < XL[1]: ax.plot([r.all_phiD], [i - .22], zorder=2.5, **M_ALL)
    if r.mat_ok and XL[0] < r.mat_phiD < XL[1]: ax.plot([r.mat_phiD], [i + .22], zorder=2.5, **M_MAT)
    if r.res_ok:
        lo, hi = r.res_phiD - Z * r.res_phiD_se, r.res_phiD + Z * r.res_phiD_se; l2, h2 = max(lo, XL[0] + 2.5), min(hi, XL[1] - 2.5)
        ax.plot([l2, h2], [i, i], color=NEUTRAL, lw=.9, solid_capstyle="butt", zorder=2); ax.plot([r.res_phiD], [i], zorder=3, **M_RES)
        if lo < l2: ax.plot([l2], [i], marker="<", ms=2.6, mfc=NEUTRAL, mec="none", ls="none", zorder=2)
        if hi > h2: ax.plot([h2], [i], marker=">", ms=2.6, mfc=NEUTRAL, mec="none", ls="none", zorder=2)
        if lo < l2 or hi > h2: NUM["truncated_intervals"].append(NAME[r.id])
    y = MB + PH - (i + .5) * PITCH
    ov.text(XNAME, y, NAME[r.id], ha="right", va="center", fontsize=FS, color=INK); ov.text(XCAT, y, r.category, ha="left", va="center", fontsize=FSS, color=MUTED)
    NUM["rows_top_to_bottom"].append(dict(id=r.id, name=NAME[r.id], panel_label=r.trait, category=r.category,
                                          three_units=dict(balance=float(r.res_phiD), se=float(r.res_phiD_se), defined=bool(r.res_ok)), all_tracks=dict(balance=float(r.all_phiD), se=float(r.all_phiD_se), defined=bool(r.all_ok)),
                                          cell_type_matched=dict(balance=float(r.mat_phiD), se=float(r.mat_phiD_se), defined=bool(r.mat_ok))))
ax.set_xlim(*XL); ax.set_ylim(n - .5, -.5); ax.set_yticks([]); ax.set_xticks([0, 30, 60, 90, 120]); ax.set_xticklabels(["0°", "30°", "60°", "90°", "120°"])
for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
ax.spines["bottom"].set_color(MUTED); ax.spines["bottom"].set_linewidth(.6)
ax.tick_params(labelsize=FSS, length=2.5, width=.6, color=MUTED, labelcolor=INK, pad=2)
ax.set_xlabel(BAL + " in the stromal/mesenchymal annotation", fontsize=FSS, color=INK, labelpad=2)
xin = lambda d: AX0 + (d - XL[0]) / (XL[1] - XL[0]) * (AX1 - AX0); yh = MB + PH + .10
ov.text(xin(0), yh, "CAGE", ha="center", va="center", fontsize=FSS, fontweight="bold", color=INK); ov.text(xin(90), yh, "DNase", ha="center", va="center", fontsize=FSS, fontweight="bold", color=INK)
ov.text(xin(45), yh, "equal", ha="center", va="center", fontsize=FSS, color=MUTED)
ov.text(XNAME, yh, "GWAS enriched in the annotation", ha="right", va="center", fontsize=FSS, color=MUTED); ov.text(XCAT, yh, "Trait category", ha="left", va="center", fontsize=FSS, color=MUTED)
yk = MB - .50; x = AX0 - .9; rd = fig.canvas.get_renderer()
for kind, s in (("res", "three units per layer (95% CI)"), ("all", "all tracks"), ("mat", "cell-type-matched")):
    if kind == "res": ov.plot([x - .03, x + .13], [yk, yk], color=NEUTRAL, lw=.9, solid_capstyle="butt"); ov.plot([x + .05], [yk], clip_on=False, **M_RES)
    else: ov.plot([x + .05], [yk], clip_on=False, **(M_ALL if kind == "all" else M_MAT))
    t_ = ov.text(x + .19, yk, s, fontsize=FSS, color=INK, ha="left", va="center"); fig.canvas.draw(); x += .19 + t_.get_window_extent(rd).width / fig.dpi + .20

fig.canvas.draw(); fb = fig.bbox
ts = [t_ for t_ in fig.findobj(matplotlib.text.Text) if t_.get_visible() and t_.get_text().strip()]; ext = [(t_, matplotlib.text.Text.get_window_extent(t_, rd)) for t_ in ts]
cl = [t_.get_text() for t_, e in ext if e.x0 < fb.x0 or e.x1 > fb.x1 or e.y0 < fb.y0 or e.y1 > fb.y1]; ovl = [(a.get_text(), b.get_text()) for i, (a, x_) in enumerate(ext) for b, y_ in ext[i + 1:] if x_.overlaps(y_)]
print("clipped:", cl, "| overlaps:", ovl); assert not cl and not ovl
fig.savefig(OUT / "figS4.png", dpi=300); fig.savefig(OUT / "figS4.pdf"); plt.close(fig)
NUM.update(n_pairs=n, n_three_units_defined=int(S.res_ok.sum()), n_cell_type_matched_defined=int(S.mat_ok.sum()), by_category=S.category.value_counts().to_dict(),
           ten_most_precise=[(NAME[r.id], round(float(r.res_phiD), 1), round(float(r.res_phiD_se), 1)) for r in S.head(10).itertuples()])
(OUT / "numbers_figS4.json").write_text(json.dumps(NUM, indent=1, ensure_ascii=False), encoding="utf-8")
print(f"saved {OUT / 'figS4.png'} ({W} x {H:.2f} in); truncated: {NUM['truncated_intervals']}"); print(NUM["ten_most_precise"])
