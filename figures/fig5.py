# -*- coding: utf-8 -*-
"""Figure 5 — 예측 효과 주석과 기존 기준: 실측 주석, 세포형 eQTL, 예측 활성도 (벤치마크 14형질).
 줄 = 벤치마크 14형질(네 패널이 줄을 같이 쓴다). 옅은 띠는 조직 묶음을 하나 걸러 나눈다. 형질마다 τ*가 가장 큰 층위의 예측 효과 주석: 키·BMD는 DNase, 나머지는 CAGE. 점의 색 = 그 층위.
 줄 안에서 위 = 빈 점, 아래 = 채운 점, 가로 구간 = 95% CI, 세로선 = 0. 네 패널의 가로 범위가 같다.
 a 실측 주석: 빈 점 = 예측 주석 하나만(w04_map72.csv), 채운 점 = 같은 조직 묶음의 실측 CAGE·DNase·H3K27ac 주석과 함께(w04_measured.csv). 별표: 14형질 FDR < 0.05.
 b 세포형 eQTL: 빈 점 = 세 층위 효과 주석을 함께 적합한 모형(E), 채운 점 = 세포형 eQTL 주석 다섯을 더한 모형(EQ) (e18_cell_type_eqtl/results/s03_values.csv). 별표: 14형질 FDR < 0.05(단측 P = tau/se; 여기서 계산).
 c 예측 활성도: 빈 점 = 효과 주석만, 옅은 점 = 활성도 상위 5% 주석과 함께, 채운 점 = 활성도 네 구간과 함께 (w04_activity.csv, w06_actbin.csv). 별표: 네 구간 모형의 14형질 FDR < 0.05.
 d 활성 요소 안: 예측 활성도 상위 10% 안에서 활성도를 맞춘 두 절반 — 빈 점 = 효과가 작은 절반, 채운 점 = 효과가 큰 절반 (w04_activity.csv). 별표: 차이의 14형질 FDR < 0.05.
 기존 결과 파일만 읽는다. 새 적합 없음(b의 P와 FDR만 여기서 낸다 → figures/out/eqtl14.csv).
실행: python figures/fig5.py      (저장소 뿌리에서; 모두 한 번에는 python figures/make_figures.py)
출력: figures/out/fig5.png, fig5.pdf, numbers_fig5.json, eqtl14.csv
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex, to_rgb
from common import R, COL, TN, LN, INK, MUTED

plt.rcParams.update({"font.family": "Arial", "font.size": 7, "text.color": INK,
                     "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
                     "pdf.fonttype": 42, "axes.unicode_minus": True})
from matplotlib import _mathtext as _mt
def _xheight_from_glyph(self, fontname, fontsize, dpi):
    return self.get_metrics(fontname, matplotlib.rcParams["mathtext.default"], "x", fontsize, dpi).iceberg
_mt.TruetypeFonts.get_xheight = _xheight_from_glyph
OUT = Path(__file__).resolve().parent / "out"; OUT.mkdir(parents=True, exist_ok=True)
FS, FSS, FS_LET, Z, ALPHA = 7.0, 6.5, 9, 1.96, .05
TAU = r"$\tau$*"
W3 = R / "e22_benchmark14" / "results"
B14 = ["EA", "N_SCZ", "N_BIP", "NEURO", "N_BMI", "N_AD", "MONO", "NEUT", "IBD", "RA", "N_LDL", "TC", "N_HEIGHT", "BMD"]
GROUPS = [("CNS_tissue", 5), ("myeloid", 3), ("lymphoid", 2), ("liver", 2), ("stromal", 2)]
LAYER = {t: ("DNASE" if t in ("N_HEIGHT", "BMD") else "CAGE") for t in B14}
NAME = {t: TN[t[2:] if t.startswith("N_") else t] for t in B14}
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n); m = 1.0
    for r in range(n - 1, -1, -1): m = min(m, p[o[r]] * n / (r + 1)); q[o[r]] = m
    return q

# ------------------------------------------------------------------ 자료
ME = pd.read_csv(W3 / "w04_measured.csv").set_index("trait").loc[B14]; AC = pd.read_csv(W3 / "w04_activity.csv").set_index("trait").loc[B14]; AB = pd.read_csv(W3 / "w06_actbin.csv").set_index("trait").loc[B14]
M72 = pd.read_csv(W3 / "w04_map72.csv"); S3 = pd.read_csv(R / "e18_cell_type_eqtl" / "results" / "s03_values.csv")
grp = [g for g, k in GROUPS for _ in range(k)]
assert list(ME.group) == grp and list(ME.layer) == [LAYER[t] for t in B14] and list(AC.layer) == list(ME.layer) and list(AB.layer) == list(ME.layer)
PA = []
for t in B14:
    r = M72[(M72.trait == t) & (M72.annot == f"L5_{ME.group[t]}_{LAYER[t]}")]; assert len(r) == 1 and abs(r.tau.iat[0] - ME.alone[t]) < 1e-9, t
    PA.append(dict(o=float(ME.alone[t]), o_se=float(r.se.iat[0]), f=float(ME.cond[t]), f_se=float(ME.cond_se[t]), star=bool(ME.q[t] < ALPHA), q=float(ME.q[t]), retained=float(ME.retained[t])))
def s3(t, model): r = S3[(S3.trait == t) & (S3.model == model) & (S3.kind == "E") & (S3.layer == LAYER[t])]; assert len(r) == 1, (t, model); return float(r.tau.iat[0]), float(r.se.iat[0])
EQ = pd.DataFrame([dict(trait=t, name=NAME[t], group=ME.group[t], layer=LAYER[t], E=s3(t, "E")[0], E_se=s3(t, "E")[1], EQ=s3(t, "EQ")[0], EQ_se=s3(t, "EQ")[1]) for t in B14]).set_index("trait")
EQ["retention"] = EQ.EQ / EQ.E; EQ["p1"] = stats.norm.sf(EQ.EQ / EQ.EQ_se); EQ["q"] = bh(EQ.p1); EQ.to_csv(OUT / "eqtl14.csv")
assert np.abs(EQ.E.to_numpy() - AC.effect_only.to_numpy()).max() < 1e-9 and np.abs(AB.effect_only.to_numpy() - AC.effect_only.to_numpy()).max() < 1e-12      # 세 표의 "효과 주석만" 값이 같은 모형이다
assert np.abs(AB.binary.to_numpy() - AC.effect.to_numpy()).max() < 1e-12
PB = [dict(o=float(EQ.E[t]), o_se=float(EQ.E_se[t]), f=float(EQ.EQ[t]), f_se=float(EQ.EQ_se[t]), star=bool(EQ.q[t] < ALPHA), q=float(EQ.q[t]), retained=float(EQ.retention[t])) for t in B14]
PC = [dict(o=float(AC.effect_only[t]), o_se=float(EQ.E_se[t]), m=float(AC.effect[t]), m_se=float(AC.effect_se[t]), m_q=float(AC.q[t]), f=float(AB.bins[t]), f_se=float(AB.bins_se[t]), star=bool(AB.q[t] < ALPHA), q=float(AB.q[t]),
           retained_top5=float(AC.retained[t]), retained_bins=float(AB.retained[t])) for t in B14]
PD_ = [dict(o=float(AC.low[t]), o_se=float(AC.low_se[t]), f=float(AC.high[t]), f_se=float(AC.high_se[t]), star=bool(AC.delta_q[t] < ALPHA), q=float(AC.delta_q[t]), delta=float(AC.delta[t]), delta_se=float(AC.delta_se[t])) for t in B14]
STARS = {k: [t for t, p in zip(B14, P) if p["star"]] for k, P in zip("abcd", (PA, PB, PC, PD_))}
STARS["c_top5"] = [t for t, p in zip(B14, PC) if p["m_q"] < ALPHA]
assert len(STARS["a"]) == 14 and STARS["c"] == ["TC", "N_HEIGHT", "BMD"] and STARS["d"] == ["N_HEIGHT", "BMD"] and len(STARS["c_top5"]) == 8      # e22_benchmark14 w04_summary.txt, w06_summary.txt와 같은 표시

# ------------------------------------------------------------------ 배치 (인치)
W = 7.2
XL, MR, GAP = 1.30, .08, .16                                   # 첫 패널의 왼쪽 끝(형질 이름 칸), 오른쪽 여백, 패널 사이
PW = (W - XL - MR - 3 * GAP) / 4
PITCH = .245; n = len(B14); PH = n * PITCH
MB, MT = .46, .64
H = MB + PH + MT
XLIM = (-2.9, 6.6); TICKS = [-2, 0, 2, 4, 6]                   # 네 패널이 같은 가로 범위
fig = plt.figure(figsize=(W, H))
ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(0, W); ov.set_ylim(0, H); ov.axis("off"); ov.set_zorder(5)
def tint(c, a=.45): r, g, b = to_rgb(c); return to_hex((1 - a + a * r, 1 - a + a * g, 1 - a + a * b))
def MK_O(c): return dict(marker="o", ms=3.8 - .9, mfc="white", mec=c, mew=.9, ls="none")
def MK_M(c): return dict(marker="o", ms=3.8 - .9, mfc=tint(c), mec=c, mew=.9, ls="none")
def MK_F(c): return dict(marker="o", ms=3.8, mfc=c, mec="white", mew=.5, ls="none")
PANELS = [("a", PA, [("o", MK_O, "predicted annotation alone"), ("f", MK_F, "with measured annotations")]),
          ("b", PB, [("o", MK_O, "effect annotations only"), ("f", MK_F, "with cell-type eQTLs")]),
          ("c", PC, [("o", MK_O, "effect annotations only"), ("m", MK_M, "with top 5% of activity"), ("f", MK_F, "with four activity bins")]),
          ("d", PD_, [("o", MK_O, "smaller predicted effects"), ("f", MK_F, "larger predicted effects")])]
edges = np.concatenate([[0], np.cumsum([k for _, k in GROUPS])])
for k, (let, P, marks) in enumerate(PANELS):
    x0 = XL + k * (PW + GAP)
    ax = fig.add_axes([x0 / W, MB / H, PW / W, PH / H])
    for j in range(0, len(GROUPS), 2): ax.axhspan(edges[j] - .5, edges[j + 1] - .5, color="#f3f3f1", lw=0, zorder=0)       # 조직 묶음을 하나 걸러 옅은 띠
    ax.axvline(0, color=MUTED, lw=.6, zorder=1)
    offs = [-.19, .19] if len(marks) == 2 else [-.27, 0.0, .27]
    for i, (t, p) in enumerate(zip(B14, P)):
        c = COL[LAYER[t]]
        for (key, sty, _), dy in zip(marks, offs):
            tau, se = p[key], p[key + "_se"]; assert XLIM[0] < tau - Z * se and tau + Z * se < XLIM[1] - .5, (let, t, key, tau, se)      # 모든 95% CI가 축 안, 별표 칸과 겹치지 않음
            ax.plot([tau - Z * se, tau + Z * se], [i + dy, i + dy], color=c, lw=.9, solid_capstyle="butt", zorder=2); ax.plot([tau], [i + dy], zorder=3, **sty(c))
        if p["star"]: ax.text(XLIM[1] - .12, i + .12, "*", ha="right", va="center", fontsize=FS, color=INK)
    ax.set_xlim(*XLIM); ax.set_ylim(n - .5, -.5); ax.set_xticks(TICKS); ax.set_yticks([])
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(MUTED); ax.spines["bottom"].set_linewidth(.6)
    ax.tick_params(labelsize=FS, length=2.5, width=.6, color=MUTED, labelcolor=INK, pad=2)
    ax.set_xlabel(TAU + " (95% CI)", fontsize=FS, color=INK, labelpad=2, math_fontfamily="stixsans")
    if k == 0:
        for i, t in enumerate(B14): ax.text(-.035, i, NAME[t], transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=FS, color=INK)
    ov.text(x0 - .14, H - .04, let, fontsize=FS_LET, fontweight="bold", ha="left", va="top", color=INK)
    for j, (_, sty, s) in enumerate(marks):                                                                  # 기호 풀이(패널 위; 표지는 먹색 — 색은 층위를 뜻한다)
        yk = H - .13 - (j + 3 - len(marks)) * .15
        ov.plot([x0 + .10], [yk], clip_on=False, **sty(INK)); ov.text(x0 + .19, yk, s, fontsize=FSS, color=INK, ha="left", va="center")
# 색 풀이: 형질마다의 층위 (이름 칸 아래)
yk = MB - .30
ov.text(.06, yk + .13, "Predicted annotation", fontsize=FSS, color=MUTED, ha="left", va="center")
x = .06
for lay, s in (("CAGE", "CAGE"), ("DNASE", "DNase")):
    ov.plot([x + .05], [yk], clip_on=False, **MK_F(COL[lay])); t_ = ov.text(x + .13, yk, s, fontsize=FSS, color=INK, ha="left", va="center")
    fig.canvas.draw(); x += .13 + t_.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi + .14

# ------------------------------------------------------------------ 점검 후 저장
fig.canvas.draw(); rd = fig.canvas.get_renderer(); fb = fig.bbox
ts = [t_ for t_ in fig.findobj(matplotlib.text.Text) if t_.get_visible() and t_.get_text().strip()]
ext = [(t_, matplotlib.text.Text.get_window_extent(t_, rd)) for t_ in ts]
assert not [t_.get_text() for t_, e in ext if e.x0 < fb.x0 or e.x1 > fb.x1 or e.y0 < fb.y0 or e.y1 > fb.y1]
assert all(6.5 <= t_.get_fontsize() <= 7 or (t_.get_fontsize() == FS_LET and t_.get_fontweight() == "bold") for t_ in ts)
assert all(to_hex(t_.get_color()) in (INK, MUTED) for t_ in ts)
ovl = [(a.get_text(), b.get_text()) for i, (a, x_) in enumerate(ext) for b, y_ in ext[i + 1:] if x_.overlaps(y_)]
assert not ovl, ovl
fig.savefig(OUT / "fig5.png", dpi=300); fig.savefig(OUT / "fig5.pdf"); plt.close(fig)
def rows(P): return {t: p for t, p in zip(B14, P)}
NUM = {"form": "rows = 14 benchmark traits (shared by the four panels; light bands separate the tissue groups); open above, filled below; horizontal interval = 95% CI; asterisk at the right edge; marker colour = layer of the annotation",
       "row_order": B14, "names": NAME, "tissue_group": {t: LN[ME.group[t]] for t in B14}, "layer": LAYER, "colours": {k: COL[k] for k in ("CAGE", "DNASE")}, "x_range": list(XLIM),
       "a": rows(PA), "a_keys": "o = predicted annotation alone (single-annotation model with baselineLD; w04_map72), f = with the measured CAGE, DNase and H3K27ac annotations of the same tissue group (w04_measured); asterisk: BH FDR < 0.05 across 14 for f > 0",
       "a_retained_median": float(np.median(ME.retained)), "a_retained_range": [float(ME.retained.min()), float(ME.retained.max())],
       "b": rows(PB), "b_keys": "o = effect annotation in the joint model of the three effect layers, all tracks (model E), f = with the five cell-type eQTL annotations added (model EQ; e18_cell_type_eqtl s03_values.csv); asterisk: BH FDR < 0.05 across 14 for f > 0 (one-sided P from tau/se)",
       "b_retained_range": [float(EQ.retention.min()), float(EQ.retention.max())], "b_retained_median": float(EQ.retention.median()), "b_max_p1": float(EQ.p1.max()), "b_n_fdr": int((EQ.q < ALPHA).sum()),
       "c": rows(PC), "c_keys": "o = effect annotations only (model E), m = with the top-5% predicted-activity annotations of the three layers (w04_activity: effect), f = with four predicted-activity bins per layer (w06_actbin: bins); asterisk: BH FDR < 0.05 across 14 for f > 0",
       "d": rows(PD_), "d_keys": "within the top 10% of predicted activity: o = activity-matched half with smaller predicted effects, f = half with larger predicted effects (w04_activity: low, high); asterisk: BH FDR < 0.05 across 14 for delta = f - o > 0",
       "asterisks": STARS, "asterisk_names": {k: [NAME[t] for t in v] for k, v in STARS.items()}}
(OUT / "numbers_fig5.json").write_text(json.dumps(NUM, indent=1, ensure_ascii=False), encoding="utf-8")
print(f"saved {OUT / 'fig5.png'} ({W} x {H:.2f} in); asterisks { {k: len(v) for k, v in STARS.items()} }")
print("eQTL retention:", {NAME[t]: round(float(EQ.retention[t]), 3) for t in B14}, "| max one-sided P", f"{EQ.p1.max():.2g}", "| FDR < 0.05:", int((EQ.q < ALPHA).sum()))
