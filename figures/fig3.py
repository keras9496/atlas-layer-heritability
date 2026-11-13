# -*- coding: utf-8 -*-
"""Figure 3 — 벤치마크 14형질의 층위.
 a 형질마다 범주의 조직 묶음에서 세 층위(CAGE, DNase, H3K27ac)를 함께 적합한 τ* ± 95% CI. 층위당 세 단위 재추출 20회 평균(w04_layers.csv).
   형질 순서: brain tissue(EA, SCZ, BIP, Neuroticism, BMI) → myeloid(AD, Monocytes, Neutrophils) → lymphoid(IBD, RA) → liver(LDL, Total chol.) → stromal(Height, BMD).
 b 4분원: 각도 = CAGE–DNase 균형 = atan2(DNase τ*, CAGE τ*) (w04_layers.csv의 phiD, phiD_se). 0° = CAGE만, 90° = DNase만, 45° 점선 = 같음.
   고리 하나 = 형질 하나(안쪽 = a의 첫 형질). 옅은 띠 = 같은 조직 묶음의 형질(하나 걸러). 중심에서의 거리는 값이 아니다. 호 = 95% CI(±1.96 s.e.).
 c 형질마다 CAGE − DNase 대비 ± 95% CI: 채운 점 = 층위당 세 단위(a와 같은 적합), 빈 점 = 세포형을 맞춘 층위(w04_controls.csv, model == matched).
 d 같은 대비를 TSS에서 1 kb 밖의 SNP(채운 점, model == distal)와 1 kb 안의 SNP(빈 점, model == proximal)에서. 전체 트랙 주석을 TSS 거리로 나눈 모형이다.
   c와 d는 줄(형질)과 가로 눈금이 같다. 옅은 띠는 b와 같은 조직 묶음 구분이다.
 자료: e22_benchmark14/results/w04_layers.csv, w04_controls.csv, w04_summary.txt. 기존 결과 파일만 읽는다. 새 통계 없음.
실행: python figures/fig3.py      (저장소 뿌리에서; 모두 한 번에는 python figures/make_figures.py)
출력: figures/out/fig3.png, fig3.pdf, numbers_fig3.json
"""
import sys, json, math, re
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))      # -I 로 실행하면 스크립트 폴더가 sys.path에 없다
import numpy as np, pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Wedge, Rectangle
from common import R as RESROOT, COL, SHORT, LAY, LN, INK, MUTED, GRID

plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "legend.fontsize": 6.5, "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6,
                     "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.color": INK, "ytick.color": INK,
                     "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
                     "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
                     "pdf.fonttype": 42, "axes.unicode_minus": True})
from matplotlib import _mathtext as _mt
def _xheight_from_glyph(self, fontname, fontsize, dpi):
    """matplotlib는 PCLT 표가 있는 글꼴(Arial)의 x-height를 약 1.9배 크게 계산한다 → 위·아래첨자가 어긋난다. 'x' 글자 높이를 재서 쓴다."""
    return self.get_metrics(fontname, matplotlib.rcParams["mathtext.default"], "x", fontsize, dpi).iceberg
_mt.TruetypeFonts.get_xheight = _xheight_from_glyph
FS, FS_NOTE, FS_LET = 7, 6.5, 9
PT = 1 / 72                                                   # 1 pt (인치)
Z = 1.96
OUT = Path(__file__).resolve().parent / "out"; OUT.mkdir(parents=True, exist_ok=True)

def chk(x, target, nd=None):
    """기록(w04_summary.txt)에 적힌 수치(target, 소수 nd자리)를 값 x가 재현하는지 확인."""
    nd = len(str(target).split(".")[1]) if nd is None and "." in str(target) else (nd or 0)
    assert abs(float(x) - target) <= .5 * 10 ** -nd + 1e-9, (x, target)

# ------------------------------------------------------------------ 자료와 기록 대조
W4 = RESROOT / "e22_benchmark14" / "results"
ORD = ["EA", "N_SCZ", "N_BIP", "NEURO", "N_BMI", "N_AD", "MONO", "NEUT", "IBD", "RA", "N_LDL", "TC", "N_HEIGHT", "BMD"]      # 벤치마크 14형질
ABB = {"EA": "EA", "N_SCZ": "SCZ", "N_BIP": "BIP", "NEURO": "Neuroticism", "N_BMI": "BMI", "N_AD": "AD", "MONO": "Monocytes", "NEUT": "Neutrophils", "IBD": "IBD", "RA": "RA",
       "N_LDL": "LDL", "TC": "Total chol.", "N_HEIGHT": "Height", "BMD": "BMD"}
TN = {"EA": "Educational attainment", "N_SCZ": "Schizophrenia", "N_BIP": "Bipolar disorder", "NEURO": "Neuroticism", "N_BMI": "BMI", "N_AD": "Alzheimer's disease", "MONO": "Monocyte count",
      "NEUT": "Neutrophil count", "IBD": "IBD", "RA": "Rheumatoid arthritis", "N_LDL": "LDL cholesterol", "TC": "Total cholesterol", "N_HEIGHT": "Height", "BMD": "Bone mineral density"}
LY = pd.read_csv(W4 / "w04_layers.csv").set_index("trait"); assert set(LY.index) == set(ORD) and LY.index.is_unique
CT = pd.read_csv(W4 / "w04_controls.csv"); assert len(CT) == 42 and set(CT.model) == {"matched", "distal", "proximal"} and CT.groupby(["trait", "model"]).size().max() == 1
CTL = {m: CT[CT.model == m].set_index("trait").loc[ORD] for m in ("matched", "distal", "proximal")}
n = len(ORD); assert n == 14
KEY = {"CAGE": "cage", "DNASE": "dnase", "H3K27ac": "h3k27ac"}
for t in ORD:                                                                                   # 열의 정의 확인(계산해서 그리지는 않는다)
    r = LY.loc[t]; assert abs(r.c1 - (r.cage - r.dnase)) < 1e-9 and abs(r.phiD - math.degrees(math.atan2(r.dnase, r.cage))) < 1e-6 and bool(r.phiD_valid), t
    assert [m for m in LAY if r[KEY[m]] == max(r.cage, r.dnase, r.h3k27ac)][0] == r.top, t
e, h, b, ad = LY.loc["EA"], LY.loc["N_HEIGHT"], LY.loc["BMD"], LY.loc["N_AD"]                    # w04_summary.txt 3절
chk(e.cage, 1.13); chk(e.cage_se, .21); chk(e.dnase, -.06); chk(e.dnase_se, .18); chk(e.h3k27ac, .26)
chk(h.cage, .88); chk(h.dnase, 1.48); chk(h.c1, -.61); chk(h.c1_se, .38); chk(h.phiD, 59.4); chk(h.phiD_se, 9.2); assert int(h.n_first) == 16 and h.top == "DNASE"
chk(b.c1, -1.21); chk(b.c1_se, .45); chk(b.phiD, 69.6); chk(b.phiD_se, 8.0); assert int(b.n_first) == 20 and b.top == "DNASE"
chk(ad.c1, 3.76); chk(ad.c1_se, 1.37); chk(ad.phiD, -8.3); chk(ad.phiD_se, 10.9)
assert all(LY.top[t] == "CAGE" for t in ORD if t not in ("N_HEIGHT", "BMD"))                    # 나머지 12형질은 CAGE의 τ*가 가장 크다
chk(CTL["matched"].c1["N_HEIGHT"], -1.03); chk(CTL["matched"].c1_se["N_HEIGHT"], .46); chk(CTL["matched"].c1["BMD"], -2.23); chk(CTL["matched"].c1_se["BMD"], .55)   # 5절
chk(CTL["distal"].c1["N_HEIGHT"], -1.47); chk(CTL["distal"].c1_se["N_HEIGHT"], .43); chk(CTL["proximal"].c1["N_HEIGHT"], -.37); chk(CTL["proximal"].c1_se["N_HEIGHT"], .66)
chk(CTL["distal"].c1["N_LDL"], 3.60); chk(CTL["proximal"].c1["N_LDL"], .01); chk(CTL["matched"].c1["N_AD"], 3.17); chk(CTL["proximal"].c1["MONO"], 3.87)
SAME = {m: int(sum((CTL[m].c1[t] > 0) == (LY.c1[t] > 0) for t in ORD)) for m in CTL}; assert SAME == {"matched": 14, "distal": 14, "proximal": 12}
lin = LY.loc[ORD, "group"].tolist(); groups = []
for i, l in enumerate(lin):
    if groups and groups[-1][0] == l: groups[-1][2] = i
    else: groups.append([l, i, i])
assert [(g[0], g[2] - g[1] + 1) for g in groups] == [("CNS_tissue", 5), ("myeloid", 3), ("lymphoid", 2), ("liver", 2), ("stromal", 2)]
def ci(r): return float(r.phiD - Z * r.phiD_se), float(r.phiD + Z * r.phiD_se)
# 그림 설명에 쓸 검정값(그림에는 넣지 않는다): w04_summary.txt의 줄을 그대로 읽는다
SUMTXT = (W4 / "w04_summary.txt").read_text(encoding="utf-8")
def grab(pat):
    m = re.search(pat, SUMTXT); assert m, pat; return m
STAT = {}
m = grab(r"이질성 \[14형질; \(C−D, C−H\)\]: Q = ([\d.]+), 자유도 (\d+), P = ([\d.e+-]+) \(F 근사 ([\d.e+-]+)\)"); STAT["layer_profile_Q_14"] = dict(Q=float(m[1]), df=int(m[2]), P=float(m[3]), P_F=float(m[4]))
m = grab(r"이질성 \[14형질; φD\]: Q = ([\d.]+), 자유도 (\d+), P = ([\d.e+-]+) \(F 근사 ([\d.e+-]+)\); 공통 φD ([\d.−-]+)°"); STAT["balance_Q_14"] = dict(Q=float(m[1]), df=int(m[2]), P=float(m[3]), P_F=float(m[4]), common_deg=float(m[5]))
m = grab(r"이질성 \[골격 둘을 뺀 12형질; φD\]: Q = ([\d.]+), 자유도 (\d+), P = ([\d.e+-]+) \(F 근사 ([\d.e+-]+)\); 공통 φD ([\d.−-]+)°"); STAT["balance_Q_12_without_height_BMD"] = dict(Q=float(m[1]), df=int(m[2]), P=float(m[3]), P_F=float(m[4]), common_deg=float(m[5]))
m = grab(r"대비 \[골격 − 나머지 12; φD\]: ([\d.]+) ± ([\d.]+) 대 ([\d.]+) ± ([\d.]+), 차이 \+([\d.]+) ± ([\d.]+), 양측 P = ([\d.e+-]+)"); STAT["balance_height_BMD_vs_other_12"] = dict(skeletal=float(m[1]), skeletal_se=float(m[2]), other=float(m[3]), other_se=float(m[4]), diff=float(m[5]), diff_se=float(m[6]), P=float(m[7]))
for lab in ("matched", "distal", "proximal"):
    m = grab(rf"\[{lab}\] Q = ([\d.]+) \(자유도 (\d+), P = ([\d.e+-]+)\); C−D의 부호가 재추출 지도와 같은 형질 (\d+)/14"); STAT[f"layer_profile_Q_{lab}"] = dict(Q=float(m[1]), df=int(m[2]), P=float(m[3]), same_sign=int(m[4]))
    assert STAT[f"layer_profile_Q_{lab}"]["same_sign"] == SAME[lab]
chk(STAT["layer_profile_Q_14"]["Q"], 104.3); chk(STAT["balance_Q_14"]["Q"], 107.3); chk(STAT["balance_Q_12_without_height_BMD"]["Q"], 13.7); chk(STAT["layer_profile_Q_proximal"]["Q"], 22.1)
RULE = dict(met=[t for t in ORD if bool(LY.resolved[t])], not_met=[t for t in ORD if not bool(LY.resolved[t])]); assert RULE["not_met"] == ["NEURO", "NEUT", "IBD", "N_HEIGHT"]

# ------------------------------------------------------------------ 배치 (인치): 위 블록(a, b)과 아래 블록(c, d)
W = 7.2
RQ, RIN = 2.30, .46                                # b: 4분원 반지름, 가장 안쪽 줄의 안쪽 끝
LW_ = (RQ - RIN) / n                               # 줄 하나의 폭
RC = {t: RIN + (i + .5) * LW_ for i, t in enumerate(ORD)}
ARCS = [(t, "main", RC[t], *ci(LY.loc[t])) for t in ORD]
YM = max(-r * math.sin(math.radians(lo)) for _, _, r, lo, hi in ARCS) + .07      # 가로축 아래 여백: 가장 낮은 95% CI 끝 + 0.07 in
XM = max(max(-r * math.cos(math.radians(hi)) for _, _, r, lo, hi in ARCS) + .07, .12)   # 세로축 왼쪽 여백(90°를 넘는 95% CI가 없으면 0.12 in)
assert .3 < YM < 1.2 and .1 < XM < .6, (YM, XM)
TOPPAD, BOTPAD, RIGHTPAD = .44, .05, .60           # 테두리 위(90°, DNase), 아래, 오른쪽(0°, CAGE)
HT = TOPPAD + RQ + YM + BOTPAD                     # 위 블록의 높이
ROWH, C_TOP, C_BELOW = .17, .40, .42               # 아래 블록: 줄 높이, 축 위(패널 글자·풀이), 축 아래(눈금·축 이름)
HBOT = C_TOP + n * ROWH + C_BELOW
H = HT + HBOT
OX, OY = W - RIGHTPAD - RQ, HBOT + BOTPAD + YM
fig = plt.figure(figsize=(W, H)); RD = fig.canvas.get_renderer()
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")          # b와 패널 글자: 자료 좌표 = 인치
def T(x, y, s, **kw):
    kw.setdefault("fontsize", FS_NOTE); kw.setdefault("color", INK); kw.setdefault("ha", "left"); kw.setdefault("va", "center")
    if "$" in s: kw.setdefault("math_fontfamily", "stixsans")
    return ax.text(x, y, s, **kw)
def tsize(s, **kw):
    """글자의 폭·높이(인치)."""
    t = T(0, 0, s, **kw); e_ = t.get_window_extent(RD); t.remove(); return e_.width / fig.dpi, e_.height / fig.dpi
def mk(x, y, z=3, **kw): ax.plot([x], [y], ls="none", zorder=z, clip_on=False, **kw)
def letter(x, y, s): ax.text(x, y, s, fontsize=FS_LET, fontweight="bold", ha="left", va="top", color=INK)
def axin(l, b_, w, h_): return fig.add_axes([l / W, b_ / H, w / W, h_ / H])
def clean(a): a.spines[["top", "right"]].set_visible(False); a.tick_params(length=2.5, color=MUTED, width=.6)
NUM = {}

# ------------------------------------------------------------------ b 4분원 (Quad·canvas·ci_arc)
BAND, SPOKE, LANE = "#f3f3f1", "#e2e2e2", "#d6d6d6"            # 조직 묶음 띠의 옅은 바탕(하나 걸러), 15° 간격 살, 줄 안내선
class Quad:
    """4분원 하나. 원점 (ox, oy)와 길이는 그림 인치. 각도 0° = 가로축(CAGE), 90° = 세로축(DNase).
    그리는 영역은 L자 모양: 4분원 + 세로축 왼쪽 띠(x ≥ ox − xm, y ≥ oy) + 가로축 아래 띠(x ≥ ox, y ≥ oy − ym)."""
    def __init__(s, ox, oy, R, rin, xm, ym):
        s.ox, s.oy, s.R, s.rin, s.xm, s.ym = ox, oy, R, rin, xm, ym
        e_ = R + .02
        s.shape = Polygon([(ox - xm, oy), (ox - xm, oy + e_), (ox + e_, oy + e_), (ox + e_, oy - ym), (ox, oy - ym), (ox, oy)], closed=True, transform=ax.transData)
    def xy(s, r, deg): a = math.radians(deg); return s.ox + r * math.cos(a), s.oy + r * math.sin(a)
    def arc(s, r, a0, a1):
        a = np.radians(np.linspace(a0, a1, max(2, int(abs(a1 - a0) / .4) + 1))); return s.ox + r * np.cos(a), s.oy + r * np.sin(a)
    def clip(s, artists):
        for a in artists: a.set_clip_path(s.shape)
    def left_y(s, r): return s.oy + math.sqrt(r * r - s.xm * s.xm)                       # 줄 가운데 선이 왼쪽 끝과 만나는 높이
    def inside(s, r, deg, pad=0.0):
        x, y = s.xy(r, deg)
        return x - pad > s.ox - s.xm and y - pad > s.oy - s.ym and r + pad < s.R and not (x - pad < s.ox and y - pad < s.oy)

def canvas(q, bands, lanes, ticklabs=(0, 30, 60, 90)):
    """조직 묶음 띠(하나 걸러 옅은 바탕), 줄 안내선, 테두리 호, 15° 간격 살, 45° 점선, 두 축(층위 색), 각도 눈금과 축 끝 이름."""
    for r0, r1 in bands:
        w = ax.add_patch(Wedge((q.ox, q.oy), r1, -90, 180, width=r1 - r0, fc=BAND, ec="none", zorder=.5)); w.set_clip_path(q.shape)
    for r in lanes: q.clip(ax.plot(*q.arc(r, -90, 180), color=LANE, lw=.4, zorder=.9))
    q.clip(ax.plot(*q.arc(q.R, -90, 180), color=MUTED, lw=.6, zorder=1.3))
    for a in (15, 30, 60, 75): ax.plot(*zip(q.xy(q.rin, a), q.xy(q.R, a)), color=SPOKE, lw=.5, zorder=1, solid_capstyle="butt")
    ax.plot(*zip(q.xy(.10, 45), q.xy(q.R, 45)), color=GRID, lw=.7, ls=(0, (3, 2)), zorder=1.1)
    ax.plot([q.ox, q.ox + q.R + 5 * PT], [q.oy, q.oy], color=COL["CAGE"], lw=1.1, zorder=2, solid_capstyle="butt")
    ax.plot([q.ox, q.ox], [q.oy, q.oy + q.R + 5 * PT], color=COL["DNASE"], lw=1.1, zorder=2, solid_capstyle="butt")
    for a in range(15, 90, 15): ax.plot(*zip(q.xy(q.R, a), q.xy(q.R + 2.5 * PT, a)), color=MUTED, lw=.6, zorder=1.3, solid_capstyle="butt")
    for a in ticklabs:
        x, y = q.xy(q.R + (12.5 if a in (0, 90) else 9) * PT, a); T(x, y, f"{a}°", ha="center", va="center")
    x, y = q.xy(q.rin * .58, 45)
    T(x, y, "equal", ha="center", va="center", rotation=45, color=MUTED, bbox=dict(fc="white", ec="none", pad=.8), zorder=1.2)
    T(q.ox, q.oy + q.R + 23 * PT, "DNase", ha="center", va="center", fontweight="bold")           # 축 끝의 이름: 세로축 = DNase, 가로축 = CAGE
    T(q.ox + q.R + 21.5 * PT, q.oy, "CAGE", ha="left", va="center", fontweight="bold")
def ci_arc(q, r, lo, hi, color, lw=.9, z=2.6):
    assert q.inside(r, lo) and q.inside(r, hi), ("95% CI가 그리는 영역 밖", r, lo, hi)
    ax.plot(*q.arc(r, lo, hi), color=color, lw=lw, zorder=z, solid_capstyle="butt")

QB = Quad(OX, OY, RQ, RIN, XM, YM)
EDGE = {g[0]: (RIN + g[1] * LW_, RIN + (g[2] + 1) * LW_) for g in groups}                    # 조직 묶음마다 띠의 안쪽·바깥쪽 반지름
canvas(QB, [EDGE[g[0]] for k, g in enumerate(groups) if k % 2 == 0], [RC[t] for t in ORD])
FILL = dict(marker="o", ms=4.0, mfc=INK, mec="white", mew=.5)
DOTS = []                                                                                     # (x pt, y pt, 반지름 pt, 이름): 겹침 점검용
for t in ORD:
    r = LY.loc[t]; lo, hi = ci(r); rr = RC[t]
    ci_arc(QB, rr, lo, hi, INK); mk(*QB.xy(rr, r.phiD), z=3.2, **FILL)
    DOTS.append((rr * 72 * math.cos(math.radians(r.phiD)), rr * 72 * math.sin(math.radians(r.phiD)), FILL["ms"] / 2, t))
XL = QB.ox - QB.xm - .05                                                                      # 줄 이름의 오른쪽 끝
for t in ORD: T(XL, QB.left_y(RC[t]), ABB[t], ha="right")
LABW = max(tsize(ABB[t])[0] for t in ORD)
x, y = QB.xy(QB.R + 15 * PT, 45)
T(x, y, "CAGE–DNase balance", ha="center", va="center", rotation=-45, fontsize=FS)
NUM["b_definition"] = ("angle = CAGE-DNase balance = atan2(DNase tau*, CAGE tau*) in degrees; 0 = CAGE only, 90 = DNase only, 45 = equal; one ring per trait, inner to outer in the order of panel a; "
                       "light bands group traits analysed in the same tissue group; distance from the centre carries no value; arcs = estimate ± 1.96 × s.e.")
NUM["b_source"] = {"file": "e22_benchmark14/results/w04_layers.csv", "columns": ["phiD", "phiD_se"]}
NUM["b_ring_order_inner_to_outer"] = [dict(key=t, label=ABB[t], name=TN[t]) for t in ORD]
NUM["b_bands"] = [dict(tissue_group=g[0], name=LN[g[0]], traits=ORD[g[1]:g[2] + 1], shaded=(k % 2 == 0)) for k, g in enumerate(groups)]
NUM["b_balance"] = {t: dict(deg=float(LY.phiD[t]), se=float(LY.phiD_se[t]), ci95=list(ci(LY.loc[t])), cage_tau=float(LY.cage[t]), dnase_tau=float(LY.dnase[t])) for t in ORD}
NUM["b_margins_in"] = dict(below_axis=YM, left_of_axis=XM, lowest_ci_deg=min(a[3] for a in ARCS), highest_ci_deg=max(a[4] for a in ARCS))

# ------------------------------------------------------------------ a 세 층위의 τ* (14형질)
X0 = .44; X1 = XL - LABW - .13                    # a 축의 왼쪽·오른쪽 끝 (b의 줄 이름 칸 왼쪽까지)
A_TOP, A_BELOW = .30, .86                         # a 위(패널 글자·범례), a 아래(기울인 눈금 글자, 묶음 괄호와 이름 두 줄)
A_H = HT - A_TOP - A_BELOW
aA = axin(X0, HBOT + A_BELOW, X1 - X0, A_H)
xs = np.arange(n); off = {"CAGE": -.26, "DNASE": 0, "H3K27ac": .26}
for m_ in LAY:
    aA.errorbar(xs + off[m_], LY.loc[ORD, KEY[m_]], yerr=Z * LY.loc[ORD, KEY[m_] + "_se"], fmt="o", ms=3.3, color=COL[m_], mec="white", mew=.45,
                ecolor=COL[m_], elinewidth=.85, capsize=0, zorder=3)
aA.axhline(0, color=MUTED, lw=.6, zorder=1)
for g in groups[1:]: aA.axvline(g[1] - .5, color=GRID, lw=.6, zorder=0)
aA.set_xlim(-.6, n - .4); aA.set_xticks(xs)
aA.set_xticklabels([ABB[t] for t in ORD], fontsize=FS_NOTE, rotation=40, ha="right", rotation_mode="anchor")
_lo = min((LY.loc[ORD, KEY[m_]] - Z * LY.loc[ORD, KEY[m_] + "_se"]).min() for m_ in LAY); _hi = max((LY.loc[ORD, KEY[m_]] + Z * LY.loc[ORD, KEY[m_] + "_se"]).max() for m_ in LAY)
YLIM_A = (-2.3, 6.0); aA.set_ylim(*YLIM_A); aA.set_yticks([-2, 0, 2, 4, 6]); assert YLIM_A[0] < _lo and _hi < YLIM_A[1], (_lo, _hi)   # 모든 95% CI가 축 안
tr = aA.get_xaxis_transform(); BY = -.43 / A_H                                                # 묶음 괄호: 축 아래 0.43 in (기울인 눈금 글자 밑)
GLAB = {l: LN[l] for l in lin}; GLAB["stromal"] = LN["stromal"].replace("/", "/\n"); assert GLAB["stromal"] == "Stromal/\nmesenchymal"
for l, i0, i1 in groups:                                                                       # 이름은 괄호의 왼쪽 끝에 맞춘다
    aA.plot([i0 - .42, i1 + .42], [BY, BY], color=MUTED, lw=.6, transform=tr, clip_on=False)
    aA.text(i0 - .42, BY - .045 / A_H, GLAB[l], ha="left", va="top", fontsize=FS_NOTE, color=MUTED, transform=tr, linespacing=1.15)
aA.text(-.66, BY - .045 / A_H, "tissue\ngroup:", ha="right", va="top", fontsize=FS_NOTE, color=MUTED, transform=tr, linespacing=1.15)
aA.set_ylabel(r"$\tau$* (± 95% CI)", labelpad=2, math_fontfamily="stixsans"); clean(aA); aA.tick_params(axis="x", length=0, pad=2)
hs = [Line2D([], [], marker="o", ms=3.3, color=COL[m_], mec="white", mew=.45, ls="none") for m_ in LAY]
aA.legend(hs, [SHORT[m_] for m_ in LAY], loc="lower right", bbox_to_anchor=(1, 1.02), ncol=3, frameon=False, fontsize=FS_NOTE, handlelength=.9,
          handletextpad=.3, columnspacing=1.2, borderaxespad=0, borderpad=0)
NUM["a_layers"] = {t: {m_: dict(tau=float(LY.loc[t, KEY[m_]]), se=float(LY.loc[t, KEY[m_] + "_se"])) for m_ in LAY} | {"tissue_group": LY.group[t], "tissue_group_name": LN[LY.group[t]], "name": TN[t], "gwas": LY.gwas[t],
                       "largest_tau_layer": LY.top[t], "resamples_ranked_first_of_20": int(LY.n_first[t]), "largest_vs_second_one_sided_P": float(LY.top_vs_second_p1[t])} for t in ORD}
NUM["a_order"] = ORD; NUM["a_source"] = "e22_benchmark14/results/w04_layers.csv" + " (three layers fitted jointly in the tissue group of the trait's category; at most three units per layer, 20 resamples averaged)"

# ------------------------------------------------------------------ c, d: CAGE − DNase 대비 (줄 = 형질, a와 같은 순서로 위에서 아래)
CX0, CW, CGAP = 1.40, 2.66, .34                   # c 축의 왼쪽 끝, 축 폭, c와 d 사이
DX0 = CX0 + CW + CGAP; assert DX0 + CW < W - .08
CY0 = C_BELOW                                      # 축의 아래 끝
PAN = {"c": dict(x0=CX0, sets=[("three units", LY.loc[ORD], True, "Three units per layer"), ("matched", CTL["matched"], False, "Cell-type-matched")]),
       "d": dict(x0=DX0, sets=[("distal", CTL["distal"], True, "> 1 kb from TSS"), ("proximal", CTL["proximal"], False, "≤ 1 kb from TSS")])}
allv = [(float(df.c1[t] - Z * df.c1_se[t]), float(df.c1[t] + Z * df.c1_se[t])) for p in PAN.values() for _, df, _, _ in p["sets"] for t in ORD]
XLIM_C = (-4.0, 7.4); assert XLIM_C[0] < min(a for a, _ in allv) and max(b_ for _, b_ in allv) < XLIM_C[1], (min(a for a, _ in allv), max(b_ for _, b_ in allv))
DY = .20                                           # 한 줄 안에서 두 표지의 세로 간격의 절반(줄 단위)
FILLc = dict(marker="o", ms=3.9, mfc=INK, mec="white", mew=.5); OPENc = dict(marker="o", ms=3.9 - .8, mfc="white", mec=INK, mew=.8)
CAX = {}
for name, p in PAN.items():
    a = axin(p["x0"], CY0, CW, n * ROWH); CAX[name] = a
    for k, g in enumerate(groups):                                                             # b와 같은 띠(조직 묶음 하나 걸러)
        if k % 2 == 0: a.add_patch(Rectangle((XLIM_C[0], g[1] - .5), XLIM_C[1] - XLIM_C[0], g[2] - g[1] + 1, fc=BAND, ec="none", zorder=0))
    a.axvline(0, color=MUTED, lw=.6, zorder=1)
    for (key, df, filled, lab), sgn in zip(p["sets"], (-1, 1)):
        for i, t in enumerate(ORD):
            v, se = float(df.c1[t]), float(df.c1_se[t]); yy = i + sgn * DY
            a.plot([v - Z * se, v + Z * se], [yy, yy], color=INK, lw=.85, zorder=2, solid_capstyle="butt")
            a.plot([v], [yy], ls="none", zorder=3, **(FILLc if filled else OPENc))
    a.set_xlim(*XLIM_C); a.set_ylim(n - .5, -.5); a.set_xticks([-4, -2, 0, 2, 4, 6])
    a.set_xticklabels(["−4", "−2", "0", "+2", "+4", "+6"]); a.set_yticks([]); clean(a); a.spines["left"].set_visible(False)
    a.set_xlabel(r"CAGE − DNase ($\tau$*)", labelpad=2, math_fontfamily="stixsans")
    # 풀이: 축 위 한 줄
    kx = p["x0"] + .04; ky = CY0 + n * ROWH + .12
    for (key, df, filled, lab) in p["sets"]:
        mk(kx + .03, ky, **(FILLc if filled else OPENc)); T(kx + .11, ky, lab); kx += .11 + tsize(lab)[0] + .22
for i, t in enumerate(ORD): T(CX0 - .08, CY0 + (n - i - .5) * ROWH, TN[t], ha="right", fontsize=FS)
NUM["c_d_definition"] = "CAGE − DNase contrast of tau* (± 1.96 s.e.) in the tissue group of the trait's category; rows in the order of panel a; shaded bands = tissue groups as in b"
NUM["c"] = {t: {"three_units": dict(c1=float(LY.c1[t]), se=float(LY.c1_se[t])), "cell_type_matched": dict(c1=float(CTL["matched"].c1[t]), se=float(CTL["matched"].c1_se[t]))} for t in ORD}
NUM["d"] = {t: {"distal_gt_1kb_from_TSS": dict(c1=float(CTL["distal"].c1[t]), se=float(CTL["distal"].c1_se[t])), "proximal_le_1kb_from_TSS": dict(c1=float(CTL["proximal"].c1[t]), se=float(CTL["proximal"].c1_se[t]))} for t in ORD}
NUM["c_d_source"] = "e22_benchmark14/results/w04_controls.csv" + " (model = matched | distal | proximal); filled circles of c = w04_layers.csv (same fits as a)"
NUM["c_d_xlim"] = list(XLIM_C); NUM["c_d_same_sign_as_three_units_of_14"] = SAME

# ------------------------------------------------------------------ 패널 글자
letter(.02, H - .05, "a"); letter(X1 + .10, H - .05, "b"); letter(.02, HBOT - .03, "c"); letter(DX0 - .22, HBOT - .03, "d")

# ------------------------------------------------------------------ 그림 설명에 쓸 검정값 (그림에는 없음)
NUM["legend_statistics"] = STAT | {"rule_15_of_20_and_P_lt_0.05": {k: [TN[t] for t in v] for k, v in RULE.items()}, "source": "e22_benchmark14/results/w04_summary.txt"}

# ------------------------------------------------------------------ 겹침·잘림 점검 후 저장
def _box(t, rd):
    """글자의 상자 네 꼭짓점(화면 좌표). 돌린 글자는 돌린 사각형으로 잰다."""
    rot = t.get_rotation() % 360
    if rot % 90 == 0:
        e_ = matplotlib.text.Text.get_window_extent(t, rd); return np.array([[e_.x0, e_.y0], [e_.x1, e_.y0], [e_.x1, e_.y1], [e_.x0, e_.y1]])
    assert t.get_rotation_mode() == "anchor" or (t.get_ha() == "center" and t.get_va() == "center"), ("돌린 글자의 상자를 잴 수 없다", t.get_text())
    t.set_rotation(0); e_ = matplotlib.text.Text.get_window_extent(t, rd); t.set_rotation(rot)
    a = np.asarray(t.get_transform().transform(t.get_unitless_position()), float); th = math.radians(rot)
    M = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
    return (np.array([[e_.x0, e_.y0], [e_.x1, e_.y0], [e_.x1, e_.y1], [e_.x0, e_.y1]]) - a) @ M.T + a
def _hit(P, Q):
    """볼록 다각형 둘이 겹치는지(분리축)."""
    for poly in (P, Q):
        for i in range(len(poly)):
            e_ = poly[(i + 1) % len(poly)] - poly[i]; nrm = np.array([-e_[1], e_[0]]); a, b_ = P @ nrm, Q @ nrm
            if a.max() <= b_.min() + 1e-6 or b_.max() <= a.min() + 1e-6: return False
    return True
def text_issues(fig):
    fig.canvas.draw(); rd = fig.canvas.get_renderer(); fb = fig.bbox
    ts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    for t in ts:
        s = t.get_text(); assert "τ" not in s and "φ" not in s and "\t" not in s, repr(s)
        if r"\tau" in s: assert t.get_math_fontfamily() == "stixsans", s
        assert r"\varphi" not in s and "lineage" not in s.lower() and "direction" not in s.lower(), s    # 그림 글자에는 기호 φD와 코드의 낱말(lineage 등)을 쓰지 않는다
        assert matplotlib.colors.to_hex(t.get_color()) in (INK, MUTED), (s, t.get_color())                 # 글자 색은 INK/MUTED만
        assert 6.5 <= t.get_fontsize() <= 7 or (t.get_fontsize() == FS_LET and t.get_fontweight() == "bold"), (s, t.get_fontsize())
    bb = [(t, _box(t, rd)) for t in ts]
    ov = [(a.get_text(), b_.get_text()) for i, (a, x_) in enumerate(bb) for b_, y_ in bb[i + 1:] if _hit(x_, y_)]
    cl = [t.get_text() for t, x_ in bb if x_[:, 0].min() < fb.x0 - 1 or x_[:, 0].max() > fb.x1 + 1 or x_[:, 1].min() < fb.y0 - 1 or x_[:, 1].max() > fb.y1 + 1]
    return ov, cl, bb
def b_dot_overlaps():
    return [(p[3], q[3]) for i, p in enumerate(DOTS) for q in DOTS[:i] if math.hypot(p[0] - q[0], p[1] - q[1]) < p[2] + q[2] + .3]
def _near(poly, pts, margin):
    c_ = poly.mean(axis=0); d = np.full(len(pts), np.inf)
    for i in range(len(poly)):
        p, e_ = poly[i], poly[(i + 1) % len(poly)] - poly[i]; nrm = np.array([-e_[1], e_[0]]) / np.hypot(*e_)
        if (c_ - p) @ nrm < 0: nrm = -nrm
        d = np.minimum(d, (pts - p) @ nrm)
    return bool((d > -margin).any())
def b_text_hits(bb):
    """그림 전체 좌표 축(ax)의 글자가 b의 점이나 95% CI 호에 닿는지(글자 상자에서 2.5 pt 안)."""
    out = []; k = fig.dpi / 72
    for t, box in bb:
        if t.axes is not ax: continue
        for name, kind, r, lo, hi in ARCS:
            a = np.radians(np.linspace(lo, hi, 200))
            pts = np.c_[(QB.ox + r * np.cos(a)) * fig.dpi, (QB.oy + r * np.sin(a)) * fig.dpi]
            if _near(box, pts, 2.5 * k): out.append((t.get_text(), f"{name} {kind}"))
    return out
def a_text_in_b(bb):
    """a의 글자가 b의 그리는 영역에 들어가는지."""
    out = []
    for t, box in bb:
        if t.axes is not aA: continue
        x1, y1 = box[:, 0].max() / fig.dpi, box[:, 1].max() / fig.dpi
        if x1 > (QB.ox - QB.xm if y1 > QB.oy else QB.ox) - .08: out.append(t.get_text())
    return out
def cd_marker_overlaps():
    """c, d에서 한 줄의 두 표지가 서로 겹치는지(가운데 사이 거리 < 두 반지름의 합)."""
    out = []
    for name, p in PAN.items():
        a = CAX[name]
        for i, t in enumerate(ORD):
            (x0_, y0_), (x1_, y1_) = [a.transData.transform((float(df.c1[t]), i + sg * DY)) for (_, df, _, _), sg in zip(p["sets"], (-1, 1))]
            if math.hypot(x0_ - x1_, y0_ - y1_) * 72 / fig.dpi < FILLc["ms"] + .2: out.append((name, t))
    return out
ov, cl, bb = text_issues(fig)
print("text overlaps:", ov); print("clipped:", cl); print("b dot overlaps:", b_dot_overlaps()); print("b text on dots or arcs:", b_text_hits(bb))
print("a text inside b:", a_text_in_b(bb)); print("c, d marker overlaps:", cd_marker_overlaps())
assert not ov and not cl and not b_dot_overlaps() and not b_text_hits(bb) and not a_text_in_b(bb) and not cd_marker_overlaps()
fig.savefig(OUT / "fig3.png", dpi=300); fig.savefig(OUT / "fig3.pdf"); plt.close(fig)
NUM["figure_size_in"] = [W, H]
(OUT / "numbers_fig3.json").write_text(json.dumps(NUM, indent=1, default=float, ensure_ascii=False), encoding="utf-8")
print("b balance:", {t: (round(v["deg"], 1), round(v["se"], 1)) for t, v in NUM["b_balance"].items()})
print(f"margins: below axis {YM:.2f} in, left of axis {XM:.2f} in | lane pitch {LW_ * 72:.1f} pt")
print(f"saved {OUT / 'fig3.png'} ({W} × {H:.2f} in)")
