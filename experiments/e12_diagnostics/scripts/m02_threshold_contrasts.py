"""실험 3m M2: 역치별 층위 대비와 95 % 신뢰구간 (PROTOCOL_exp3m_2026-10-17, 기술; 새 실행 없음).
 모형: 전체 트랙 3층위 동시 모형, 역치 0.1·0.5·1 % (3g G3_*) + 5 % (3f T3all). 정렬 블록 delete 값(3j)으로 SE.
 재현 확인: 역치별 Q(공동 jackknife)가 3j j03_main_summary.csv의 (다) 값과 같아야 한다.
출력: results/m02_threshold_contrasts.csv, results/m02_tests.txt, figures/m02_threshold_contrasts.png
"""
import sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
sys.path.insert(0, str(Path(__file__).parent)); from m_common import *
EXP = Path(__file__).resolve().parents[1]; (EXP/"results").mkdir(exist_ok=True); (EXP/"figures").mkdir(exist_ok=True)
TH = [("G3_L0p1", "0.1%"), ("G3_L0p5", "0.5%"), ("G3_L1p0", "1%"), ("T3all", "5%")]
J3 = pd.read_csv(EXP.parent/"e09_joint_jackknife"/"results"/"j03_main_summary.csv")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
rows = []
for model, lab in TH:
    PR = {}
    for t in ALL:
        L = layers(t, [model]); c, D = prof(t, [model]); PR[t] = (c, D)
        top = max(LAY, key=lambda m: L[m][0])
        r = dict(trait=t, threshold=lab, top=top, **{f"tau_{m}": L[m][0] for m in LAY}, **{f"se_{m}": jse(L[m][1]) for m in LAY})
        for k, nm in enumerate(["c1", "c2"]):
            se = jse(D[:, k]); r.update({nm: c[k], f"{nm}_se": se, f"{nm}_lo": c[k] - 1.96 * se, f"{nm}_hi": c[k] + 1.96 * se, f"{nm}_p": 2 * stats.norm.sf(abs(c[k] / se))})
        rows.append(r)
    Q, df, p, pf = Qtest(PR, MAP)
    ref = J3[J3.test == f"[{lab.replace('%', ' %')}] Q"]
    if len(ref): assert abs(ref.stat_c.iloc[0] - Q) < 1e-6, (lab, Q, ref.stat_c.iloc[0])
    P(f"[{lab}] Q {Q:.1f} (자유도 {df}), P {p:.2g}, F 근사 {pf:.2g}" + (" — 3j 재현" if len(ref) else ""))
D = pd.DataFrame(rows); D.to_csv(EXP/"results"/"m02_threshold_contrasts.csv", index=False)
# 형질별 요약: 대비의 부호와 구간이 0을 넘는지, 1위 층위
P("\n형질별 (c1 = CAGE − DNase; 괄호 안 95 % CI; * = 구간이 0을 포함하지 않음) | 1위 층위")
for t in ALL:
    s = D[D.trait == t].set_index("threshold")
    f = lambda th, nm: f"{s.loc[th, nm]:+.2f} ({s.loc[th, nm + '_lo']:+.2f}, {s.loc[th, nm + '_hi']:+.2f}){'*' if s.loc[th, nm + '_lo'] > 0 or s.loc[th, nm + '_hi'] < 0 else ''}"
    P(f"  {KO[t]:6s} c1: " + " | ".join(f"{th} {f(th, 'c1')}" for _, th in TH) + "  || 1위: " + " ".join(s.loc[th, "top"] for _, th in TH))
P("\n형질별 (c2 = CAGE − H3K27ac)")
for t in ALL:
    s = D[D.trait == t].set_index("threshold")
    f = lambda th, nm: f"{s.loc[th, nm]:+.2f} ({s.loc[th, nm + '_lo']:+.2f}, {s.loc[th, nm + '_hi']:+.2f}){'*' if s.loc[th, nm + '_lo'] > 0 or s.loc[th, nm + '_hi'] < 0 else ''}"
    P(f"  {KO[t]:6s} c2: " + " | ".join(f"{th} {f(th, 'c2')}" for _, th in TH))
# 역치 사이 부호 일관성(기술): c1 부호가 4개 역치에서 같은 형질, 구간이 0을 넘는 역치 수
P("\n일관성 요약 (c1):")
for t in ALL:
    s = D[D.trait == t]; sg = np.sign(s.c1).unique(); sig = int(((s.c1_lo > 0) | (s.c1_hi < 0)).sum())
    P(f"  {KO[t]:6s} 부호 {'일정' if len(sg) == 1 else '바뀜'} ({'+' if s.c1.iloc[-1] > 0 else '−'} at 5%), 구간이 0 밖인 역치 {sig}/4")
(EXP/"results"/"m02_tests.txt").write_text("\n".join(log), encoding="utf-8")
# 그림: 형질별 행, 역치별 점 + 95 % CI, 두 패널
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
COL = ["#c6dbef", "#6baed6", "#2171b5", "#08306b"]   # 한 색상 연속(얕은 → 짙은 = 엄격 → 넓은 역치)
fig, axs = plt.subplots(1, 2, figsize=(7.2, 5.2), sharey=True)
order = ALL[::-1]
for ax, nm, ttl in zip(axs, ["c1", "c2"], ["CAGE − DNase", "CAGE − H3K27ac"]):
    for j, (_, th) in enumerate(TH):
        s = D[D.threshold == th].set_index("trait").loc[order]; y = np.arange(len(order)) + (j - 1.5) * 0.18
        ax.errorbar(s[nm], y, xerr=[s[nm] - s[nm + "_lo"], s[nm + "_hi"] - s[nm]], fmt="o", ms=3.2, lw=1, color=COL[j], label=f"top {th}")
    ax.axvline(0, color="#888", lw=.8); ax.set_title(ttl, fontsize=9); ax.set_xlabel("Difference in τ* (95% CI)", fontsize=8); ax.tick_params(labelsize=7)
    for b in [2.5]: ax.axhline(b, color="#ddd", lw=.6)
axs[0].set_yticks(range(len(order))); axs[0].set_yticklabels([EN[t] for t in order], fontsize=7)
axs[1].legend(fontsize=7, frameon=False, loc="lower right"); fig.tight_layout(); fig.savefig(EXP/"figures"/"m02_threshold_contrasts.png", dpi=200); P("그림 저장")
