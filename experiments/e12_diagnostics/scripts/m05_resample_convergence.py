"""실험 3m M5: 재추출 수렴 (PROTOCOL_exp3m_2026-10-17, 기술; 새 실행 없음).
 3f 재추출 20회(T3_01–T3_20, 정렬 블록 npz 점 추정치 = 원래 점 추정치)에서:
  - 재추출 수 n = 1–20에 따른 τ* 누적 평균과 1위 층위 빈도
  - 재추출 단위 부트스트랩(2,000회, 시드 20261017)으로 20회 평균의 SE(트랙 선택에 따른 흔들림)와 평균 기준 1위 층위가 유지되는 비율
  - 비교: 표본 추출에 따른 불확실성(블록 jackknife SE, 20회 평균 delete 값)
 재추출 확대(20회 초과)는 하지 않는다.
출력: results/m05_convergence.csv, m05_tests.txt, figures/m05_convergence.png
"""
import sys
import numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent)); from m_common import *
EXP = Path(__file__).resolve().parents[1]
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
rng = np.random.default_rng(20261017); NB = 2000; rows = []; CUM = {}
for t in ALL:
    T = np.array([[load(t, s)[m][0] for m in LAY] for s in SEEDS])                 # 20 × 3
    jk = layers(t, SEEDS); jse_m = np.array([jse(jk[m][1]) for m in LAY])
    cum = np.cumsum(T, 0) / np.arange(1, 21)[:, None]; CUM[t] = cum
    top20 = int(np.argmax(cum[-1])); topk = T.argmax(1)
    stable_from = next(n for n in range(1, 21) if all(np.argmax(cum[k]) == top20 for k in range(n - 1, 20)))
    B = T[rng.integers(0, 20, size=(NB, 20))].mean(1)                               # NB × 3
    bse = B.std(0, ddof=1); keep = float((B.argmax(1) == top20).mean())
    d10 = cum[-1] - cum[9]
    rows.append(dict(trait=t, top=LAY[top20], n_first=int((topk == top20).sum()), stable_from_n=stable_from, boot_keep_top=keep,
                     **{f"mean_{m}": cum[-1][i] for i, m in enumerate(LAY)}, **{f"bootse_{m}": bse[i] for i, m in enumerate(LAY)},
                     **{f"jkse_{m}": jse_m[i] for i, m in enumerate(LAY)}, **{f"change10to20_{m}": d10[i] for i, m in enumerate(LAY)}))
    P(f"{KO[t]:6s} 1위 {LAY[top20]:7s} (20회 중 {int((topk == top20).sum())}회) | 평균 기준 1위가 n ≥ {stable_from}부터 고정 | 부트스트랩 1위 유지 {keep:.1%} | "
      f"재추출 SE C/D/H {bse[0]:.2f}/{bse[1]:.2f}/{bse[2]:.2f} vs jackknife SE {jse_m[0]:.2f}/{jse_m[1]:.2f}/{jse_m[2]:.2f} | 10→20회 변화 {d10[0]:+.2f}/{d10[1]:+.2f}/{d10[2]:+.2f}")
D = pd.DataFrame(rows); D.to_csv(EXP/"results"/"m05_convergence.csv", index=False)
r = (D[[f"bootse_{m}" for m in LAY]].to_numpy() / D[[f"jkse_{m}" for m in LAY]].to_numpy())
P(f"\n재추출 SE / jackknife SE: 중앙값 {np.median(r):.2f}, 최대 {r.max():.2f} | 10→20회 평균 변화 |Δ| 최대 {D[[f'change10to20_{m}' for m in LAY]].abs().to_numpy().max():.2f}")
(EXP/"results"/"m05_tests.txt").write_text("\n".join(log), encoding="utf-8")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
COL = {"CAGE": "#2a6fb0", "DNASE": "#d08a2c", "H3K27ac": "#5b9a5b"}
fig, axs = plt.subplots(4, 4, figsize=(7.2, 6.4), sharex=True); axs = axs.ravel()
for ax, t in zip(axs, ALL):
    for i, m in enumerate(LAY): ax.plot(range(1, 21), CUM[t][:, i], color=COL[m], lw=1.4, label={"DNASE": "DNase"}.get(m, m))
    ax.set_title(EN[t], fontsize=7); ax.tick_params(labelsize=6); ax.axhline(0, color="#bbb", lw=.6)
for ax in axs[len(ALL):]: ax.axis("off")
axs[len(ALL)].legend(*axs[0].get_legend_handles_labels(), fontsize=7, frameon=False, loc="center")
fig.supxlabel("Number of resamples averaged", fontsize=8); fig.supylabel("Cumulative mean τ*", fontsize=8); fig.tight_layout(); fig.savefig(EXP/"figures"/"m05_convergence.png", dpi=200)
P("그림 저장")
