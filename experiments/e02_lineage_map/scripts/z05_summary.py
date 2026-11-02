"""실험 3c z05: 형질 × 계통 τ* 지도 그림 + 결과 문서 초안 (판정은 z04 출력 그대로; 해석은 결과 확인 후 사람이 추가).
 그림: M1(농축)·M2(계통 특이) 두 패널, 같은 발산 색 척도(τ* 0 = 회색, 양 = 빨강, 음 = 파랑, ±4 에서 자름). q < 0.05 칸에만 값과 테두리.
출력: figures/exp3c_trait_lineage_map.png, results/z05_map_tau_star.csv (표 보기), results/RESULTS_EXP3c_<날짜>.md
"""
import datetime
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; (EXP/"figures").mkdir(exist_ok=True)
T = pd.read_csv(EXP/"results"/"z04_sldsc.csv")
LIN = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
LAB = {"CNS_tissue": "Brain tissue", "neural_cells": "Neural cells", "myeloid": "Myeloid", "lymphoid": "Lymphoid", "liver": "Liver", "GI": "GI", "cardiovascular": "Cardiovascular",
       "skeletal_muscle": "Skeletal muscle", "adipose": "Adipose", "stromal": "Stromal", "reproductive_endocrine": "Reprod./endocrine", "epithelial_barrier": "Epithelial/barrier"}
TR = ["AD", "EA", "SCZ", "BIP", "NEURO", "BMI", "MONO", "IBD", "RA", "HEIGHT", "LDL"]
TLAB = {"AD": "Alzheimer's disease", "EA": "Education years", "SCZ": "Schizophrenia", "BIP": "Bipolar disorder", "NEURO": "Neuroticism", "BMI": "BMI",
        "MONO": "Monocyte count", "IBD": "IBD", "RA": "Rheumatoid arthritis", "HEIGHT": "Height", "LDL": "LDL cholesterol"}
# 발산 색: 기준 팔레트의 blue 램프 ↔ red, 가운데 중립 회색
cmap = LinearSegmentedColormap.from_list("div", ["#104281", "#3987e5", "#9ec5f4", "#f0efec", "#f3aaa7", "#e34948", "#9e2324"])
norm = TwoSlopeNorm(vmin=-4, vcenter=0, vmax=4)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": "#c3c2b7", "text.color": "#0b0b0b", "axes.labelcolor": "#52514e", "xtick.color": "#52514e", "ytick.color": "#0b0b0b"})
fig, axes = plt.subplots(1, 2, figsize=(15, 6.2), sharey=True, facecolor="#fcfcfb")
tabs = {}
for ax, (m, title) in zip(axes, [("M1", "Enrichment map (baselineLD + lineage)"), ("M2", "Lineage-specific map (+ Atlas-all annotation)")]):
    d = T[T.model == m]; Z = d.pivot(index="trait", columns="lineage", values="tau_star").reindex(index=TR, columns=LIN)
    Q = d.pivot(index="trait", columns="lineage", values="q_bh").reindex(index=TR, columns=LIN); tabs[m] = Z
    ax.set_facecolor("#fcfcfb"); ax.imshow(Z.clip(-4, 4).to_numpy(), cmap=cmap, norm=norm, aspect="auto")
    for i in range(len(TR)):
        for j in range(len(LIN)):
            if Q.iat[i, j] < .05:
                v = Z.iat[i, j]; ax.add_patch(plt.Rectangle((j - .5 + .04, i - .5 + .04), .92, .92, fill=False, lw=1.2, ec="#0b0b0b"))
                ax.text(j, i, f"{v:+.1f}", ha="center", va="center", fontsize=7.5, color="#ffffff" if abs(v) > 2.6 else "#0b0b0b")
    ax.set_xticks(range(len(LIN))); ax.set_xticklabels([LAB[g] for g in LIN], rotation=45, ha="right")
    ax.set_yticks(range(len(TR))); ax.set_yticklabels([TLAB[t] for t in TR])
    for y in [0.5, 5.5]: ax.axhline(y, color="#fcfcfb", lw=9)   # AD | 뇌 형질 | 비뇌 형질 구분 (표면색 간격)
    ax.set_xticks(np.arange(-.5, len(LIN)), minor=True); ax.set_yticks(np.arange(-.5, len(TR)), minor=True)
    ax.grid(which="minor", color="#fcfcfb", lw=2); ax.tick_params(which="minor", length=0); ax.tick_params(length=0)
    ax.set_title(title, loc="left", fontsize=10.5, color="#0b0b0b", pad=8)
    for s in ax.spines.values(): s.set_visible(False)
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes, fraction=.025, pad=.015); cb.set_label("τ* (clipped at ±4)"); cb.outline.set_visible(False)
fig.text(.01, -.06, "Boxed cells with values: BH q < 0.05 across 132 trait × lineage tests within each map. Rows: AD | brain-related traits | other traits.", fontsize=8, color="#52514e")
fig.savefig(EXP/"figures"/"exp3c_trait_lineage_map.png", dpi=200, bbox_inches="tight", facecolor="#fcfcfb")
pd.concat({m: z for m, z in tabs.items()}, names=["model", "trait"]).round(3).to_csv(EXP/"results"/"z05_map_tau_star.csv")
today = datetime.date.today().isoformat()
txt = (EXP/"results"/"z04_sldsc.txt").read_text(encoding="utf-8")
(EXP/"results"/f"RESULTS_EXP3c_{today}.md").write_text(
    f"# 실험 3c 결과: 형질 × 계통 유전율 지도 ({today})\n\n자동 생성(z05). 판정은 PROTOCOL_exp3c의 α를 그대로 적용했다. 해석은 결과 확인 후 아래 '해석' 절에 추가한다.\n\n"
    f"## z04 출력\n```\n{txt}\n```\n\n![](../figures/exp3c_trait_lineage_map.png)\n\n- 표 보기: `results/z05_map_tau_star.csv` (τ\\*), `results/z04_sldsc.csv` (전체 수치)\n\n## 해석\n(결과 확인 후 작성)\n", encoding="utf-8")
print("z05 완료")
