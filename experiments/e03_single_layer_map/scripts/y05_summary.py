"""실험 3d y05: 그림 + 결과 문서 초안 (있는 단계만 그린다; 판정은 y04 출력 그대로, 해석은 결과 확인 후 따로 적는다).
 그림 1 (단계 B 후): D1 모달리티 지도 — 형질 11행 × (계통 12 × 모달리티 6)열, 발산 색(τ* 0 = 회색), q < 0.05 칸에 점
 그림 2 (단계 A 후): 19쌍 D2 모달리티 프로필 — 막대 = D2 τ* ± SE, 속 빈 마름모 = D2k(트랙 맞춤). 모달리티 이름을 x축에 직접 표시
 그림 3 (단계 A 후): D0 형질 × 모달리티 지도
출력: figures/exp3d_*.png, results/y05_*.csv (표 보기), results/RESULTS_EXP3d_<날짜>.md
"""
import datetime
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; R = EXP/"results"; F = EXP/"figures"; F.mkdir(exist_ok=True)
LIN = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
LAB = {"CNS_tissue": "Brain tissue", "neural_cells": "Neural cells", "myeloid": "Myeloid", "lymphoid": "Lymphoid", "liver": "Liver", "GI": "GI", "cardiovascular": "Cardiovasc.",
       "skeletal_muscle": "Skel. muscle", "adipose": "Adipose", "stromal": "Stromal", "reproductive_endocrine": "Reprod./endo.", "epithelial_barrier": "Epithelial"}
MODS = ["DNASE", "H3K27ac", "H3K4me1", "H3K4me3", "CAGE", "RNA"]; MLAB = {"DNASE": "DNase", "H3K27ac": "H3K27ac", "H3K4me1": "H3K4me1", "H3K4me3": "H3K4me3", "CAGE": "CAGE", "RNA": "RNA"}
MCOL = dict(zip(MODS, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]))   # 검증 통과(대비 경고 → 이름 직접 표시 + CSV)
TR = ["AD", "EA", "SCZ", "BIP", "NEURO", "BMI", "MONO", "IBD", "RA", "HEIGHT", "LDL"]
TLAB = {"AD": "Alzheimer's disease", "EA": "Education years", "SCZ": "Schizophrenia", "BIP": "Bipolar disorder", "NEURO": "Neuroticism", "BMI": "BMI",
        "MONO": "Monocyte count", "IBD": "IBD", "RA": "Rheumatoid arthritis", "HEIGHT": "Height", "LDL": "LDL cholesterol"}
cmap = LinearSegmentedColormap.from_list("div", ["#104281", "#3987e5", "#9ec5f4", "#f0efec", "#f3aaa7", "#e34948", "#9e2324"]); norm = TwoSlopeNorm(vmin=-4, vcenter=0, vmax=4)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": "#c3c2b7", "text.color": "#0b0b0b", "axes.labelcolor": "#52514e", "xtick.color": "#52514e", "ytick.color": "#0b0b0b"})
SURF = "#fcfcfb"; made = []
if (R/"y04_A.csv").exists():
    A = pd.read_csv(R/"y04_A.csv")
    # 그림 3: D0
    Z = A[A.model == "D0"].pivot(index="trait", columns="mod", values="tau_star").reindex(index=TR, columns=MODS); Z.round(3).to_csv(R/"y05_D0_tau_star.csv")
    fig, ax = plt.subplots(figsize=(6.2, 5.6), facecolor=SURF); ax.set_facecolor(SURF)
    ax.imshow(Z.clip(-4, 4).to_numpy(), cmap=cmap, norm=norm, aspect="auto")
    for i in range(len(TR)):
        for j in range(len(MODS)): ax.text(j, i, f"{Z.iat[i, j]:+.1f}", ha="center", va="center", fontsize=7.5, color="#ffffff" if abs(Z.iat[i, j]) > 2.6 else "#0b0b0b")
    ax.set_xticks(range(6)); ax.set_xticklabels([MLAB[m] for m in MODS], rotation=30, ha="right"); ax.set_yticks(range(len(TR))); ax.set_yticklabels([TLAB[t] for t in TR])
    ax.set_xticks(np.arange(-.5, 6), minor=True); ax.set_yticks(np.arange(-.5, len(TR)), minor=True); ax.grid(which="minor", color=SURF, lw=2); ax.tick_params(which="both", length=0)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_title("Modality enrichment, all lineages pooled (τ*)", loc="left", fontsize=10.5, pad=8)
    fig.savefig(F/"exp3d_D0_trait_modality.png", dpi=200, bbox_inches="tight", facecolor=SURF); plt.close(fig); made.append("exp3d_D0_trait_modality.png")
    # 그림 2: 19쌍 프로필
    pairs = sorted(A[A.model == "D2k"][["trait", "unit"]].drop_duplicates().values.tolist(), key=lambda p: (TR.index(p[0]), LIN.index(p[1])))   # 지도와 같은 형질 순서
    nc = 5; nr = int(np.ceil(len(pairs) / nc)); fig, axes = plt.subplots(nr, nc, figsize=(15, 2.7 * nr), facecolor=SURF, sharey=True); axes = axes.ravel()
    yl = A[A.model.isin(["D2", "D2k"])]; ymax = np.nanmax(np.abs(yl.tau_star + np.sign(yl.tau_star) * yl.tau_star_se)); ymax = min(max(ymax, 2), 8)
    rows = []
    for ax, (t, g) in zip(axes, pairs):
        a = A[(A.model == "D2") & (A.trait == t) & (A.unit == g)].set_index("mod").reindex(MODS); k = A[(A.model == "D2k") & (A.trait == t) & (A.unit == g)].set_index("mod").reindex(MODS)
        ax.set_facecolor(SURF); x = np.arange(6)
        ax.bar(x, a.tau_star, width=.7, color=[MCOL[m] for m in MODS], edgecolor=SURF, linewidth=2)
        ax.errorbar(x, a.tau_star, yerr=a.tau_star_se, fmt="none", ecolor="#52514e", elinewidth=1, capsize=0)
        ax.scatter(x, k.tau_star, marker="D", s=26, facecolor=SURF, edgecolor="#0b0b0b", linewidth=1, zorder=3)
        ax.axhline(0, color="#c3c2b7", lw=1); ax.set_ylim(-ymax, ymax)
        ax.set_xticks(x); ax.set_xticklabels([MLAB[m] for m in MODS], rotation=55, ha="right", fontsize=7.5); ax.tick_params(length=0)
        ax.set_title(f"{TLAB[t]} · {LAB[g]}", loc="left", fontsize=8.5)
        for s in ["top", "right"]: ax.spines[s].set_visible(False)
        for m in MODS: rows.append(dict(trait=t, lineage=g, mod=m, D2=a.tau_star[m], D2_se=a.tau_star_se[m], D2_p=a.tau_p1[m], D2k=k.tau_star[m]))
    for ax in axes[len(pairs):]: ax.axis("off")
    axes[0].set_ylabel("τ* (joint, same lineage)")
    fig.text(.01, -.01, "Bars: τ* ± SE with the six modalities of one lineage fitted jointly (D2). Open diamonds: same model with ≤3 tracks per modality (D2k, track-count matched).", fontsize=8, color="#52514e")
    fig.tight_layout(); fig.savefig(F/"exp3d_D2_modality_profiles.png", dpi=200, bbox_inches="tight", facecolor=SURF); plt.close(fig); made.append("exp3d_D2_modality_profiles.png")
    pd.DataFrame(rows).round(4).to_csv(R/"y05_D2_profiles.csv", index=False)
if (R/"y04_B.csv").exists():
    Bt = pd.read_csv(R/"y04_B.csv"); cols = [(g, m) for g in LIN for m in MODS]
    Z = Bt.pivot_table(index="trait", columns=["unit", "mod"], values="tau_star").reindex(index=TR, columns=pd.MultiIndex.from_tuples(cols))
    Q = Bt.pivot_table(index="trait", columns=["unit", "mod"], values="q_bh").reindex(index=TR, columns=pd.MultiIndex.from_tuples(cols)); Z.round(3).to_csv(R/"y05_D1_map_tau_star.csv")
    fig, ax = plt.subplots(figsize=(18, 6.4), facecolor=SURF); ax.set_facecolor(SURF)
    ax.imshow(Z.clip(-4, 4).to_numpy(), cmap=cmap, norm=norm, aspect="auto")
    yy, xx = np.where(Q.to_numpy() < .05); ax.scatter(xx, yy, s=10, color="#0b0b0b", zorder=3)
    ax.set_yticks(range(len(TR))); ax.set_yticklabels([TLAB[t] for t in TR])
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([MLAB[m] for _, m in cols], rotation=90, fontsize=6.5)
    for i, g in enumerate(LIN):
        ax.text(i * 6 + 2.5, -0.9, LAB[g], ha="center", va="bottom", fontsize=8.5, color="#0b0b0b")
        if i: ax.axvline(i * 6 - .5, color=SURF, lw=5)
    for y in [0.5, 5.5]: ax.axhline(y, color=SURF, lw=7)
    ax.set_xticks(np.arange(-.5, len(cols)), minor=True); ax.set_yticks(np.arange(-.5, len(TR)), minor=True); ax.grid(which="minor", color=SURF, lw=1); ax.tick_params(which="both", length=0)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_title("Trait × lineage × modality heritability map (τ*, single annotation per cell)", loc="left", fontsize=11, pad=26)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=.015, pad=.01); cb.set_label("τ* (clipped at ±4)"); cb.outline.set_visible(False)
    fig.text(.01, -.03, "Dots: BH q < 0.05 across 792 cells. Rows: AD | brain-related traits | other traits.", fontsize=8, color="#52514e")
    fig.savefig(F/"exp3d_D1_modality_map.png", dpi=200, bbox_inches="tight", facecolor=SURF); plt.close(fig); made.append("exp3d_D1_modality_map.png")
today = datetime.date.today().isoformat()
txt = "\n\n".join(f"### {p.name}\n```\n{p.read_text(encoding='utf-8')}\n```" for p in [R/"y04_A.txt", R/"y04_B.txt"] if p.exists())
(R/f"RESULTS_EXP3d_{today}.md").write_text(
    f"# 실험 3d 결과: 형질 × 계통 × 모달리티 유전율 지도 ({today})\n\n자동 생성(y05). 판정은 PROTOCOL_exp3d의 α를 그대로 적용했다. 해석은 결과 확인 후 아래 '해석' 절에 추가한다.\n\n"
    f"## y04 출력\n{txt}\n\n## 그림\n" + "\n".join(f"![](../figures/{m})" for m in made) +
    "\n\n- 표 보기: `results/y05_*.csv`, 전체 수치 `results/y04_*.csv`\n\n## 해석\n(결과 확인 후 작성)\n", encoding="utf-8")
print("y05 완료:", made)
