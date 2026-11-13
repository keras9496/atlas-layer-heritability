"""실험 3 x05: x04 결과 요약 표 + 그림 (PROTOCOL_exp3 주 가설 판정은 x04 기준을 그대로 옮김; 해석은 따로 적는다).
 그림: (A) 세포형 τ* (단일·동시·AVI 조건부) (B) 맥락 × 모달리티 τ* 열지도 (C) 임계값 민감도 농축도 (D) 세포형 농축도
출력: results/RESULTS_EXP3_<날짜>.md, figures/exp3_sldsc.png
"""
import datetime
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; RES = EXP/"results"; FIG = EXP/"figures"; FIG.mkdir(exist_ok=True)
plt.rcParams["font.family"] = "Malgun Gothic"; plt.rcParams["axes.unicode_minus"] = False
T = pd.read_csv(RES/"x04_sldsc.csv"); TXT = (RES/"x04_sldsc.txt").read_text(encoding="utf-8")
CTX = ["myeloid", "neuron", "astrocyte", "liver", "nonliver"]; KO = {"myeloid": "골수계", "neuron": "뉴런", "astrocyte": "성상세포", "liver": "간", "nonliver": "기타 조직"}
MODS = ["DNASE", "H3K27ac", "H3K4me1", "H3K4me3", "CAGE", "RNA"]
g = lambda a, n: T[(T.analysis == a) & (T.annot == n)].iloc[0]
fig, ax = plt.subplots(2, 2, figsize=(13, 10))
# A
for j, (a, lab, col) in enumerate([("single|baselineLD", "단일(baselineLD)", "C0"), ("joint|ct5set", "5개 동시", "C1"), ("single|baselineLD+avi5", "+AVI 조건부", "C2")]):
    r = [g(a, f"ct5_{c}") for c in CTX]; x = np.arange(5) + (j - 1) * .25
    ax[0, 0].errorbar(x, [v.tau_star for v in r], yerr=[1.96 * v.tau_star_se for v in r], fmt="o", color=col, label=lab, capsize=3)
ax[0, 0].axhline(0, color="k", lw=.8); ax[0, 0].set_xticks(range(5)); ax[0, 0].set_xticklabels([KO[c] for c in CTX]); ax[0, 0].set_ylabel("τ* (95% CI)")
ax[0, 0].set_title("A. 세포형 상위 5 % 주석의 표준화 효과 τ*"); ax[0, 0].legend(fontsize=9)
# B
H = np.array([[g("single|baselineLD", f"mod_{c}_{m}").tau_star for m in MODS] for c in CTX]); Pm = np.array([[g("single|baselineLD", f"mod_{c}_{m}").tau_p1 for m in MODS] for c in CTX])
v = np.nanmax(np.abs(H)); im = ax[0, 1].imshow(H, cmap="RdBu_r", vmin=-v, vmax=v)
for i in range(5):
    for k in range(6): ax[0, 1].text(k, i, f"{H[i, k]:+.2f}" + ("**" if Pm[i, k] < .05 / 30 else "*" if Pm[i, k] < .05 else ""), ha="center", va="center", fontsize=8)
ax[0, 1].set_xticks(range(6)); ax[0, 1].set_xticklabels(MODS); ax[0, 1].set_yticks(range(5)); ax[0, 1].set_yticklabels([KO[c] for c in CTX])
ax[0, 1].set_title("B. 모달리티별 τ* (** Bonferroni 통과, * 단측 p<0.05)"); plt.colorbar(im, ax=ax[0, 1], shrink=.8)
# C
for j, t in enumerate(["ct1", "ct5", "ct10"]):
    r = [g("single|baselineLD", f"{t}_{c}") for c in CTX]; x = np.arange(5) + (j - 1) * .25
    ax[1, 0].errorbar(x, [v.enrichment for v in r], yerr=[1.96 * v.enrich_se for v in r], fmt="o", label={"ct1": "상위 1 %", "ct5": "상위 5 %", "ct10": "상위 10 %"}[t], capsize=3)
ax[1, 0].axhline(1, color="k", lw=.8); ax[1, 0].set_xticks(range(5)); ax[1, 0].set_xticklabels([KO[c] for c in CTX]); ax[1, 0].set_ylabel("농축도 (h² 비율 / SNP 비율)")
ax[1, 0].set_title("C. 임계값 민감도 (농축도)"); ax[1, 0].legend(fontsize=9)
# D
r = [g("single|baselineLD", n) for n in ["avi5"] + [f"ctc_{c}" for c in CTX]]
ax[1, 1].errorbar(range(6), [v.tau_star for v in r], yerr=[1.96 * v.tau_star_se for v in r], fmt="s", capsize=3, color="C3"); ax[1, 1].axhline(0, color="k", lw=.8)
ax[1, 1].set_xticks(range(6)); ax[1, 1].set_xticklabels(["AVI 상위 5 %"] + [f"{KO[c]}\n(연속)" for c in CTX]); ax[1, 1].set_title("D. AVI 대조와 연속 주석 τ*")
fig.suptitle("실험 3: Atlas 변이 효과 주석 × Bellenguez S-LDSC (baselineLD v2.2 조건부, APOE 제외)"); fig.tight_layout(); fig.savefig(FIG/"exp3_sldsc.png", dpi=150)
cols = ["analysis", "annot", "prop_snps", "prop_h2", "enrichment", "enrich_se", "enrich_p", "tau_star", "tau_star_se", "tau_p1"]
fmt = lambda d: d[cols].to_markdown(index=False, floatfmt=".3g")
md = [f"# 실험 3 결과: Atlas 주석 S-LDSC ({datetime.date.today()})", "", "자동 생성(x05). 판정은 PROTOCOL_exp3의 α를 그대로 적용. 해석은 아래 '해석' 절에 추가한다.", "",
      "## x04 출력 (주 가설 판정 포함)", "```", TXT, "```", "", "![](../figures/exp3_sldsc.png)", "",
      "## 세포형 주석 (단일 / 동시 / AVI 조건부)", fmt(T[T.annot.str.startswith("ct5_")]), "",
      "## 모달리티 주석", fmt(T[T.annot.str.startswith("mod_")].sort_values("tau_p1")), "",
      "## 민감도 (상위 1·10 %, 연속, AVI)", fmt(T[T.annot.str.match(r"ct1_|ct10_|ctc_|avi5")]), "", "## 해석", "(작성 예정)"]
(RES/f"RESULTS_EXP3_{datetime.date.today()}.md").write_text("\n".join(md), encoding="utf-8"); print("x05 완료")
