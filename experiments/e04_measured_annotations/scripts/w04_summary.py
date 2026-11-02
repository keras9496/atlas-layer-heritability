"""실험 3e w04: 판정·표·그림 (PROTOCOL_exp3e_2026-10-08).
 E-H1 (주): J1 의 Atlas L5_g_CAGE τ* > 0, 단측 p, 형질 11개 BH, q < 0.05 형질 수 → 읽는 규칙 (≥ 6 / 1–5 / 0).  E-H1-AD: AD 골수계 J1 단측 α 0.05.
 기술: 유지율 τ*(J1)/τ*(D1), S_m 같은 모달리티 비교와 차이, 실측 감쇠 M1_m → S_m, 민감도 J1k(분모 3d-A D1k)·J1p.
입력: results/w03_sldsc.csv, w03_diffs.csv, 실험3d results/y04_B.csv (D1), y06_eqtl_compare.csv (D1k)
출력: results/w04_EH1.csv, w04_Sm.csv, w04_summary.txt, figures/exp3e_conditional.png
"""
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Malgun Gothic"; plt.rcParams["axes.unicode_minus"] = False
EXP = Path(__file__).resolve().parents[1]; R3D = EXP.parent/"e03_single_layer_map"/"results"; (EXP/"figures").mkdir(exist_ok=True)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); q = p[o] * len(p) / np.arange(1, len(p) + 1); q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty_like(q); out[o] = np.minimum(q, 1); return out
TR = ["AD", "MONO", "EA", "SCZ", "BIP", "NEURO", "BMI", "IBD", "RA", "LDL", "HEIGHT"]; KO = {"AD": "AD", "MONO": "단핵구 수", "EA": "교육 연수", "SCZ": "조현병", "BIP": "양극성장애",
      "NEURO": "신경증", "BMI": "BMI", "IBD": "IBD", "RA": "RA", "LDL": "LDL", "HEIGHT": "키"}
T = pd.read_csv(EXP/"results"/"w03_sldsc.csv"); Dd = pd.read_csv(EXP/"results"/"w03_diffs.csv")
D1 = pd.read_csv(R3D/"y04_B.csv"); D1k = pd.read_csv(R3D/"y06_eqtl_compare.csv")
def get(t, model, annot):
    r = T[(T.trait == t) & (T.model == model) & (T.annot == annot)]; assert len(r) == 1, (t, model, annot); return r.iloc[0]
# ---------- E-H1 ----------
rows = []
for t in TR:
    g = T[T.trait == t].lineage.iloc[0]; a = f"L5_{g}_CAGE"; j = get(t, "J1", a)
    d1 = D1[(D1.trait == t) & (D1.model == "D1") & (D1.annot == a)].iloc[0]
    dk = D1k[(D1k.trait == t) & (D1k.model == "D1k") & (D1k.atlas == f"L5k_{g}_CAGE")].iloc[0]; jk = get(t, "J1k", f"L5k_{g}_CAGE"); jp = get(t, "J1p", a)
    rec = dict(trait=t, lineage=g, D1_tau=d1.tau_star, D1_se=d1.tau_star_se, J1_tau=j.tau_star, J1_se=j.tau_star_se, J1_p1=j.tau_p1, retention=j.tau_star / d1.tau_star,
               J1k_tau=jk.tau_star, J1k_se=jk.tau_star_se, J1k_p1=jk.tau_p1, D1k_tau=dk.atlas_tau, retention_k=jk.tau_star / dk.atlas_tau, J1p_tau=jp.tau_star, J1p_se=jp.tau_star_se, J1p_p1=jp.tau_p1)
    for m in ["CAGE", "DNASE", "H3K27ac"]:
        b = get(t, "BM", f"M_{g}_{m}"); x = get(t, "J1", f"M_{g}_{m}"); rec[f"BM_{m}_tau"] = b.tau_star; rec[f"BM_{m}_se"] = b.tau_star_se; rec[f"J1_M_{m}_tau"] = x.tau_star
    rows.append(rec)
E = pd.DataFrame(rows); E["J1_q"] = bh(E.J1_p1); E["J1k_q"] = bh(E.J1k_p1); E["J1p_q"] = bh(E.J1p_p1)
E.to_csv(EXP/"results"/"w04_EH1.csv", index=False)
n = int((E.J1_q < .05).sum()); passed = E[E.J1_q < .05].trait.tolist()
P("=== E-H1 (주): 실측 계통 CAGE·DNase·H3K27ac를 넣은 모형(J1)의 Atlas CAGE τ* ===")
P("형질        계통        D1 τ* (단독)     J1 τ* (조건부)    단측 p     q      유지율 | 같은 모형 실측 τ*: CAGE / DNase / H3K27ac")
for r in E.itertuples():
    P(f"{KO[r.trait]:9s}  {r.lineage:10s}  {r.D1_tau:+5.2f} ± {r.D1_se:.2f}   {r.J1_tau:+5.2f} ± {r.J1_se:.2f}   {r.J1_p1:.2g}   {r.J1_q:.2g}   {r.retention * 100:4.0f} % | "
      f"{r.J1_M_CAGE_tau:+.2f} / {r.J1_M_DNASE_tau:+.2f} / {r.J1_M_H3K27ac_tau:+.2f}")
rule = ("6/11 이상 → 예측 CAGE 효과는 형질에 맞는 조직의 실측 CAGE·DNase·H3K27ac 주석을 넘는 유전율 정보를 준다" if n >= 6 else
        "1–5/11 → 일부 형질에서만 실측을 넘는다" if n >= 1 else "0/11 → Atlas CAGE 농축은 대부분 조직 활성 위치로 설명된다")
P(f"\nE-H1: q < 0.05 형질 {n}/11 ({', '.join(KO[t] for t in passed) or '없음'}) → 읽는 규칙: {rule}")
ad = E[E.trait == "AD"].iloc[0]; P(f"E-H1-AD: AD 골수계 J1 τ* {ad.J1_tau:+.2f} ± {ad.J1_se:.2f}, 단측 p {ad.J1_p1:.3g} → {'통과' if ad.J1_p1 < .05 else '불통과'} (α 0.05)")
P(f"유지율 τ*(J1)/τ*(D1): 중앙값 {E.retention.median() * 100:.0f} % (범위 {E.retention.min() * 100:.0f}–{E.retention.max() * 100:.0f} %) → "
  f"{'대부분 실측과 독립' if E.retention.median() >= .5 else '상당 부분 겹침'} (기준 50 %)")
P(f"\n민감도: J1k(트랙 맞춤) q < 0.05 {int((E.J1k_q < .05).sum())}/11, 유지율 중앙값 {E.retention_k.median() * 100:.0f} % (분모 3d-A D1k) | "
  f"J1p(창 없는 피크) q < 0.05 {int((E.J1p_q < .05).sum())}/11, τ* 중앙값 {E.J1p_tau.median():+.2f}")
# ---------- S_m ----------
srows = []
for t in TR:
    g = E[E.trait == t].lineage.iloc[0]
    for m in ["CAGE", "DNASE", "H3K27ac"]:
        a = get(t, f"S_{m}", f"L5_{g}_{m}"); x = get(t, f"S_{m}", f"M_{g}_{m}"); x1 = get(t, f"M1_{m}", f"M_{g}_{m}"); d = Dd[(Dd.trait == t) & (Dd.model == f"S_{m}")].iloc[0]
        d1 = D1[(D1.trait == t) & (D1.model == "D1") & (D1.annot == f"L5_{g}_{m}")].iloc[0]
        srows.append(dict(trait=t, lineage=g, mod=m, atlas_D1=d1.tau_star, atlas_S=a.tau_star, atlas_S_se=a.tau_star_se, atlas_S_p1=a.tau_p1, meas_M1=x1.tau_star, meas_M1_se=x1.tau_star_se,
                          meas_S=x.tau_star, meas_S_se=x.tau_star_se, diff=d["diff"], diff_se=d.se, diff_p2=d.p2, atlas_ret=a.tau_star / d1.tau_star, meas_att=(x.tau_star - x1.tau_star) / abs(x1.tau_star)))
S = pd.DataFrame(srows); S["atlas_S_q"] = bh(S.atlas_S_p1); S.to_csv(EXP/"results"/"w04_Sm.csv", index=False)
P("\n=== 같은 모달리티 비교 (S_m = baselineLD + 실측 M_g_m + Atlas L5_g_m) ===")
P("형질        모달리티  Atlas τ* (D1→S)        실측 τ* (단독→S)        Atlas − 실측 (p)    실측 감쇠")
for r in S.itertuples():
    P(f"{KO[r.trait]:9s}  {r.mod:8s}  {r.atlas_D1:+5.2f} → {r.atlas_S:+5.2f} ± {r.atlas_S_se:.2f}   {r.meas_M1:+5.2f} → {r.meas_S:+5.2f} ± {r.meas_S_se:.2f}   "
      f"{r.diff:+5.2f} ({r.diff_p2:.2g})   {r.meas_att * 100:+4.0f} %")
for m in ["CAGE", "DNASE", "H3K27ac"]:
    s = S[S["mod"] == m]
    P(f"  {m:8s}: Atlas τ* > 0 (S, BH q < 0.05 / 33칸) {int((s.atlas_S_q < .05).sum())}/11 | Atlas > 실측 점추정 {int((s['diff'] > 0).sum())}/11, 차이 p < 0.05 "
      f"(Atlas 큼 {int(((s['diff'] > 0) & (s.diff_p2 < .05)).sum())}, 실측 큼 {int(((s['diff'] < 0) & (s.diff_p2 < .05)).sum())}) | 실측 감쇠 중앙값 {s.meas_att.median() * 100:+.0f} %")
# ---------- 그림 ----------
fig, ax = plt.subplots(1, 2, figsize=(12, 5.2), gridspec_kw=dict(width_ratios=[1.1, 1]))
y = np.arange(len(E))[::-1]
for k, (col, se, lab, c, dy) in enumerate([("D1_tau", "D1_se", "Atlas CAGE 단독 (3d D1)", "#9aa3ad", .17), ("J1_tau", "J1_se", "Atlas CAGE | 실측 CAGE·DNase·H3K27ac (J1)", "#1f4e8c", -.17)]):
    ax[0].errorbar(E[col], y + dy, xerr=1.96 * E[se], fmt="o", color=c, ms=6, capsize=0, lw=1.6, label=lab)
ax[0].axvline(0, color="#555", lw=.8); ax[0].set_yticks(y); ax[0].set_yticklabels([f"{KO[t]} · {g}" for t, g in zip(E.trait, E.lineage)])
for yy, q in zip(y, E.J1_q): ax[0].text(ax[0].get_xlim()[1], yy - .17, "  *" if q < .05 else "", va="center", fontsize=12, color="#1f4e8c")
ax[0].set_xlabel("τ* (± 95 % CI)"); ax[0].set_title("A. 실측 계통 주석 조건부 후 Atlas CAGE", loc="left", fontsize=11); ax[0].legend(loc="lower right", fontsize=8, frameon=False)
mods = ["CAGE", "DNASE", "H3K27ac"]; xs = np.arange(len(mods))
for i, m in enumerate(mods):
    s = S[S["mod"] == m]; jit = np.random.default_rng(i).uniform(-.05, .05, len(s)); xa, xm = i - .15 + jit, i + .15 + jit
    for a, b, u, v in zip(s.atlas_S, s.meas_S, xa, xm): ax[1].plot([u, v], [a, b], color="#d0d0d0", lw=.8, zorder=1)
    ax[1].scatter(xa, s.atlas_S, color="#1f4e8c", s=28, label="Atlas 예측 (S_m)" if i == 0 else None, zorder=3)
    ax[1].scatter(xm, s.meas_S, color="#d9822b", s=28, label="실측 (S_m)" if i == 0 else None, zorder=3)
ax[1].axhline(0, color="#555", lw=.8); ax[1].set_xticks(xs); ax[1].set_xticklabels(mods); ax[1].set_ylabel("τ* (같은 모형)")
ax[1].set_title("B. 같은 모달리티: 예측 대 실측 (형질 11개)", loc="left", fontsize=11); ax[1].legend(fontsize=8, frameon=False)
for a in ax: a.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig(EXP/"figures"/"exp3e_conditional.png", dpi=200)
(EXP/"results"/"w04_summary.txt").write_text("\n".join(log), encoding="utf-8")
