"""실험 3d y04: S-LDSC (PROTOCOL_exp3d_2026-10-02). 순차 실행. 호출은 실험 3c z04와 같다.
 단계 A: D0 (L5_all_m 단일, 66) → D2 (계통 g 의 6 모달리티 동시, 132) → D2k (트랙 맞춤 동시, 19쌍)
 단계 B: D1 (L5_{g,m} 단일, 792) — 모달리티 지도
 판정(단계 A): D-H1 AD 골수계 CAGE (D2, α 0.05), D-H2 같은 것 D2k, D-H3 19쌍 1위 모달리티 집중도(다항 균등 귀무, 최다 범주 정확 검정)
 재실행 안전: .results 가 있으면 건너뛴다.  사용: python y04_sldsc.py A | B
출력: exp3d/sldsc/<형질>__<모형>__<주석>.*, results/y04_A.csv, y04_A_diffs.csv, y04_A.txt (단계 A), results/y04_B.csv, y04_B.txt (단계 B)
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess, os, re, sys, time, math
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); Y = Path(f"{DATA}/exp3d"); AN = Y/"annot"; SL = Y/"sldsc"; SL.mkdir(exist_ok=True)
BV = Path(f"{DATA}/tools/ldsc_bv"); B = X/"GRCh38"; SS = Path(f"{DATA}/exp3b/sumstats")
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LIN = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
MODS = ["DNASE", "H3K27ac", "H3K4me1", "H3K4me3", "CAGE", "RNA"]
TRAITS = {"AD": X/"bellenguez_noAPOE.sumstats.gz", "EA": SS/"PASS_Years_of_Education2.sumstats", "SCZ": SS/"PASS_Schizophrenia.sumstats", "BIP": SS/"PASS_Bipolar_Disorder.sumstats",
          "NEURO": SS/"PASS_Neuroticism.sumstats", "BMI": SS/"PASS_BMI1.sumstats", "MONO": Path(f"{DATA}/exp3b/monocyte_count.sumstats.gz"), "IBD": SS/"PASS_IBD.sumstats",
          "RA": SS/"PASS_Rheumatoid_Arthritis.sumstats", "HEIGHT": SS/"PASS_Height1.sumstats", "LDL": SS/"PASS_LDL.sumstats"}
PAIRS = [("AD", "myeloid"), ("AD", "lymphoid"), ("BIP", "myeloid"), ("EA", "CNS_tissue"), ("HEIGHT", "adipose"), ("HEIGHT", "skeletal_muscle"), ("HEIGHT", "stromal"),
         ("HEIGHT", "cardiovascular"), ("IBD", "lymphoid"), ("IBD", "myeloid"), ("IBD", "liver"), ("LDL", "liver"), ("LDL", "GI"), ("LDL", "adipose"), ("MONO", "myeloid"),
         ("MONO", "lymphoid"), ("NEURO", "CNS_tissue"), ("RA", "lymphoid"), ("SCZ", "CNS_tissue")]   # 실험 3c M1 q < 0.05 (프로토콜 고정)
SD = pd.read_csv(Y/"annot_sd_common.csv", index_col=0).sd_common
M_TOT = sum(float(open(B/"baselineLD_v2.2"/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
def run(t, model, a):
    out = SL/f"{t}__{model}__{a}"
    if not Path(f"{out}.results").exists():
        ref = ",".join([str(B/"baselineLD_v2.2"/"baselineLD."), str(AN/a/f"{a}.")])
        cmd = ["python", "-W", "ignore", "-c", "import sys; sys.argv=['ldsc']+sys.argv[1:]; from ldscore.ldsc import main; main()",
               "--h2", str(TRAITS[t]), "--ref-ld-chr", ref, "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."), "--overlap-annot",
               "--frqfile-chr", str(B/"plink_files"/"1000G.EUR.hg38."), "--print-coefficients", "--print-delete-vals", "--out", str(out)]
        for attempt in range(3):                       # LDSC occasionally fails to read the LD score files ("Error parsing reference panel LD Score")
            subprocess.run(cmd, cwd=BV, env=env, capture_output=True, text=True)   # on an otherwise intact annotation; a plain retry succeeds, so retry up to 3 times
            if Path(f"{out}.results").exists(): break
            time.sleep(30)
        if not Path(f"{out}.results").exists(): raise SystemExit(f"{t} {model} {a} 실패: {open(str(out) + '.log', encoding='utf-8').read()[-1200:]}")
    R = pd.read_csv(f"{out}.results", sep="\t"); lg = open(f"{out}.log", encoding="utf-8").read()
    h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", lg).group(1)); return R, h2, out
def rec(t, model, unit, annot, R, h2, k):
    r = R.iloc[k]; f = SD[annot] * M_TOT / h2
    return dict(trait=t, model=model, unit=unit, annot=annot, h2_obs=h2, prop_h2=r["Prop._h2"], enrichment=r.Enrichment, enrich_p=r.Enrichment_p,
                tau_star=r.Coefficient * f, tau_star_se=r.Coefficient_std_error * f, tau_p1=stats.norm.sf(r.Coefficient / r.Coefficient_std_error))
def bh(p):
    p = np.asarray(p); o = np.argsort(p); q = p[o] * len(p) / np.arange(1, len(p) + 1); q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty_like(q); out[o] = np.minimum(q, 1); return out
def joint(t, model, g, prefix):
    R, h2, out = run(t, model, f"{prefix}_{g}"); k0 = len(R) - 6; rows = []; diffs = []
    names = [f"{'L5k' if prefix == 'Jk' else 'L5'}_{g}_{m}" for m in MODS]
    for j, m in enumerate(MODS): rows.append(rec(t, model, g, names[j], R, h2, k0 + j) | {"mod": m})
    D = pd.read_csv(f"{out}.part_delete", sep=r"\s+", header=None).to_numpy(); nb = D.shape[0]
    for i in range(6):
        for j in range(i + 1, 6):
            d = (D[:, k0 + i] * SD[names[i]] - D[:, k0 + j] * SD[names[j]]) * M_TOT / h2; est = rows[i]["tau_star"] - rows[j]["tau_star"]
            se = np.sqrt((nb - 1) / nb * ((d - d.mean()) ** 2).sum()); diffs.append(dict(trait=t, model=model, lineage=g, a=MODS[i], b=MODS[j], diff=est, se=se, p2=2 * stats.norm.sf(abs(est / se))))
    return rows, diffs
t0 = time.time()
if STAGE == "A":
    rows, diffs = [], []
    jobs = [("D0", t, m) for t in TRAITS for m in MODS] + [("D2", t, g) for t in TRAITS for g in LIN] + [("D2k", t, g) for t, g in PAIRS]
    P(f"단계 A: S-LDSC {len(jobs)}회, 순차")
    for i, (model, t, x) in enumerate(jobs):
        if model == "D0": R, h2, _ = run(t, "D0", f"L5_all_{x}"); rows.append(rec(t, "D0", "all", f"L5_all_{x}", R, h2, len(R) - 1) | {"mod": x})
        else: r_, d_ = joint(t, model, x, "J" if model == "D2" else "Jk"); rows += r_; diffs += d_
        if (i + 1) % 11 == 0: P(f"  {i + 1}/{len(jobs)} ({(time.time() - t0) / 60:.0f}분)")
    T = pd.DataFrame(rows); Dd = pd.DataFrame(diffs); T.to_csv(EXP/"results"/"y04_A.csv", index=False); Dd.to_csv(EXP/"results"/"y04_A_diffs.csv", index=False)
    P("\n=== D0 형질 × 모달리티 (계통 전체 L5_all_m, 단일) τ* ===")
    P("형질    " + " ".join(f"{m:>9s}" for m in MODS))
    for t in TRAITS: P(f"{t:7s} " + " ".join(f"{T[(T.model == 'D0') & (T.trait == t) & (T['mod'] == m)].tau_star.iloc[0]:+9.2f}" for m in MODS))
    def prof(model, t, g): return T[(T.model == model) & (T.trait == t) & (T.unit == g)].set_index("mod").reindex(MODS)
    P("\n=== D2 19쌍: 같은 계통 6 모달리티 동시 모형 τ* (괄호 D2k) ===")
    P("쌍                         " + " ".join(f"{m:>15s}" for m in MODS) + "  1위(D2/D2k)")
    top, topk = [], []
    for t, g in PAIRS:
        a, k = prof("D2", t, g), prof("D2k", t, g); top.append(a.tau_star.idxmax()); topk.append(k.tau_star.idxmax())
        P(f"{t + '|' + g:26s} " + " ".join(f"{a.tau_star[m]:+6.2f}({k.tau_star[m]:+5.2f})" + ("*" if a.tau_p1[m] < .05 else " ") for m in MODS) + f"  {top[-1]}/{topk[-1]}")
    a = prof("D2", "AD", "myeloid").loc["CAGE"]; k = prof("D2k", "AD", "myeloid").loc["CAGE"]
    P(f"\n=== D-H1 AD 골수계 CAGE (D2, α 0.05): τ* {a.tau_star:+.3f} ± {a.tau_star_se:.3f}, p {a.tau_p1:.2e} → {'통과' if a.tau_p1 < .05 else '불통과'}")
    P(f"=== D-H2 AD 골수계 CAGE (D2k 트랙 맞춤, α 0.05): τ* {k.tau_star:+.3f} ± {k.tau_star_se:.3f}, p {k.tau_p1:.2e} → {'통과' if k.tau_p1 < .05 else '불통과'}")
    def max_count_p(counts, n, K=6):   # P(최다 범주 ≥ 관측) under 다항 균등, 정확 계산
        c = max(counts)
        def p_le(cmax):
            poly = [1 / math.factorial(i) for i in range(cmax + 1)]; acc = [1.0]
            for _ in range(K): acc = np.convolve(acc, poly)[: n + 1]
            return acc[n] * math.factorial(n) / K ** n if len(acc) > n else 0.0
        return 1 - p_le(c - 1)
    for lab, tp in [("D2", top), ("D2k", topk)]:
        cnt = pd.Series(tp).value_counts().reindex(MODS, fill_value=0); p = max_count_p(cnt.tolist(), len(tp))
        P(f"=== D-H3 1위 모달리티 분포 ({lab}): {cnt.to_dict()} | 최다 {cnt.idxmax()} {cnt.max()}/19, 정확 p {p:.3f}" + (" → 층위가 형질·계통마다 다르다" if p >= .05 else f" → 공통 운반 층위 후보: {cnt.idxmax()}"))
    epi, tx = ["DNASE", "H3K27ac", "H3K4me1", "H3K4me3"], ["CAGE", "RNA"]
    P("\n=== 접근성·히스톤 대 전사 (D2, 19쌍; 양의 τ* 합 비율) ===")
    for t, g in PAIRS:
        a = prof("D2", t, g).tau_star.clip(lower=0); s = a.sum(); P(f"  {t + '|' + g:26s} 접근성·히스톤 {a[epi].sum() / s if s > 0 else float('nan'):.2f} / 전사 {a[tx].sum() / s if s > 0 else float('nan'):.2f}")
    (EXP/"results"/"y04_A.txt").write_text("\n".join(log), encoding="utf-8")
else:
    rows = []; jobs = [(t, g, m) for t in TRAITS for g in LIN for m in MODS]
    P(f"단계 B: D1 S-LDSC {len(jobs)}회, 순차")
    for i, (t, g, m) in enumerate(jobs):
        a = f"L5_{g}_{m}"; R, h2, _ = run(t, "D1", a); rows.append(rec(t, "D1", g, a, R, h2, len(R) - 1) | {"mod": m})
        if (i + 1) % 72 == 0: P(f"  {i + 1}/{len(jobs)} ({(time.time() - t0) / 60:.0f}분)")
    T = pd.DataFrame(rows); T["q_bh"] = bh(T.tau_p1); T.to_csv(EXP/"results"/"y04_B.csv", index=False)
    P("\n=== D1 형질별 상위 5칸 (계통 × 모달리티, 단일; * q < 0.05) ===")
    for t in TRAITS:
        d = T[T.trait == t].sort_values("tau_star", ascending=False).head(5)
        P(f"  {t:7s} " + " > ".join(f"{r.unit}·{r['mod']} {r.tau_star:+.2f}{'*' if r.q_bh < .05 else ''}" for _, r in d.iterrows()))
    P(f"q < 0.05 칸: {int((T.q_bh < .05).sum())}/{len(T)}")
    (EXP/"results"/"y04_B.txt").write_text("\n".join(log), encoding="utf-8")
