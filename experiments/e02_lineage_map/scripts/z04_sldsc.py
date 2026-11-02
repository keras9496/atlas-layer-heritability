"""실험 3c z04: 형질 11개 × 계통 12개 S-LDSC (PROTOCOL_exp3c_2026-09-28). 순차 실행(사용자 선호). 실험 3 x04와 같은 호출.
 M0: baselineLD + L5_all / M1: baselineLD + L5_g / M2: baselineLD + L5_all + L5_g  → 11 + 132 + 132 = 275회
 τ* = τ · sd(주석, 흔한 SNP) · M_5_50 / h²_obs. 단측 p(τ > 0). BH FDR 은 모형마다 132칸에서.
 판정: M-H1 (AD myeloid, M1·M2, α 0.05/12), M-H2 (뇌 형질 5개 CNS_tissue, M1, 형질마다 α 0.05/12; 과반 통과 기준), M-H3 (AD CNS_tissue, M1, α 0.05/12)
 재실행 안전: .results 가 있으면 건너뛴다.
출력: exp3c/sldsc/<형질>__<모형>__<계통>.*, results/z04_sldsc.csv, results/z04_sldsc.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess, os, re, time
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); Y = Path(f"{DATA}/exp3c"); AN = Y/"annot"; SL = Y/"sldsc"; SL.mkdir(exist_ok=True)
BV = Path(f"{DATA}/tools/ldsc_bv"); B = X/"GRCh38"; SS = Path(f"{DATA}/exp3b/sumstats")
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LIN = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
TRAITS = {"AD": X/"bellenguez_noAPOE.sumstats.gz", "EA": SS/"PASS_Years_of_Education2.sumstats", "SCZ": SS/"PASS_Schizophrenia.sumstats", "BIP": SS/"PASS_Bipolar_Disorder.sumstats",
          "NEURO": SS/"PASS_Neuroticism.sumstats", "BMI": SS/"PASS_BMI1.sumstats", "MONO": Path(f"{DATA}/exp3b/monocyte_count.sumstats.gz"), "IBD": SS/"PASS_IBD.sumstats",
          "RA": SS/"PASS_Rheumatoid_Arthritis.sumstats", "HEIGHT": SS/"PASS_Height1.sumstats", "LDL": SS/"PASS_LDL.sumstats"}
BRAIN = ["EA", "SCZ", "BIP", "NEURO", "BMI"]
SD = pd.read_csv(Y/"annot_sd_common.csv", index_col=0).sd_common
M_TOT = sum(float(open(B/"baselineLD_v2.2"/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
def run(t, model, g):
    out = SL/f"{t}__{model}__{g}"
    if not Path(f"{out}.results").exists():
        ann = {"M0": ["L5_all"], "M1": [f"L5_{g}"], "M2": ["L5_all", f"L5_{g}"]}[model]
        ref = ",".join([str(B/"baselineLD_v2.2"/"baselineLD.")] + [str(AN/a/f"{a}.") for a in ann])
        cmd = ["python", "-W", "ignore", "-c", "import sys; sys.argv=['ldsc']+sys.argv[1:]; from ldscore.ldsc import main; main()",
               "--h2", str(TRAITS[t]), "--ref-ld-chr", ref, "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."), "--overlap-annot",
               "--frqfile-chr", str(B/"plink_files"/"1000G.EUR.hg38."), "--print-coefficients", "--print-delete-vals", "--out", str(out)]
        subprocess.run(cmd, cwd=BV, env=env, capture_output=True, text=True)
        if not Path(f"{out}.results").exists(): raise SystemExit(f"{t} {model} {g} 실패: {open(str(out) + '.log', encoding='utf-8').read()[-1200:]}")
    R = pd.read_csv(f"{out}.results", sep="\t"); lg = open(f"{out}.log", encoding="utf-8").read()
    h2, h2se = map(float, re.search(r"Total Observed scale h2: ([\d.eE-]+) \(([\d.eE-]+)\)", lg).groups())
    r = R.iloc[-1]; a = "L5_all" if model == "M0" else f"L5_{g}"; f = SD[a] * M_TOT / h2
    return dict(trait=t, model=model, lineage=g, h2_obs=h2, h2_se=h2se, prop_h2=r["Prop._h2"], enrichment=r.Enrichment, enrich_se=r.Enrichment_std_error,
                enrich_p=r.Enrichment_p, tau_star=r.Coefficient * f, tau_star_se=r.Coefficient_std_error * f, tau_z=r.Coefficient / r.Coefficient_std_error,
                tau_p1=stats.norm.sf(r.Coefficient / r.Coefficient_std_error))
jobs = [(t, "M0", "all") for t in TRAITS] + [(t, m, g) for m in ["M1", "M2"] for t in TRAITS for g in LIN]
P(f"S-LDSC {len(jobs)}회, 순차 실행")
rows = []; t0 = time.time()
for i, j in enumerate(jobs):
    rows.append(run(*j))
    if (i + 1) % 11 == 0: P(f"  {i + 1}/{len(jobs)} ({(time.time() - t0) / 60:.0f}분)")
T = pd.DataFrame(rows)
def bh(p):
    p = np.asarray(p); o = np.argsort(p); q = p[o] * len(p) / np.arange(1, len(p) + 1); q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty_like(q); out[o] = np.minimum(q, 1); return out
T["q_bh"] = np.nan
for m in ["M1", "M2"]: k = T.model == m; T.loc[k, "q_bh"] = bh(T.loc[k, "tau_p1"])
T.to_csv(EXP/"results"/"z04_sldsc.csv", index=False)
def cell(t, m, g): return T[(T.trait == t) & (T.model == m) & (T.lineage == g)].iloc[0]
for m, name in [("M1", "농축 지도 (baselineLD + L5_g)"), ("M2", "계통 특이 지도 (baselineLD + L5_all + L5_g)")]:
    P(f"\n=== {m} {name}: τ* (* q < 0.05) ===")
    P("형질    " + " ".join(f"{g[:9]:>9s}" for g in LIN))
    for t in TRAITS:
        P(f"{t:7s} " + " ".join(f"{cell(t, m, g).tau_star:+8.2f}{'*' if cell(t, m, g).q_bh < .05 else ' '}" for g in LIN))
P("\n=== M0 L5_all (Atlas 전체) ===")
for t in TRAITS: r = cell(t, "M0", "all"); P(f"  {t:7s} h²_obs {r.h2_obs:.4f} | 농축 {r.enrichment:.2f} | τ* {r.tau_star:+.2f} ± {r.tau_star_se:.2f} (p {r.tau_p1:.1e})")
P("\n=== 형질별 계통 순위 (M1 τ*, 상위 3) ===")
for t in TRAITS:
    d = T[(T.trait == t) & (T.model == "M1")].sort_values("tau_star", ascending=False).head(3)
    P(f"  {t:7s} " + " > ".join(f"{r.lineage} {r.tau_star:+.2f} (p {r.tau_p1:.1e})" for _, r in d.iterrows()))
a = 0.05 / 12
P("\n=== M-H1 AD myeloid (α 0.05/12) ===")
for m in ["M1", "M2"]: r = cell("AD", m, "myeloid"); P(f"  {m}: τ* {r.tau_star:+.3f} ± {r.tau_star_se:.3f}, p {r.tau_p1:.2e} → {'통과' if r.tau_p1 < a else '불통과'}")
P("=== M-H2 뇌 형질 CNS_tissue (M1, 형질마다 α 0.05/12; 과반 기준) ===")
npass = 0
for t in BRAIN:
    r = cell(t, "M1", "CNS_tissue"); ok = r.tau_p1 < a; npass += ok
    P(f"  {t:6s} τ* {r.tau_star:+.3f} ± {r.tau_star_se:.3f}, p {r.tau_p1:.2e} → {'통과' if ok else '불통과'} | 참고 neural_cells τ* {cell(t, 'M1', 'neural_cells').tau_star:+.3f} (p {cell(t, 'M1', 'neural_cells').tau_p1:.1e})")
P(f"  → {npass}/5 통과: {'과반 → Atlas 뇌 조직 주석 작동' if npass >= 3 else '과반 미달 → M-H3 해석 안 함'}")
r = cell("AD", "M1", "CNS_tissue")
P(f"=== M-H3 AD CNS_tissue (M1): τ* {r.tau_star:+.3f} ± {r.tau_star_se:.3f}, p {r.tau_p1:.2e} → {'통과' if r.tau_p1 < a else '불통과'}" + ("" if npass >= 3 else " (M-H2 과반 미달이라 해석하지 않음)"))
(EXP/"results"/"z04_sldsc.txt").write_text("\n".join(log), encoding="utf-8")
