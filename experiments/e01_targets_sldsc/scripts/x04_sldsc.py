"""실험 3 x04: S-LDSC (brainvisa/ldscore 포크 회귀; baselineLD v2.2 hg38 조건부, --overlap-annot, --frqfile MAF≥5 %, 가중치 hm3_noMHC).
 단일 주석 51개 + ct5set(세포형 5개 동시; H2) + AVI 조건부(avi5 + ct5_<c>).
 τ* = τ · sd(주석, 흔한 SNP) · M_5_50(전체) / h²_obs. H1·H3: τ > 0 단측 z. H2: part_delete 로 τ*_골수계 − τ*_뉴런 jackknife SE.
출력: <DATA>/exp3/sldsc/<이름>.*, results/x04_sldsc.csv, x04_sldsc.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess, os, re, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); AN = X/"annot"; SL = X/"sldsc"; SL.mkdir(exist_ok=True)
BV = Path(f"{DATA}/tools/ldsc_bv"); B = X/"GRCh38"
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
CTX = ["myeloid", "neuron", "astrocyte", "liver", "nonliver"]
SD = pd.read_csv(X/"annot_sd_common.csv", index_col=0).sd_common
M_TOT = sum(float(open(B/"baselineLD_v2.2"/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
def run(name, extra=()):
    out = SL/name
    if not Path(f"{out}.results").exists():
        ref = ",".join([str(B/"baselineLD_v2.2"/"baselineLD.")] + [str(AN/e/f"{e}.") for e in extra] + [str(AN/name.split("__")[0]/f"{name.split('__')[0]}.")])
        cmd = ["python", "-W", "ignore", "-c", "import sys; sys.argv=['ldsc']+sys.argv[1:]; from ldscore.ldsc import main; main()",
               "--h2", str(X/"bellenguez_noAPOE.sumstats.gz"), "--ref-ld-chr", ref, "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."), "--overlap-annot",
               "--frqfile-chr", str(B/"plink_files"/"1000G.EUR.hg38."), "--print-coefficients", "--print-delete-vals", "--out", str(out)]
        subprocess.run(cmd, cwd=BV, env=env, capture_output=True, text=True)
        if not Path(f"{out}.results").exists(): raise SystemExit(f"{name} 실패: {open(str(out) + '.log').read()[-1500:]}")
    R = pd.read_csv(f"{out}.results", sep="\t"); h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log").read()).group(1))
    return R, h2, out
rows = []
def record(label, annot, R, h2, k):
    r = R.iloc[k]; tau, se = r.Coefficient, r.Coefficient_std_error; ts = tau * SD[annot] * M_TOT / h2
    rows.append(dict(analysis=label, annot=annot, prop_snps=r["Prop._SNPs"], prop_h2=r["Prop._h2"], enrichment=r.Enrichment, enrich_se=r.Enrichment_std_error, enrich_p=r.Enrichment_p,
                     tau=tau, tau_se=se, tau_star=ts, tau_star_se=se * SD[annot] * M_TOT / h2, tau_z=tau / se, tau_p1=stats.norm.sf(tau / se), h2_obs=h2))
names = [f"ct5_{c}" for c in CTX] + [f"mod_{c}_{m}" for c in CTX for m in ["DNASE", "H3K27ac", "H3K4me1", "H3K4me3", "CAGE", "RNA"]] + ["avi5"] + \
        [f"{t}_{c}" for t in ["ct1", "ct10", "ctc"] for c in CTX]
for n in names:
    R, h2, _ = run(n); record("single|baselineLD", n, R, h2, len(R) - 1); r = rows[-1]
    P(f"{n:26s} 농축 {r['enrichment']:5.2f} ± {r['enrich_se']:.2f} | τ* {r['tau_star']:+.3f} ± {r['tau_star_se']:.3f} (단측 p {r['tau_p1']:.2e})")
# H2: 세포형 5개 동시
R, h2, out = run("ct5set"); k0 = len(R) - 5
for j, c in enumerate(CTX): record("joint|ct5set", f"ct5_{c}", R, h2, k0 + j)
D = pd.read_csv(f"{out}.part_delete", sep=r"\s+", header=None).to_numpy(); nb = D.shape[0]
tsd = lambda j: D[:, k0 + j] * SD[f"ct5_{CTX[j]}"] * M_TOT / h2
diffs = {}
for c2 in ["neuron", "astrocyte", "nonliver", "liver"]:
    d = tsd(0) - tsd(CTX.index(c2)); est = rows[-5]["tau_star"] - rows[-5 + CTX.index(c2)]["tau_star"]; se = np.sqrt((nb - 1) / nb * ((d - d.mean()) ** 2).sum())
    diffs[c2] = (est, se, stats.norm.sf(est / se)); P(f"[H2 동시 모형] τ*_골수계 − τ*_{c2}: {est:+.3f} ± {se:.3f}, 단측 p {stats.norm.sf(est / se):.4f}")
# AVI 조건부
for c in CTX:
    R, h2, _ = run(f"ct5_{c}__avi", extra=("avi5",)); record("single|baselineLD+avi5", f"ct5_{c}", R, h2, len(R) - 1)
T = pd.DataFrame(rows); T.to_csv(EXP/"results"/"x04_sldsc.csv", index=False)
P("\n=== E3-H1 (세포형, baselineLD 조건부; α 0.01 단측) ===")
for c in CTX:
    r = T[(T.analysis == "single|baselineLD") & (T.annot == f"ct5_{c}")].iloc[0]
    P(f"  {c:9s} 농축 {r.enrichment:.2f} (p {r.enrich_p:.1e}) | τ* {r.tau_star:+.3f} ± {r.tau_star_se:.3f}, p {r.tau_p1:.2e} → {'통과' if r.tau_p1 < .01 else '불통과'}")
P(f"=== E3-H2 τ*_골수계 − τ*_뉴런 (동시 모형): {diffs['neuron'][0]:+.3f} ± {diffs['neuron'][1]:.3f}, p {diffs['neuron'][2]:.4f} → {'통과' if diffs['neuron'][2] < .05 else '불통과'}")
P("=== E3-H3 모달리티 (α 0.05/30 단측) ===")
mm = T[(T.analysis == "single|baselineLD") & T.annot.str.startswith("mod_")]
for _, r in mm.sort_values("tau_p1").iterrows(): P(f"  {r.annot:26s} τ* {r.tau_star:+.3f} ± {r.tau_star_se:.3f}, p {r.tau_p1:.2e}{'  ← 통과' if r.tau_p1 < .05 / 30 else ''} | 농축 {r.enrichment:.2f}")
(EXP/"results"/"x04_sldsc.txt").write_text("\n".join(log), encoding="utf-8")
