"""실험 3b x07: S-LDSC 강건성·주석 타당성 (PROTOCOL_exp3b_2026-09-28). 실험 3 x04 와 같은 호출(baselineLD v2.2 조건부).
 분석 1: 대형 좌위 제외(noGWS; 부: noBIN1) → E3S-1(ct5_myeloid τ* > 0, α 0.05), E3S-2(동시 모형 골수계 − 뉴런 > 0, α 0.05)
         + 실험 3 주 실행(ct5_myeloid)의 블록 하나 빼기 τ* (회귀 SNP 순서 재구성으로 블록 구간 복원)
 분석 2: APOE 포함(withAPOE), 기술
 분석 3: 대조 형질 6개 × (ct5 5개 + ct5set) → C-H1(골수계 형질 3개, α 0.05/3), C-H2(뉴런 형질 2개, α 0.05/2)
 재실행 안전: .results 가 있으면 다시 돌리지 않는다. 요약통계가 없는 형질은 건너뛴다.
출력: <DATA>/exp3b/sldsc/<세트>__<주석>.*, results/x07_3b.csv, results/x07_3b.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess, os, re
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
from concurrent.futures import ThreadPoolExecutor
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); Y = Path(f"{DATA}/exp3b"); AN = X/"annot"; SL = Y/"sldsc"; SL.mkdir(exist_ok=True, parents=True)
BV = Path(f"{DATA}/tools/ldsc_bv"); B = X/"GRCh38"; SS = Y/"sumstats"
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
CTX = ["myeloid", "neuron", "astrocyte", "liver", "nonliver"]
SD = pd.read_csv(X/"annot_sd_common.csv", index_col=0).sd_common
M_TOT = sum(float(open(B/"baselineLD_v2.2"/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
SETS = {  # 이름: (요약통계, 역할)
    "noGWS": (Y/"bellenguez_noAPOE_noGWS.sumstats.gz", "AD 대형 좌위 제외"),
    "noBIN1": (Y/"bellenguez_noAPOE_noBIN1.sumstats.gz", "AD BIN1 제외"),
    "withAPOE": (Y/"bellenguez_withAPOE.sumstats.gz", "AD APOE 포함"),
    "MONO": (Y/"monocyte_count.sumstats.gz", "골수계 대조"), "IBD": (SS/"PASS_IBD.sumstats", "골수계 대조"), "RA": (SS/"PASS_Rheumatoid_Arthritis.sumstats", "골수계 대조"),
    "SCZ": (SS/"PASS_Schizophrenia.sumstats", "뉴런 대조"), "EA": (SS/"PASS_Years_of_Education2.sumstats", "뉴런 대조"), "HEIGHT": (SS/"PASS_Height1.sumstats", "기준"),
}
AD_SETS = ["noGWS", "noBIN1", "withAPOE"]; EXTRA = ["mod_myeloid_CAGE", "mod_myeloid_DNASE"]
def job_list(s): return [f"ct5_{c}" for c in CTX] + ["ct5set"] + (EXTRA if s in AD_SETS else [])
def run(s, annot):
    out = SL/f"{s}__{annot}"
    if not Path(f"{out}.results").exists():
        ref = ",".join([str(B/"baselineLD_v2.2"/"baselineLD."), str(AN/annot/f"{annot}.")])
        cmd = ["python", "-W", "ignore", "-c", "import sys; sys.argv=['ldsc']+sys.argv[1:]; from ldscore.ldsc import main; main()",
               "--h2", str(SETS[s][0]), "--ref-ld-chr", ref, "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."), "--overlap-annot",
               "--frqfile-chr", str(B/"plink_files"/"1000G.EUR.hg38."), "--print-coefficients", "--print-delete-vals", "--out", str(out)]
        subprocess.run(cmd, cwd=BV, env=env, capture_output=True, text=True)
        if not Path(f"{out}.results").exists(): return (s, annot, f"실패: {open(str(out) + '.log', encoding='utf-8').read()[-800:]}")
    return (s, annot, "ok")
jobs = [(s, a) for s in SETS if SETS[s][0].exists() for a in job_list(s)]
skipped = [s for s in SETS if not SETS[s][0].exists()]
P(f"실행 {len(jobs)}개 (6개 병렬). 요약통계 없어 건너뜀: {skipped or '없음'}")
with ThreadPoolExecutor(6) as ex:
    for s, a, st in ex.map(lambda j: run(*j), jobs):
        if st != "ok": P(f"  {s} {a} {st}")

rows = []
def read(s, annot):
    out = SL/f"{s}__{annot}"; R = pd.read_csv(f"{out}.results", sep="\t"); t = open(f"{out}.log", encoding="utf-8").read()
    h2, h2se = map(float, re.search(r"Total Observed scale h2: ([\d.eE-]+) \(([\d.eE-]+)\)", t).groups())
    nreg = int(re.findall(r"\((\d+) SNPs remain\)", t)[-1]); return R, h2, h2se, nreg, out
def record(s, label, annot, R, h2, h2se, nreg, k):
    r = R.iloc[k]; tau, se = r.Coefficient, r.Coefficient_std_error; f = SD[annot] * M_TOT / h2
    rows.append(dict(set=s, role=SETS[s][1], analysis=label, annot=annot, n_reg_snps=nreg, h2_obs=h2, h2_se=h2se, prop_h2=r["Prop._h2"], enrichment=r.Enrichment,
                     enrich_se=r.Enrichment_std_error, enrich_p=r.Enrichment_p, tau_star=tau * f, tau_star_se=se * f, tau_z=tau / se, tau_p1=stats.norm.sf(tau / se)))
diffs = []
for s in SETS:
    if s in skipped: continue
    for a in job_list(s):
        if a == "ct5set": continue
        R, h2, h2se, nreg, _ = read(s, a); record(s, "single", a, R, h2, h2se, nreg, len(R) - 1)
    R, h2, h2se, nreg, out = read(s, "ct5set"); k0 = len(R) - 5
    for j, c in enumerate(CTX): record(s, "joint", f"ct5_{c}", R, h2, h2se, nreg, k0 + j)
    D = pd.read_csv(f"{out}.part_delete", sep=r"\s+", header=None).to_numpy(); nb = D.shape[0]
    J = {r["annot"]: r for r in rows if r["set"] == s and r["analysis"] == "joint"}
    for c2 in ["neuron", "astrocyte", "liver", "nonliver"]:
        j2 = CTX.index(c2); d = (D[:, k0] * SD["ct5_myeloid"] - D[:, k0 + j2] * SD[f"ct5_{c2}"]) * M_TOT / h2
        est = J["ct5_myeloid"]["tau_star"] - J[f"ct5_{c2}"]["tau_star"]; se = np.sqrt((nb - 1) / nb * ((d - d.mean()) ** 2).sum())
        diffs.append(dict(set=s, contrast=f"myeloid-{c2}", est=est, se=se, p1=stats.norm.sf(est / se)))
T = pd.DataFrame(rows); Dd = pd.DataFrame(diffs)
T.to_csv(EXP/"results"/"x07_3b.csv", index=False); Dd.to_csv(EXP/"results"/"x07_3b_diffs.csv", index=False)

# 실험 3 주 결과(비교 기준)
base = pd.read_csv(EXP/"results"/"x04_sldsc.csv"); b0 = base[(base.analysis == "single|baselineLD")].set_index("annot")
def one(s, a, lab="single"): return T[(T.set == s) & (T.annot == a) & (T.analysis == lab)].iloc[0]
def dif(s, c2="neuron"): return Dd[(Dd.set == s) & (Dd.contrast == f"myeloid-{c2}")].iloc[0]
P("\n=== 분석 1: 대형 좌위 제외 (사후 민감도) ===")
for s in ["noGWS", "noBIN1"]:
    if s in skipped: continue
    r0 = one(s, "ct5_myeloid"); P(f"[{s}] 회귀 SNP {r0.n_reg_snps:,}, h²_obs {r0.h2_obs:.4f} ({r0.h2_se:.4f})")
    for a in [f"ct5_{c}" for c in CTX] + EXTRA:
        r = one(s, a); o = b0.loc[a]
        P(f"  {a:20s} 농축 {r.enrichment:6.2f} ± {r.enrich_se:.2f} | τ* {r.tau_star:+.3f} ± {r.tau_star_se:.3f} (p {r.tau_p1:.2e}) | 실험 3: τ* {o.tau_star:+.3f}, 농축 {o.enrichment:.2f} → 비율 {r.tau_star / o.tau_star:.2f}")
    for c2 in ["neuron", "astrocyte", "liver", "nonliver"]:
        d = dif(s, c2); P(f"  [동시] 골수계 − {c2:9s}: {d.est:+.3f} ± {d.se:.3f}, p {d.p1:.4f}")
if "noGWS" not in skipped:
    r = one("noGWS", "ct5_myeloid"); d = dif("noGWS")
    P(f"=== E3S-1 대형 좌위 제외 후 골수계 τ* {r.tau_star:+.3f} ± {r.tau_star_se:.3f}, p {r.tau_p1:.2e} → {'통과' if r.tau_p1 < .05 else '불통과'}")
    P(f"=== E3S-2 대형 좌위 제외 후 골수계 − 뉴런 {d.est:+.3f} ± {d.se:.3f}, p {d.p1:.4f} → {'통과' if d.p1 < .05 else '불통과'}")

# 블록 하나 빼기 (실험 3 주 실행 ct5_myeloid)
P("\n=== 블록 영향 (실험 3 주 실행 ct5_myeloid, 블록 하나 빼기 τ*) ===")
o = X/"sldsc"/"ct5_myeloid"; t = open(f"{o}.log", encoding="utf-8").read(); h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", t).group(1)); nreg = int(re.findall(r"\((\d+) SNPs remain\)", t)[-1])
PD = pd.read_csv(f"{o}.part_delete", sep=r"\s+", header=None).to_numpy(); f = SD["ct5_myeloid"] * M_TOT / h2
loo = PD[:, -1] * f; full = b0.loc["ct5_myeloid"].tau_star
ss = set(pd.read_csv(X/"bellenguez_noAPOE.sumstats.gz", sep="\t", usecols=["SNP"]).SNP)
ref = pd.concat([pd.read_csv(B/"baselineLD_v2.2"/f"baselineLD.{c}.l2.ldscore.gz", sep="\t", usecols=["CHR", "SNP", "BP"]) for c in range(1, 23)])
w = set(pd.concat([pd.read_csv(B/"weights"/f"weights.hm3_noMHC.{c}.l2.ldscore.gz", sep="\t", usecols=["SNP"]) for c in range(1, 23)]).SNP)
reg = ref[ref.SNP.isin(ss) & ref.SNP.isin(w)].drop_duplicates("SNP").reset_index(drop=True)
P(f"회귀 SNP 재구성 {len(reg):,} (로그 {nreg:,}) → {'일치' if len(reg) == nreg else '불일치: 블록 구간은 근사'}")
sep = np.floor(np.linspace(0, len(reg), PD.shape[0] + 1)).astype(int)
P(f"블록 하나 빼기 τ*: 최소 {loo.min():+.3f}, 중앙 {np.median(loo):+.3f}, 최대 {loo.max():+.3f} (전체 {full:+.3f}, SE {b0.loc['ct5_myeloid'].tau_star_se:.3f}); 모든 블록에서 > 0: {bool((loo > 0).all())}")
for i in np.argsort(loo)[:5]:
    a, z = reg.iloc[sep[i]], reg.iloc[sep[i + 1] - 1]
    P(f"  블록 {i:3d}: τ* {loo[i]:+.3f} (Δ {loo[i] - full:+.3f}) | chr{a.CHR}:{a.BP / 1e6:.1f}–" + (f"{z.BP / 1e6:.1f} Mb" if a.CHR == z.CHR else f"chr{z.CHR}:{z.BP / 1e6:.1f} Mb"))
pd.DataFrame({"block": range(len(loo)), "loo_tau_star": loo, "chr_start": [reg.iloc[sep[i]].CHR for i in range(len(loo))], "bp_start": [reg.iloc[sep[i]].BP for i in range(len(loo))],
              "chr_end": [reg.iloc[sep[i + 1] - 1].CHR for i in range(len(loo))], "bp_end": [reg.iloc[sep[i + 1] - 1].BP for i in range(len(loo))]}).to_csv(EXP/"results"/"x07_3b_block_loo.csv", index=False)

P("\n=== 분석 2: APOE 포함 (기술) ===")
if "withAPOE" not in skipped:
    r0 = one("withAPOE", "ct5_myeloid"); P(f"회귀 SNP {r0.n_reg_snps:,}, h²_obs {r0.h2_obs:.4f} ({r0.h2_se:.4f})")
    for a in [f"ct5_{c}" for c in CTX] + EXTRA:
        r = one("withAPOE", a); P(f"  {a:20s} 농축 {r.enrichment:6.2f} ± {r.enrich_se:.2f} | τ* {r.tau_star:+.3f} ± {r.tau_star_se:.3f} (p {r.tau_p1:.2e}) | 실험 3 τ* {b0.loc[a].tau_star:+.3f}")
    d = dif("withAPOE"); P(f"  [동시] 골수계 − 뉴런: {d.est:+.3f} ± {d.se:.3f}, p {d.p1:.4f}")

P("\n=== 분석 3: 대조 형질 (단일 주석, baselineLD 조건부) ===")
tr = [s for s in ["MONO", "IBD", "RA", "SCZ", "EA", "HEIGHT"] if s not in skipped]
P("형질     h²_obs  | " + " | ".join(f"{c:>18s}" for c in CTX))
for s in tr:
    P(f"{s:7s} {one(s, 'ct5_myeloid').h2_obs:.4f} | " + " | ".join(f"{one(s, f'ct5_{c}').tau_star:+.2f} (p {one(s, f'ct5_{c}').tau_p1:.1e})" for c in CTX))
P("동시 모형 골수계 − 뉴런: " + "; ".join(f"{s} {dif(s).est:+.2f} ± {dif(s).se:.2f}" for s in tr))
my = [s for s in ["MONO", "IBD", "RA"] if s in tr]; nu = [s for s in ["SCZ", "EA"] if s in tr]
P("=== C-H1 골수계 주석 타당성 (α 0.05/3) ===")
for s in my: r = one(s, "ct5_myeloid"); P(f"  {s:5s} τ* {r.tau_star:+.3f} ± {r.tau_star_se:.3f}, p {r.tau_p1:.2e} → {'통과' if r.tau_p1 < .05 / 3 else '불통과'}")
if len(my) < 3: P(f"  (미실행: {sorted(set(['MONO', 'IBD', 'RA']) - set(my))})")
P("=== C-H2 뉴런 주석 타당성 (α 0.05/2) ===")
for s in nu: r = one(s, "ct5_neuron"); P(f"  {s:5s} τ* {r.tau_star:+.3f} ± {r.tau_star_se:.3f}, p {r.tau_p1:.2e} → {'통과' if r.tau_p1 < .05 / 2 else '불통과'}")
(EXP/"results"/"x07_3b.txt").write_text("\n".join(log), encoding="utf-8")
