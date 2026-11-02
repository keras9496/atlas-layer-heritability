"""실험 3d-A y06: Atlas CAGE·RNA 대 eQTL·분자 QTL·FANTOM5 인핸서 τ* 비교 (PROTOCOL_exp3d_A_eqtl_2026-10-04, 사후 기술).
 모형(형질마다): B0 = baselineLD만(새 실행) / D1_CAGE, D1_RNA = 기존 3d D1 / D1k_CAGE, D1k_RNA = 트랙 맞춤 단일(새 실행). 순차.
 τ* = τ · sd · M_5_50 / h²_obs. baseline 주석 sd 는 MAF ≥ 5 % 참조 SNP에서 annot 파일로 계산한다.
 같은 모형 안의 Atlas − eQTL 차이 SE 는 part_delete 로 계산한다.
출력: results/y06_eqtl_compare.csv, results/y06_eqtl_compare.txt, exp3d/baseline_sd_common.csv
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess, os, re
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); Y = Path(f"{DATA}/exp3d"); AN = Y/"annot"; SL = Y/"sldsc"
BV = Path(f"{DATA}/tools/ldsc_bv"); B = X/"GRCh38"; SS = Path(f"{DATA}/exp3b/sumstats")
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
TRAITS = {"AD": X/"bellenguez_noAPOE.sumstats.gz", "EA": SS/"PASS_Years_of_Education2.sumstats", "SCZ": SS/"PASS_Schizophrenia.sumstats", "BIP": SS/"PASS_Bipolar_Disorder.sumstats",
          "NEURO": SS/"PASS_Neuroticism.sumstats", "BMI": SS/"PASS_BMI1.sumstats", "MONO": Path(f"{DATA}/exp3b/monocyte_count.sumstats.gz"), "IBD": SS/"PASS_IBD.sumstats",
          "RA": SS/"PASS_Rheumatoid_Arthritis.sumstats", "HEIGHT": SS/"PASS_Height1.sumstats", "LDL": SS/"PASS_LDL.sumstats"}
LINE = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
        "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal"}   # 프로토콜 고정
CMP = ["GTEx_eQTL_MaxCPP", "BLUEPRINT_H3K27acQTL_MaxCPP", "BLUEPRINT_H3K4me1QTL_MaxCPP", "BLUEPRINT_DNA_methylation_MaxCPP", "Enhancer_Andersson", "TSS_Hoffman", "Promoter_UCSC"]
M_TOT = sum(float(open(B/"baselineLD_v2.2"/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
# baseline 주석 sd (MAF ≥ 5 %)
sdf = Y/"baseline_sd_common.csv"
if not sdf.exists():
    acc = []
    for c in range(1, 23):
        a = pd.read_csv(B/"baselineLD_v2.2"/f"baselineLD.{c}.annot.gz", sep="\t", usecols=["SNP"] + CMP)
        f = pd.read_csv(B/"plink_files"/f"1000G.EUR.hg38.{c}.frq", sep=r"\s+", usecols=["SNP", "MAF"])
        acc.append(a.merge(f, on="SNP")[lambda d: d.MAF >= .05][CMP])
    pd.concat(acc).std(ddof=0).rename("sd_common").to_csv(sdf)
SD = pd.concat([pd.read_csv(sdf, index_col=0).sd_common, pd.read_csv(Y/"annot_sd_common.csv", index_col=0).sd_common])
def run(t, model, annot):
    out = SL/f"{t}__{model}__{annot or 'baseline'}"
    if not Path(f"{out}.results").exists():
        ref = ",".join([str(B/"baselineLD_v2.2"/"baselineLD.")] + ([str(AN/annot/f"{annot}.")] if annot else []))
        cmd = ["python", "-W", "ignore", "-c", "import sys; sys.argv=['ldsc']+sys.argv[1:]; from ldscore.ldsc import main; main()",
               "--h2", str(TRAITS[t]), "--ref-ld-chr", ref, "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."), "--overlap-annot",
               "--frqfile-chr", str(B/"plink_files"/"1000G.EUR.hg38."), "--print-coefficients", "--print-delete-vals", "--out", str(out)]
        subprocess.run(cmd, cwd=BV, env=env, capture_output=True, text=True)
        if not Path(f"{out}.results").exists(): raise SystemExit(f"{t} {model} 실패: {open(str(out) + '.log', encoding='utf-8').read()[-1200:]}")
    return out
rows = []
for t in TRAITS:
    g = LINE[t]
    for model, annot in [("B0", None), ("D1", f"L5_{g}_CAGE"), ("D1", f"L5_{g}_RNA"), ("D1k", f"L5k_{g}_CAGE"), ("D1k", f"L5k_{g}_RNA")]:
        out = run(t, model, annot)
        R = pd.read_csv(f"{out}.results", sep="\t"); R["name"] = R.Category.str.replace(r"L2_\d+$", "", regex=True)
        h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log", encoding="utf-8").read()).group(1))
        D = pd.read_csv(f"{out}.part_delete", sep=r"\s+", header=None).to_numpy(); nb = D.shape[0]
        def ts(n):
            i = R.index[R.name == n][0]; f = SD[n] * M_TOT / h2
            return R.Coefficient[i] * f, R.Coefficient_std_error[i] * f, D[:, i] * f
        rec = dict(trait=t, lineage=g, model=model, atlas=annot or "", h2_obs=h2)
        for n in CMP: v, s, _ = ts(n); rec[f"{n}"] = v; rec[f"{n}_se"] = s
        if annot:
            v, s, da = ts(annot); rec["atlas_tau"] = v; rec["atlas_se"] = s
            _, _, de = ts("GTEx_eQTL_MaxCPP"); d = da - de; rec["atlas_minus_eqtl"] = v - rec["GTEx_eQTL_MaxCPP"]
            rec["diff_se"] = np.sqrt((nb - 1) / nb * ((d - d.mean()) ** 2).sum()); rec["diff_p2"] = 2 * stats.norm.sf(abs(rec["atlas_minus_eqtl"] / rec["diff_se"]))
        rows.append(rec)
T = pd.DataFrame(rows); T.to_csv(EXP/"results"/"y06_eqtl_compare.csv", index=False)
P("=== 형질별: 같은 모형 안 τ* (± SE) ===")
P("형질    계통        모형  Atlas 주석            Atlas τ*        eQTL τ*         Atlas−eQTL (p)      FANTOM5 인핸서  B0 eQTL τ*")
for t in TRAITS:
    b0 = T[(T.trait == t) & (T.model == "B0")].iloc[0]
    for _, r in T[(T.trait == t) & (T.model != "B0")].iterrows():
        P(f"{t:7s} {r.lineage:11s} {r.model:4s}  {r.atlas:20s} {r.atlas_tau:+5.2f} ± {r.atlas_se:.2f}   {r.GTEx_eQTL_MaxCPP:+5.2f} ± {r.GTEx_eQTL_MaxCPP_se:.2f}   "
          f"{r.atlas_minus_eqtl:+5.2f} ({r.diff_p2:.2g})   {r.Enhancer_Andersson:+5.2f}         {b0.GTEx_eQTL_MaxCPP:+5.2f}")
P("\n=== 요약 (CAGE, 프로토콜의 읽는 규칙) ===")
for model in ["D1", "D1k"]:
    c = T[(T.model == model) & T.atlas.str.endswith("CAGE")]
    big = c[(c.atlas_minus_eqtl > 0) & (c.diff_p2 < .05)]
    P(f"  {model}: Atlas CAGE τ* > eQTL τ* (같은 모형, 차이 p < 0.05) {len(big)}/{len(c)} 형질: {', '.join(big.trait)}"
      f" | 점추정만 큰 형질 {int((c.atlas_minus_eqtl > 0).sum())}/{len(c)}")
P("\n=== eQTL 감쇠: B0 → +Atlas CAGE (D1) ===")
for t in TRAITS:
    b0 = T[(T.trait == t) & (T.model == "B0")].iloc[0].GTEx_eQTL_MaxCPP; c1 = T[(T.trait == t) & (T.model == "D1") & T.atlas.str.endswith("CAGE")].iloc[0].GTEx_eQTL_MaxCPP
    P(f"  {t:7s} eQTL τ* {b0:+.2f} → {c1:+.2f} ({(c1 - b0) / abs(b0) * 100 if abs(b0) > 1e-9 else float('nan'):+.0f} %)")
(EXP/"results"/"y06_eqtl_compare.txt").write_text("\n".join(log), encoding="utf-8")
