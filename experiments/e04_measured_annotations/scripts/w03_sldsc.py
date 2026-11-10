"""실험 3e w03: S-LDSC (PROTOCOL_exp3e_2026-10-08). 형질 11개 × 모형 10개 = 110회, 순차. 호출은 실험 3d y06과 같다.
 모형(형질 t, 계통 g): BM = M_g 셋 / J1 = BM + L5_g_CAGE (주) / S_m = M_g_m + L5_g_m / M1_m = M_g_m / J1k = BM + L5k_g_CAGE / J1p = BMp + L5_g_CAGE
 τ* = τ · sd(MAF ≥ 5 %) · M_5_50 / h²_obs.  S_m 의 Atlas − 실측 차이 SE 는 part_delete 로 계산한다.
 시험: W03_DRYRUN=1 → AD 의 Z를 섞은 가짜 형질로 J1 한 번만 (결과 아님, 파싱·형식 확인용). exp3e/test/ 에 쓴다.
 재실행 안전: .results 가 있으면 건너뛴다.
출력: exp3e/sldsc/<형질>__<모형>.*, results/w03_sldsc.csv, results/w03_diffs.csv, results/w03_sldsc.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess, os, re, time
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); Y = Path(f"{DATA}/exp3d"); W = Path(f"{DATA}/exp3e")
BV = Path(f"{DATA}/tools/ldsc_bv"); B = X/"GRCh38"; SS = Path(f"{DATA}/exp3b/sumstats")
DRY = os.environ.get("W03_DRYRUN") == "1"; SL = W/"test"/"sldsc" if DRY else W/"sldsc"; SL.mkdir(parents=True, exist_ok=True)
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
TRAITS = {"AD": X/"bellenguez_noAPOE.sumstats.gz", "EA": SS/"PASS_Years_of_Education2.sumstats", "SCZ": SS/"PASS_Schizophrenia.sumstats", "BIP": SS/"PASS_Bipolar_Disorder.sumstats",
          "NEURO": SS/"PASS_Neuroticism.sumstats", "BMI": SS/"PASS_BMI1.sumstats", "MONO": Path(f"{DATA}/exp3b/monocyte_count.sumstats.gz"), "IBD": SS/"PASS_IBD.sumstats",
          "RA": SS/"PASS_Rheumatoid_Arthritis.sumstats", "HEIGHT": SS/"PASS_Height1.sumstats", "LDL": SS/"PASS_LDL.sumstats"}
LINE = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
        "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal"}   # 3d-A와 같다 (프로토콜 고정)
MODS = ["CAGE", "DNASE", "H3K27ac"]
if DRY:
    fake = W/"test"/"AD_shuffledZ.sumstats.gz"
    if not fake.exists():
        d = pd.read_csv(TRAITS["AD"], sep="\t"); d["Z"] = np.random.default_rng(20261008).permutation(d.Z.to_numpy()); d.to_csv(fake, sep="\t", index=False, compression="gzip")
    TRAITS = {"FAKE": fake}; LINE = {"FAKE": "myeloid"}
SD = pd.concat([pd.read_csv(W/"measured_sd_common.csv", index_col=0).sd_common, pd.read_csv(Y/"annot_sd_common.csv", index_col=0).sd_common])
M_TOT = sum(float(open(B/"baselineLD_v2.2"/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
def ref(name): return str((Y if name.startswith("L5") else W)/"annot"/name/f"{name}.")
def models(g):
    m = {"BM": [f"BM_{g}"], "J1": [f"BM_{g}", f"L5_{g}_CAGE"], "J1k": [f"BM_{g}", f"L5k_{g}_CAGE"], "J1p": [f"BMp_{g}", f"L5_{g}_CAGE"]}
    for x in MODS: m[f"S_{x}"] = [f"M_{g}_{x}", f"L5_{g}_{x}"]; m[f"M1_{x}"] = [f"M_{g}_{x}"]
    return {k: m[k] for k in (["J1"] if DRY else ["BM", "J1", "S_CAGE", "S_DNASE", "S_H3K27ac", "M1_CAGE", "M1_DNASE", "M1_H3K27ac", "J1k", "J1p"])}
def run(t, model, refs):
    out = SL/f"{t}__{model}"
    if not Path(f"{out}.results").exists():
        cmd = ["python", "-W", "ignore", "-c", "import sys; sys.argv=['ldsc']+sys.argv[1:]; from ldscore.ldsc import main; main()",
               "--h2", str(TRAITS[t]), "--ref-ld-chr", ",".join([str(B/"baselineLD_v2.2"/"baselineLD.")] + [ref(r) for r in refs]), "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."),
               "--overlap-annot", "--frqfile-chr", str(B/"plink_files"/"1000G.EUR.hg38."), "--print-coefficients", "--print-delete-vals", "--out", str(out)]
        subprocess.run(cmd, cwd=BV, env=env, capture_output=True, text=True)
        if not Path(f"{out}.results").exists(): raise SystemExit(f"{t} {model} 실패: {open(str(out) + '.log', encoding='utf-8').read()[-1500:]}")
    return out
rows, diffs = [], []; t0 = time.time(); jobs = [(t, k, v) for t in TRAITS for k, v in models(LINE[t]).items()]
P(f"S-LDSC {len(jobs)}회, 순차" + (" — 시험 실행(Z 섞음, 결과 아님)" if DRY else ""))
for i, (t, model, refs) in enumerate(jobs):
    out = run(t, model, refs); g = LINE[t]
    R = pd.read_csv(f"{out}.results", sep="\t"); R["name"] = R.Category.str.replace(r"L2_\d+$", "", regex=True)
    h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log", encoding="utf-8").read()).group(1))
    D = pd.read_csv(f"{out}.part_delete", sep=r"\s+", header=None).to_numpy(); nb = D.shape[0]
    added = [n for n in R.name if n.startswith(("M_", "P_", "L5_", "L5k_"))]
    assert len(added) == sum(1 if not r.startswith(("BM_", "BMp_")) else 3 for r in refs), f"{t} {model}: 추가 주석 수가 맞지 않다 {added}"
    dv = {}
    for n in added:
        k = int(R.index[R.name == n][0]); f = SD[n] * M_TOT / h2; dv[n] = D[:, k] * f
        rows.append(dict(trait=t, lineage=g, model=model, annot=n, kind="atlas" if n.startswith("L5") else "measured", h2_obs=h2, prop_h2=R["Prop._h2"][k], enrichment=R.Enrichment[k],
                         tau_star=R.Coefficient[k] * f, tau_star_se=R.Coefficient_std_error[k] * f, tau_p1=stats.norm.sf(R.Coefficient[k] / R.Coefficient_std_error[k]),
                         tau_p2=2 * stats.norm.sf(abs(R.Coefficient[k] / R.Coefficient_std_error[k]))))
    if model.startswith("S_"):
        a, m = [n for n in added if n.startswith("L5")][0], [n for n in added if n.startswith("M_")][0]
        ea = next(r["tau_star"] for r in rows[-2:] if r["annot"] == a); em = next(r["tau_star"] for r in rows[-2:] if r["annot"] == m)
        d = dv[a] - dv[m]; est = ea - em
        se = np.sqrt((nb - 1) / nb * ((d - d.mean()) ** 2).sum()); diffs.append(dict(trait=t, lineage=g, model=model, atlas=a, measured=m, diff=est, se=se, p2=2 * stats.norm.sf(abs(est / se))))
    el = (time.time() - t0) / 60; P(f"  [{i + 1}/{len(jobs)}] {t} {model} 완료 ({el:.0f}분, 남은 약 {el / (i + 1) * (len(jobs) - i - 1):.0f}분)")
T = pd.DataFrame(rows); Dd = pd.DataFrame(diffs)
tag = "_dryrun" if DRY else ""
T.to_csv(EXP/"results"/f"w03_sldsc{tag}.csv", index=False); Dd.to_csv(EXP/"results"/f"w03_diffs{tag}.csv", index=False)
P("\n" + T[["trait", "model", "annot", "tau_star", "tau_star_se", "tau_p1"]].to_string(index=False, float_format=lambda v: f"{v:.3g}"))
(EXP/"results"/f"w03_sldsc{tag}.txt").write_text("\n".join(log), encoding="utf-8")
