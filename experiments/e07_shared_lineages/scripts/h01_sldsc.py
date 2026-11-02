"""실험 3h h01: 계통 고정 비교 S-LDSC (PROTOCOL_exp3h_2026-10-13). 순차 실행. 호출은 3f v02와 같다.
 T3all_{g} = baselineLD + 3d L5_{g}_{CAGE,DNASE,H3K27ac} (전체 트랙, 3-층위 동시). 형질 × 농축 계통(계통 지도 단측 P < 0.05) 가운데 자기 계통이 아닌 것만 새로 돌린다(자기 계통은 3f T3all).
 단계: dry = Z를 섞은 가짜 형질로 지방 계통 1회(형식 확인, 결과 아님) / main = 33회
 재실행 안전: .results 가 있으면 건너뛴다.  사용: python h01_sldsc.py dry | main
출력: exp3h/sldsc[/test]/<형질>__T3all_<계통>.{results,log,part_delete}, results/h01_<단계>.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess, os, sys, time, glob
import pandas as pd
from pathlib import Path
from scipy.stats import norm
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; R0 = Path(__file__).resolve().parents[2]
X = Path(f"{DATA}/exp3"); Y = Path(f"{DATA}/exp3d"); F = Path(f"{DATA}/exp3f"); H = Path(f"{DATA}/exp3h")
BV = Path(f"{DATA}/tools/ldsc_bv"); B = X/"GRCh38"; SS = Path(f"{DATA}/exp3b/sumstats")
SL = H/"test"/"sldsc" if STAGE == "dry" else H/"sldsc"; SL.mkdir(parents=True, exist_ok=True); (EXP/"results").mkdir(exist_ok=True)
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
TRAITS = {"AD": X/"bellenguez_noAPOE.sumstats.gz", "MONO": Path(f"{DATA}/exp3b/monocyte_count.sumstats.gz"), "EA": SS/"PASS_Years_of_Education2.sumstats",
          "SCZ": SS/"PASS_Schizophrenia.sumstats", "BIP": SS/"PASS_Bipolar_Disorder.sumstats", "NEURO": SS/"PASS_Neuroticism.sumstats", "BMI": SS/"PASS_BMI1.sumstats",
          "IBD": SS/"PASS_IBD.sumstats", "RA": SS/"PASS_Rheumatoid_Arthritis.sumstats", "LDL": SS/"PASS_LDL.sumstats", "HEIGHT": SS/"PASS_Height1.sumstats",
          "BMD": F/"sumstats"/"BMD.sumstats.gz", "NEUT": F/"sumstats"/"NEUT.sumstats.gz", "TC": F/"sumstats"/"TC.sumstats.gz"}
OWN = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
       "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal", "BMD": "stromal", "NEUT": "myeloid", "TC": "liver"}

def enriched():
    """농축 계통: 계통 지도 τ*의 단측 P < 0.05. 지도 형질은 3c M1, 검증 형질은 3f val M1 (마지막 행 = 계통 주석)."""
    m = pd.read_csv(R0/"e02_lineage_map"/"results"/"z04_sldsc.csv"); m = m[m.model == "M1"]
    E = {t: sorted(g.lineage[g.tau_p1 < .05]) for t, g in m.groupby("trait")}
    for t in ["BMD", "NEUT", "TC"]:
        E[t] = sorted(f.split("__M1_")[1][:-8] for f in glob.glob(str(F/"sldsc"/f"{t}__M1_*.results"))
                      if norm.sf(pd.read_csv(f, sep="\t")["Coefficient_z-score"].iloc[-1]) < .05)
    return E
E = enriched()
if STAGE == "dry":
    fake = F/"test"/"AD_shuffledZ.sumstats.gz"; assert fake.exists()
    jobs = [("FAKE", "adipose")]; TRAITS = {"FAKE": fake}
else:
    jobs = [(t, g) for t in TRAITS for g in E[t] if g != OWN[t]]
    for t in TRAITS: P(f"  {t:7s} 농축 계통 {E[t]}  (자기 계통 {OWN[t]} {'포함' if OWN[t] in E[t] else '없음'})")
    assert len(jobs) == 33, len(jobs)   # 프로토콜 표와 같아야 한다
P(f"단계 {STAGE}: S-LDSC {len(jobs)}회, 순차" + (" — 시험(Z 섞음, 결과 아님)" if STAGE == "dry" else ""))
t0 = time.time()
for i, (t, g) in enumerate(jobs):
    out = SL/f"{t}__T3all_{g}"
    if not Path(f"{out}.results").exists():
        refs = [str(Y/"annot"/f"L5_{g}_{m}"/f"L5_{g}_{m}.") for m in ["CAGE", "DNASE", "H3K27ac"]]
        cmd = ["python", "-W", "ignore", "-c", "import sys; sys.argv=['ldsc']+sys.argv[1:]; from ldscore.ldsc import main; main()",
               "--h2", str(TRAITS[t]), "--ref-ld-chr", ",".join([str(B/"baselineLD_v2.2"/"baselineLD.")] + refs), "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."),
               "--overlap-annot", "--frqfile-chr", str(B/"plink_files"/"1000G.EUR.hg38."), "--print-coefficients", "--print-delete-vals", "--out", str(out)]
        subprocess.run(cmd, cwd=BV, env=env, capture_output=True, text=True)
        if not Path(f"{out}.results").exists(): raise SystemExit(f"{t} {g} 실패: {open(str(out) + '.log', encoding='utf-8').read()[-1500:]}")
    el = (time.time() - t0) / 60; P(f"  [{i + 1}/{len(jobs)}] {t} {g} 완료 ({el:.0f}분, 남은 약 {el / (i + 1) * (len(jobs) - i - 1):.0f}분)")
(EXP/"results"/f"h01_{STAGE}.txt").write_text("\n".join(log), encoding="utf-8")
