"""실험 3f v02: S-LDSC (PROTOCOL_exp3f_2026-10-09). 순차 실행. 호출은 실험 3d y06·3e w03과 같다.
 T3_{s} = baselineLD + T3_{s}_{g} (CAGE, DNASE, H3K27ac 재추출 주석 동시), s = 01..20.   T3all = baselineLD + 3d L5_{g}_{CAGE,DNASE,H3K27ac} (전체 트랙, 참조)
 M1_{g} (검증 단계만) = baselineLD + 3c L5_{g} (계통 지도, 12 계통).
 단계: dry = Z를 섞은 가짜 형질로 T3 시드 01·02 + T3all (결과 아님, 형식 확인) / main = 형질 11개 × (T3 20 + T3all) = 231회 / val = 검증 형질 3개 × (M1 12 + T3 20 + T3all) = 99회
 재실행 안전: .results 가 있으면 건너뛴다.  사용: python v02_sldsc.py dry | main | val
출력: exp3f/sldsc[/test]/<형질>__<모형>.{results,log,part_delete}, results/v02_<단계>.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess, os, sys, time
import numpy as np, pandas as pd
from pathlib import Path
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); C3 = Path(f"{DATA}/exp3c"); Y = Path(f"{DATA}/exp3d"); F = Path(f"{DATA}/exp3f")
BV = Path(f"{DATA}/tools/ldsc_bv"); B = X/"GRCh38"; SS = Path(f"{DATA}/exp3b/sumstats")
SL = F/"test"/"sldsc" if STAGE == "dry" else F/"sldsc"; SL.mkdir(parents=True, exist_ok=True)
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
TRAITS = {"AD": X/"bellenguez_noAPOE.sumstats.gz", "MONO": Path(f"{DATA}/exp3b/monocyte_count.sumstats.gz"), "EA": SS/"PASS_Years_of_Education2.sumstats",
          "SCZ": SS/"PASS_Schizophrenia.sumstats", "BIP": SS/"PASS_Bipolar_Disorder.sumstats", "NEURO": SS/"PASS_Neuroticism.sumstats", "BMI": SS/"PASS_BMI1.sumstats",
          "IBD": SS/"PASS_IBD.sumstats", "RA": SS/"PASS_Rheumatoid_Arthritis.sumstats", "LDL": SS/"PASS_LDL.sumstats", "HEIGHT": SS/"PASS_Height1.sumstats"}
LINE = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
        "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal"}   # 3d-A·3e와 같다 (프로토콜 고정)
VTRAITS = {"BMD": F/"sumstats"/"BMD.sumstats.gz", "NEUT": F/"sumstats"/"NEUT.sumstats.gz", "TC": F/"sumstats"/"TC.sumstats.gz"}
VLINE = {"BMD": "stromal", "NEUT": "myeloid", "TC": "liver"}                          # 프로토콜 3단계 표 (예측 계통)
LIN12 = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
SEEDS = [f"{s:02d}" for s in range(1, 21)]
if STAGE == "dry":
    fake = F/"test"/"AD_shuffledZ.sumstats.gz"
    if not fake.exists():
        d = pd.read_csv(TRAITS["AD"], sep="\t"); d["Z"] = np.random.default_rng(20261009).permutation(d.Z.to_numpy()); d.to_csv(fake, sep="\t", index=False, compression="gzip")
    TRAITS, LINE = {"FAKE": fake}, {"FAKE": "myeloid"}; SEEDS = SEEDS[:2]
elif STAGE == "val":
    TRAITS, LINE = VTRAITS, VLINE
def refs(model, g):
    if model == "T3all": return [str(Y/"annot"/f"L5_{g}_{m}"/f"L5_{g}_{m}.") for m in ["CAGE", "DNASE", "H3K27ac"]]
    if model.startswith("M1_"): h = model[3:]; return [str(C3/"annot"/f"L5_{h}"/f"L5_{h}.")]
    s = model.split("_")[1]; return [str(F/"annot"/f"T3_{s}_{g}"/f"T3_{s}_{g}.")]
jobs = []
for t in TRAITS:
    g = LINE[t]
    if STAGE == "val": jobs += [(t, f"M1_{h}") for h in LIN12]
    jobs += [(t, f"T3_{s}") for s in SEEDS] + [(t, "T3all")]
P(f"단계 {STAGE}: S-LDSC {len(jobs)}회, 순차" + (" — 시험(Z 섞음, 결과 아님)" if STAGE == "dry" else ""))
t0 = time.time()
for i, (t, model) in enumerate(jobs):
    out = SL/f"{t}__{model}"
    if not Path(f"{out}.results").exists():
        cmd = ["python", "-W", "ignore", "-c", "import sys; sys.argv=['ldsc']+sys.argv[1:]; from ldscore.ldsc import main; main()",
               "--h2", str(TRAITS[t]), "--ref-ld-chr", ",".join([str(B/"baselineLD_v2.2"/"baselineLD.")] + refs(model, LINE[t])), "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."),
               "--overlap-annot", "--frqfile-chr", str(B/"plink_files"/"1000G.EUR.hg38."), "--print-coefficients", "--print-delete-vals", "--out", str(out)]
        subprocess.run(cmd, cwd=BV, env=env, capture_output=True, text=True)
        if not Path(f"{out}.results").exists(): raise SystemExit(f"{t} {model} 실패: {open(str(out) + '.log', encoding='utf-8').read()[-1500:]}")
    el = (time.time() - t0) / 60; P(f"  [{i + 1}/{len(jobs)}] {t} {model} 완료 ({el:.0f}분, 남은 약 {el / (i + 1) * (len(jobs) - i - 1):.0f}분)")
(EXP/"results"/f"v02_{STAGE}.txt").write_text("\n".join(log), encoding="utf-8")
