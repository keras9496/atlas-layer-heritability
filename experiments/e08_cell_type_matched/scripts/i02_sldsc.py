"""실험 3i i02: 세포형 일치 S-LDSC (PROTOCOL_exp3i_2026-10-13). 순차 실행. 호출은 3f v02·3h h01과 같다.
 M3_{g} = baselineLD + M5_{g}_{CAGE,DNASE,H3K27ac} (세포형 일치, 3-층위 동시). 형질 14개 × 자기 계통 = 14회.
 단계: dry = Z를 섞은 가짜 형질로 골수계 1회(형식 확인, 결과 아님) / main = 14회
 재실행 안전: .results 가 있으면 건너뛴다.  사용: python i02_sldsc.py dry | main
출력: exp3i/sldsc[/test]/<형질>__M3_<계통>.{results,log,part_delete}, results/i02_<단계>.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess, os, sys, time
from pathlib import Path
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); F = Path(f"{DATA}/exp3f"); I = Path(f"{DATA}/exp3i")
BV = Path(f"{DATA}/tools/ldsc_bv"); B = X/"GRCh38"; SS = Path(f"{DATA}/exp3b/sumstats")
SL = I/"test"/"sldsc" if STAGE == "dry" else I/"sldsc"; SL.mkdir(parents=True, exist_ok=True); (EXP/"results").mkdir(exist_ok=True)
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
TRAITS = {"AD": X/"bellenguez_noAPOE.sumstats.gz", "MONO": Path(f"{DATA}/exp3b/monocyte_count.sumstats.gz"), "EA": SS/"PASS_Years_of_Education2.sumstats",
          "SCZ": SS/"PASS_Schizophrenia.sumstats", "BIP": SS/"PASS_Bipolar_Disorder.sumstats", "NEURO": SS/"PASS_Neuroticism.sumstats", "BMI": SS/"PASS_BMI1.sumstats",
          "IBD": SS/"PASS_IBD.sumstats", "RA": SS/"PASS_Rheumatoid_Arthritis.sumstats", "LDL": SS/"PASS_LDL.sumstats", "HEIGHT": SS/"PASS_Height1.sumstats",
          "BMD": F/"sumstats"/"BMD.sumstats.gz", "NEUT": F/"sumstats"/"NEUT.sumstats.gz", "TC": F/"sumstats"/"TC.sumstats.gz"}
OWN = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
       "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal", "BMD": "stromal", "NEUT": "myeloid", "TC": "liver"}
AN = I/"annot"   # 시험(dry)도 본 주석을 쓴다 (S-LDSC는 22개 염색체가 모두 필요하다; 3g와 같다)
if STAGE == "dry": TRAITS = {"FAKE": F/"test"/"AD_shuffledZ.sumstats.gz"}; OWN = {"FAKE": "myeloid"}
jobs = list(TRAITS)
P(f"단계 {STAGE}: S-LDSC {len(jobs)}회, 순차" + (" — 시험(Z 섞음, 결과 아님)" if STAGE == "dry" else ""))
t0 = time.time()
for i, t in enumerate(jobs):
    g = OWN[t]; out = SL/f"{t}__M3_{g}"
    if not Path(f"{out}.results").exists():
        cmd = ["python", "-W", "ignore", "-c", "import sys; sys.argv=['ldsc']+sys.argv[1:]; from ldscore.ldsc import main; main()",
               "--h2", str(TRAITS[t]), "--ref-ld-chr", ",".join([str(B/"baselineLD_v2.2"/"baselineLD."), str(AN/f"M3_{g}"/f"M3_{g}.")]), "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."),
               "--overlap-annot", "--frqfile-chr", str(B/"plink_files"/"1000G.EUR.hg38."), "--print-coefficients", "--print-delete-vals", "--out", str(out)]
        subprocess.run(cmd, cwd=BV, env=env, capture_output=True, text=True)
        if not Path(f"{out}.results").exists(): raise SystemExit(f"{t} {g} 실패: {open(str(out) + '.log', encoding='utf-8').read()[-1500:]}")
    el = (time.time() - t0) / 60; P(f"  [{i + 1}/{len(jobs)}] {t} {g} 완료 ({el:.0f}분, 남은 약 {el / (i + 1) * (len(jobs) - i - 1):.0f}분)")
(EXP/"results"/f"i02_{STAGE}.txt").write_text("\n".join(log), encoding="utf-8")
