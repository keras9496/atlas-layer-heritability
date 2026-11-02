"""실험 3g u02: S-LDSC (PROTOCOL_exp3g_2026-10-11). 순차 실행. 호출은 실험 3f v02와 같다.
 G3_{q} = baselineLD + G3_{q}_{g} (CAGE, DNASE, H3K27ac 동시, 전체 트랙, 상위 q), q = L0p1 · L0p5 · L1p0.  5 % 참조 = 3f T3all (새로 돌리지 않음).
 단계: dry = 3f의 Z를 섞은 가짜 형질, 골수계, 3회 (결과 아님, 형식 확인) / main = 본 형질 11개 + 검증 형질 3개 × 3 = 42회
 재실행 안전: .results 가 있으면 건너뛴다.  사용: python u02_sldsc.py dry | main
출력: exp3g/sldsc[/test]/<형질>__G3_<q>.{results,log,part_delete}, results/u02_<단계>.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess, os, sys, time
from pathlib import Path
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); F = Path(f"{DATA}/exp3f"); G = Path(f"{DATA}/exp3g")
BV = Path(f"{DATA}/tools/ldsc_bv"); B = X/"GRCh38"; SS = Path(f"{DATA}/exp3b/sumstats")
SL = G/"test"/"sldsc" if STAGE == "dry" else G/"sldsc"; SL.mkdir(parents=True, exist_ok=True)
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
TRAITS = {"AD": X/"bellenguez_noAPOE.sumstats.gz", "MONO": Path(f"{DATA}/exp3b/monocyte_count.sumstats.gz"), "EA": SS/"PASS_Years_of_Education2.sumstats",
          "SCZ": SS/"PASS_Schizophrenia.sumstats", "BIP": SS/"PASS_Bipolar_Disorder.sumstats", "NEURO": SS/"PASS_Neuroticism.sumstats", "BMI": SS/"PASS_BMI1.sumstats",
          "IBD": SS/"PASS_IBD.sumstats", "RA": SS/"PASS_Rheumatoid_Arthritis.sumstats", "LDL": SS/"PASS_LDL.sumstats", "HEIGHT": SS/"PASS_Height1.sumstats",
          "BMD": F/"sumstats"/"BMD.sumstats.gz", "NEUT": F/"sumstats"/"NEUT.sumstats.gz", "TC": F/"sumstats"/"TC.sumstats.gz"}
LINE = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
        "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal", "BMD": "stromal", "NEUT": "myeloid", "TC": "liver"}   # 3f와 같다
QS = ["L0p1", "L0p5", "L1p0"]
if STAGE == "dry": TRAITS, LINE = {"FAKE": F/"test"/"AD_shuffledZ.sumstats.gz"}, {"FAKE": "myeloid"}
jobs = [(t, q) for t in TRAITS for q in QS]
P(f"단계 {STAGE}: S-LDSC {len(jobs)}회, 순차" + (" — 시험(Z 섞음, 결과 아님)" if STAGE == "dry" else ""))
t0 = time.time()
for i, (t, q) in enumerate(jobs):
    out = SL/f"{t}__G3_{q}"; g = LINE[t]
    if not Path(f"{out}.results").exists():
        refs = [str(B/"baselineLD_v2.2"/"baselineLD."), str(G/"annot"/f"G3_{q}_{g}"/f"G3_{q}_{g}.")]
        cmd = ["python", "-W", "ignore", "-c", "import sys; sys.argv=['ldsc']+sys.argv[1:]; from ldscore.ldsc import main; main()",
               "--h2", str(TRAITS[t]), "--ref-ld-chr", ",".join(refs), "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."),
               "--overlap-annot", "--frqfile-chr", str(B/"plink_files"/"1000G.EUR.hg38."), "--print-coefficients", "--print-delete-vals", "--out", str(out)]
        subprocess.run(cmd, cwd=BV, env=env, capture_output=True, text=True)
        if not Path(f"{out}.results").exists(): raise SystemExit(f"{t} {q} 실패: {open(str(out) + '.log', encoding='utf-8').read()[-1500:]}")
    el = (time.time() - t0) / 60; P(f"  [{i + 1}/{len(jobs)}] {t} G3_{q} 완료 ({el:.0f}분, 남은 약 {el / (i + 1) * (len(jobs) - i - 1):.0f}분)")
(EXP/"results"/f"u02_{STAGE}.txt").write_text("\n".join(log), encoding="utf-8")
