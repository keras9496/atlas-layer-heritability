"""실험 3 x03: Atlas 값 → 주석(PROTOCOL_exp3) → 염색체별 LD 점수(ldsc --l2, 1 cM, 1000G EUR hg38) → 주석별 파일로 분리.
 주석 51개: ct5_/ct1_/ct10_<c> (세포형 S_c 상위 5/1/10 %), ctc_<c> (S_c 백분위 연속), mod_<c>_<m> (모달리티 상위 5 %), avi5.
 MAF < 5 % 또는 REF/ALT 미결정 SNP = 0.
LD 점수: abdenlab/ldsc-python3 (pandas·Windows 패치; chr22 에서 공식 baselineLD base L2 와 완전 일치 확인)
출력: <DATA>/exp3/annot/<이름>/<이름>.<chr>.{annot.gz,l2.ldscore.gz,l2.M,l2.M_5_50}, annot/ct5set/ (세포형 5개 동시), results/x03_annot.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess, os
import numpy as np, pandas as pd
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); PF = X/"GRCh38"/"plink_files"; TEST = os.environ.get("EXP3_TEST") == "1"; AN = X/"test"/"annot" if TEST else X/"annot"; AN.mkdir(exist_ok=True, parents=True); CHRS = [int(os.environ.get("EXP3_TEST_CHR", 22))] if TEST else range(1, 23)   # 시험: chr22 만, 별도 폴더
LDSC = Path(f"{DATA}/tools/ldsc/ldsc.py")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
CTX = ["myeloid", "neuron", "astrocyte", "liver", "nonliver"]; MODS = ["DNASE", "H3K27ac", "H3K4me1", "H3K4me3", "CAGE", "RNA"]; CORE = ["DNASE", "H3K27ac", "H3K4me1", "CAGE"]
T = pd.read_parquet(X/"targets.parquet"); A = pd.concat([pd.read_parquet(f) for f in sorted((X/"atlas").glob("part_*.parquet"))], ignore_index=True).drop_duplicates("atlas_variant")
T = T.merge(A, on="atlas_variant", how="left"); P(f"대상 {len(T):,}, Atlas 값 있음 {int(T['DNASE_myeloid'].notna().sum()):,}")
ann = pd.DataFrame({"rsid": T.rsid})
for c in CTX:
    S = pd.concat([T[f"{m}_{c}"].astype(float).rank(pct=True) for m in CORE], axis=1).mean(1); Sp = S.rank(pct=True)
    for th, nm in [(.95, "ct5"), (.99, "ct1"), (.90, "ct10")]: ann[f"{nm}_{c}"] = (Sp >= th).astype(np.int8)
    ann[f"ctc_{c}"] = Sp.fillna(0).astype(np.float32)
    for m in MODS: ann[f"mod_{c}_{m}"] = (T[f"{m}_{c}"].astype(float).rank(pct=True) >= .95).astype(np.int8)
ann["avi5"] = (T.avi.astype(float).rank(pct=True) >= .95).astype(np.int8)
NAMES = [c for c in ann.columns if c != "rsid"]; P(f"주석 {len(NAMES)}개; 상위 5 % 세포형 주석 간 겹침(Jaccard):")
J = pd.DataFrame([[((ann[f"ct5_{a}"] & ann[f"ct5_{b}"]).sum() / max((ann[f"ct5_{a}"] | ann[f"ct5_{b}"]).sum(), 1)) for b in CTX] for a in CTX], index=CTX, columns=CTX).round(2); P(J.to_string())
ann = ann.drop_duplicates("rsid").set_index("rsid")
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
for c in CHRS:
    out_all = AN/"_all"/f"all.{c}"; (AN/"_all").mkdir(exist_ok=True)
    if not Path(f"{out_all}.l2.ldscore.gz").exists():
        b = pd.read_csv(PF/f"1000G.EUR.hg38.{c}.bim", sep="\t", header=None, names=["CHR", "SNP", "CM", "BP", "a1", "a2"])
        full = ann.reindex(b.SNP).fillna(0).reset_index(drop=True); full[[n for n in NAMES if not n.startswith("ctc_")]] = full[[n for n in NAMES if not n.startswith("ctc_")]].astype(np.int8)
        full.to_csv(f"{out_all}.annot.gz", sep="\t", index=False, compression="gzip")
        r = subprocess.run(["python", "-W", "ignore", str(LDSC), "--l2", "--bfile", str(PF/f"1000G.EUR.hg38.{c}"), "--ld-wind-cm", "1", "--annot", f"{out_all}.annot.gz", "--thin-annot",
                            "--out", str(out_all), "--print-snps", str(X/"hm3_no_MHC.list.txt")], cwd=LDSC.parent, env=env, capture_output=True, text=True)
        if not Path(f"{out_all}.l2.ldscore.gz").exists(): raise SystemExit(f"chr{c} l2 실패\n{r.stdout[-1500:]}")
    # 주석별 분리
    L = pd.read_csv(f"{out_all}.l2.ldscore.gz", sep="\t"); Aall = pd.read_csv(f"{out_all}.annot.gz", sep="\t")
    M = open(f"{out_all}.l2.M").read().split(); M5 = open(f"{out_all}.l2.M_5_50").read().split()
    groups = {n: [n] for n in NAMES}; groups["ct5set"] = [f"ct5_{x}" for x in CTX]
    for g, cols in groups.items():
        d = AN/g; d.mkdir(exist_ok=True); idx = [NAMES.index(x) for x in cols]
        L[["CHR", "SNP", "BP"] + [f"{x}L2" for x in cols]].to_csv(d/f"{g}.{c}.l2.ldscore.gz", sep="\t", index=False, compression="gzip")
        Aall[cols].to_csv(d/f"{g}.{c}.annot.gz", sep="\t", index=False, compression="gzip")
        open(d/f"{g}.{c}.l2.M", "w").write("\t".join(M[i] for i in idx)); open(d/f"{g}.{c}.l2.M_5_50", "w").write("\t".join(M5[i] for i in idx))
    P(f"chr{c} 완료")
# τ* 계산용: 흔한 SNP 에서 주석 SD
sd = ann[NAMES].std(ddof=0); sd.rename("sd_common").to_csv(AN.parent/"annot_sd_common.csv" if TEST else X/"annot_sd_common.csv")
if not TEST: (EXP/"results"/"x03_annot.txt").write_text("\n".join(log), encoding="utf-8")
