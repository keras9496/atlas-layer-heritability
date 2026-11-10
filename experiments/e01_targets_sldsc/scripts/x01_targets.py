"""실험 3 x01: S-LDSC hg38 참조(1000G EUR) 의 MAF ≥ 5 % SNP → REF/ALT 결정(1000G 고심도 hg38 pvar 와 위치·대립유전자 매칭) → Atlas 조회 목록.
출력: <DATA>/exp3/targets.parquet (chr, pos, rsid, ref, alt, atlas_variant), results/x01_targets.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess
import numpy as np, pandas as pd
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; D = Path(DATA); X = D/"exp3"; PF = X/"GRCh38"/"plink_files"; REF = D/"ref_1000G"
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
rows = []
for c in range(1, 23):
    b = pd.read_csv(PF/f"1000G.EUR.hg38.{c}.bim", sep="\t", header=None, names=["chr", "rsid", "cm", "pos", "a1", "a2"], dtype={"chr": str})
    f = pd.read_csv(PF/f"1000G.EUR.hg38.{c}.frq", sep=r"\s+")
    assert (f.SNP.to_numpy() == b.rsid.to_numpy()).all()
    b["maf"] = f.MAF.to_numpy(); rows.append(b[b.maf >= .05])
T = pd.concat(rows, ignore_index=True); P(f"참조 MAF≥5 % SNP {len(T):,}")
# 1000G 고심도 hg38 pvar (전체 조상, 이중대립 SNV) 에서 REF/ALT
pos = T[["chr", "pos"]].assign(e=T.pos, n="x"); pos.to_csv(X/"target_pos.txt", sep="\t", header=False, index=False)
if not (X/"hc.pvar").exists():
    subprocess.run([str(D/"tools"/"plink2.exe"), "--pgen", str(REF/"all_hg38.pgen"), "--pvar", str(REF/"all_hg38_noannot.pvar.zst"), "--psam", str(REF/"hg38_corrected.psam"),
                    "--extract", "range", str(X/"target_pos.txt"), "--snps-only", "just-acgt", "--make-just-pvar", "cols=", "--out", str(X/"hc"), "--threads", "16", "--silent"], check=True)
H = pd.read_csv(X/"hc.pvar", sep="\t", comment="#", header=None, names=["chr", "pos", "id", "ref", "alt"], usecols=range(5), dtype={"chr": str})
H = H[H.alt.str.len() == 1]
m = T.merge(H[["chr", "pos", "ref", "alt"]], on=["chr", "pos"], how="left")
ok = ((m.a1 == m.ref) & (m.a2 == m.alt)) | ((m.a1 == m.alt) & (m.a2 == m.ref))
m = m[ok].drop_duplicates(["chr", "pos", "rsid"], keep=False)                    # 같은 위치 다중 ALT 중 대립유전자가 맞는 한 행만
P(f"REF/ALT 결정 {m.rsid.nunique():,} ({m.rsid.nunique() / len(T):.1%}); 미결정은 Atlas 주석 0 으로 둔다")
m["atlas_variant"] = "chr" + m.chr + ":" + m.pos.astype(str) + ":" + m.ref + ">" + m.alt
m[["chr", "pos", "rsid", "ref", "alt", "maf", "atlas_variant"]].reset_index(drop=True).to_parquet(X/"targets.parquet", index=False)
(EXP/"results"/"x01_targets.txt").write_text("\n".join(log), encoding="utf-8")
