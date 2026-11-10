"""실험 3j j01: 형질 공통 jackknife 블록 경계 (PROTOCOL_exp3j_2026-10-14).
 기준 SNP = baselineLD v2.2 LD 점수 SNP ∩ 회귀 가중치(HapMap3, MHC 제외) SNP, (염색체, 위치) 순서 (ldsc가 회귀 SNP를 정렬하는 순서와 같다).
 ldsc와 같은 식 np.floor(np.linspace(0, N, 201))로 200등분하고 각 블록 첫 SNP의 위치를 경계로 고정한다. 형질마다 회귀 SNP를 이 경계로 나눈다(j02).
출력: exp3j/blocks.npz (key = 염색체 × 10¹⁰ + 위치; 경계 199개), exp3j/snp_key.parquet (SNP → key), results/j01_blocks.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import numpy as np, pandas as pd
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; B = Path(f"{DATA}/exp3/GRCh38"); J = Path(f"{DATA}/exp3j")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
rd = lambda f: pd.read_csv(f, sep="\t", usecols=["CHR", "SNP", "BP"])
bl = pd.concat([rd(B/"baselineLD_v2.2"/f"baselineLD.{c}.l2.ldscore.gz") for c in range(1, 23)]).drop_duplicates("SNP")
w = pd.concat([rd(B/"weights"/f"weights.hm3_noMHC.{c}.l2.ldscore.gz") for c in range(1, 23)]).drop_duplicates("SNP")
P(f"baselineLD SNP {len(bl):,}, 가중치 SNP {len(w):,}")
m = bl.merge(w[["SNP", "CHR", "BP"]], on="SNP", suffixes=("", "_w"))
P(f"공통 {len(m):,}; 위치가 다른 SNP {int(((m.CHR != m.CHR_w) | (m.BP != m.BP_w)).sum())}개 (경계와 형질 SNP 위치는 baselineLD 위치를 쓴다)")
m["key"] = m.CHR.astype(np.int64) * 10**10 + m.BP.astype(np.int64); m = m.sort_values("key").reset_index(drop=True)
N = len(m); sep = np.floor(np.linspace(0, N, 201)).astype(int); bounds = m.key.to_numpy()[sep[1:-1]]
assert (np.diff(bounds) > 0).all()
np.savez(J/"blocks.npz", bounds=bounds, n_ref=N)
bl.assign(key=bl.CHR.astype(np.int64) * 10**10 + bl.BP.astype(np.int64))[["SNP", "key"]].to_parquet(J/"snp_key.parquet", index=False)
st = m.groupby(np.searchsorted(bounds, m.key.to_numpy(), side="right")).agg(n=("SNP", "size"), chr0=("CHR", "first"), chr1=("CHR", "last"), bp0=("BP", "first"), bp1=("BP", "last"))
span = np.where(st.chr0 == st.chr1, (st.bp1 - st.bp0) / 1e6, np.nan)
P(f"블록 200개: 기준 SNP {st.n.min():,}–{st.n.max():,}개, 길이 중앙값 {np.nanmedian(span):.1f} Mb (염색체를 넘는 블록 {int(np.isnan(span).sum())}개)")
(EXP/"results"/"j01_blocks.txt").write_text("\n".join(log), encoding="utf-8")
