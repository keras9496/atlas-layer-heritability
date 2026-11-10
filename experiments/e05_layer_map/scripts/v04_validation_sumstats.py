"""실험 3f v04: 검증 형질 요약통계 준비 (PROTOCOL_exp3f_2026-10-09 3단계). 실험 3b x06의 단핵구 수와 같은 절차.
 GWAS Catalog 조화 파일 → rsID·효과/비효과 대립유전자·beta/se → Z, SNV만, HapMap3 no-MHC → munge.
 형식 두 가지: 옛 조화 형식(hm_rsid, hm_effect_allele, hm_beta 등) / 새 형식(rsid, effect_allele, beta 등). se는 원 열(standard_error)을 쓴다.
 N: 파일에 표본 수 열이 믿을 만하지 않으면 연구 보고 N(상수)을 쓴다. 골밀도 파일의 'n' 열은 실제로는 다른 P값(P.I)이라 쓰지 않는다.
출력: exp3f/sumstats/{BMD,NEUT,TC}.sumstats.gz, results/v04_validation_sumstats.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import os, subprocess, gzip
import numpy as np, pandas as pd
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); F = Path(f"{DATA}/exp3f"); RAW = F/"sumstats"/"raw"
LDSC = f"{DATA}/tools/ldsc"; env = dict(os.environ, PYTHONIOENCODING="utf-8")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
HM3 = set(pd.read_csv(X/"hm3_no_MHC.list.txt", header=None)[0])
SRC = {"BMD": ("30598549-GCST006979-EFO_0009270.h.tsv.gz", 426824, "Morris 2019 eBMD, GCST006979"),
       "NEUT": ("32888494-GCST90002398-EFO_0004833.h.tsv.gz", 408112, "Vuckovic 2020 호중구 수, GCST90002398"),
       "TC": ("GCST90500025.h.tsv.gz", 434646, "UKB 총콜레스테롤, GCST90500025")}
def munge(inp, out):
    subprocess.run(["python", "-W", "ignore", "munge_sumstats.py", "--sumstats", str(inp), "--signed-sumstats", "Z,0", "--out", str(out)], cwd=LDSC, check=True, env=env, capture_output=True)
    t = open(f"{out}.log", encoding="utf-8").read(); P(f"  munge {Path(out).name}: {t[t.find('Writing'):].splitlines()[0]}"); P("  " + " | ".join(l.split(' - INFO - ')[-1] for l in t.splitlines() if "Mean chi^2" in l or "Max chi^2" in l))
for k, (fn, N, lab) in SRC.items():
    f = RAW/fn
    with gzip.open(f, "rt") as fh: hdr = fh.readline().rstrip("\n").split("\t")
    if "hm_rsid" in hdr: c = dict(rs="hm_rsid", ea="hm_effect_allele", oa="hm_other_allele", b="hm_beta", se="standard_error", p="p_value", eaf="hm_effect_allele_frequency")
    else: c = dict(rs="rsid", ea="effect_allele", oa="other_allele", b="beta", se="standard_error", p="p_value", eaf="effect_allele_frequency")
    if c["rs"] not in hdr and "variant_id" in hdr: c["rs"] = "variant_id"
    M = pd.read_csv(f, sep="\t", usecols=list(c.values()), low_memory=False).rename(columns={v: kk for kk, v in c.items()})
    n0 = len(M); M = M[M.rs.isin(HM3) & (M.ea.astype(str).str.len() == 1) & (M.oa.astype(str).str.len() == 1)].copy()
    M["b"] = pd.to_numeric(M.b, errors="coerce"); M["se"] = pd.to_numeric(M.se, errors="coerce"); M = M[(M.se > 0) & M.b.notna()].drop_duplicates("rs")
    M["Z"] = M.b / M.se; M["N"] = N
    out = M.rename(columns={"rs": "SNP", "ea": "A1", "oa": "A2", "p": "P", "eaf": "FRQ"})[["SNP", "A1", "A2", "Z", "P", "N", "FRQ"]]
    P(f"[{k}] {lab}: 전체 {n0:,}행 → HapMap3 no-MHC SNV {len(out):,}개 (HapMap3 목록 {len(HM3):,}개의 {len(out) / len(HM3) * 100:.0f} %), N {N:,}")
    tmp = F/"sumstats"/f"{k}_in_hm3.tsv"; out.to_csv(tmp, sep="\t", index=False); munge(tmp, F/"sumstats"/k)
(EXP/"results"/"v04_validation_sumstats.txt").write_text("\n".join(log), encoding="utf-8")
