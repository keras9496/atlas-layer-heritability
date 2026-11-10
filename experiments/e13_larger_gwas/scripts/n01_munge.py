"""실험 3n n01: 최신 GWAS 요약통계 정리 (PROTOCOL_exp3n_2026-10-21, 2026-10-21 고정). 3f v04와 같은 절차.
 rsID·효과/비효과 대립유전자·효과 크기/SE → Z, SNV만, HapMap3 no-MHC → munge (--signed-sumstats Z,0).
 효과 크기: 키·BMI·LDL β, 조현병 BETA(log OR), 양극성장애 OR → log OR. LDL은 GC 보정 전 EFFECT_SIZE/SE.
 N: 키·BMI·LDL은 N 열, 조현병 NEFF, 양극성장애 2 × Neff_half. FRQ: 파일의 효과 대립유전자 빈도(조현병은 대조군 FCON).
 P: 파일의 P 열(LDL은 GC 보정 전 pvalue)을 munge에 넘긴다(3f v04와 같다; munge는 P로 χ²를, Z로 방향을 정한다). P가 0으로 기록된 SNP 수를 적는다.
 확인: HapMap3 no-MHC 목록의 80 % 이상이 남아야 한다(아니면 멈춘다).
출력: exp3n/sumstats/{HEIGHT,BMI,SCZ,BIP,LDL}.sumstats.gz, results/n01_munge.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import os, subprocess, gzip
import numpy as np, pandas as pd
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); N = Path(f"{DATA}/exp3n"); RAW = N/"raw"; OUT = N/"sumstats"; OUT.mkdir(exist_ok=True)
LDSC = f"{DATA}/tools/ldsc"; env = dict(os.environ, PYTHONIOENCODING="utf-8")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
HM3 = set(pd.read_csv(X/"hm3_no_MHC.list.txt", header=None)[0])
# 형질 → (파일, 구분자, 건너뛸 '##' 줄 여부, 열 대응, 효과 크기 종류, N 규칙, 설명)
SRC = {
 "HEIGHT": ("HEIGHT_Yengo2022_EUR.gz", "\t", dict(rs="RSID", ea="EFFECT_ALLELE", oa="OTHER_ALLELE", b="BETA", se="SE", p="P", n="N", frq="EFFECT_ALLELE_FREQ"), "beta", "Yengo 2022 키 EUR (GIANT, UKB 포함, 23andMe 제외)"),
 "BMI": ("BMI_Yengo2018_UPDATED.txt.gz", "\t", dict(rs="SNP", ea="Tested_Allele", oa="Other_Allele", b="BETA", se="SE", p="P", n="N", frq="Freq_Tested_Allele_in_HRS"), "beta", "Yengo 2018 BMI (GIANT + UKB)"),
 "SCZ": ("SCZ_PGC3_wave3_european.vcf.tsv.gz", "\t", dict(rs="ID", ea="A1", oa="A2", b="BETA", se="SE", p="PVAL", n="NEFF", frq="FCON"), "beta", "Trubetskoy 2022 조현병 PGC3 유럽인"),
 "BIP": ("BIP_bip2024_eur_no23andMe.gz", r"\s+", dict(rs="SNP", ea="A1", oa="A2", b="OR", se="SE", p="P", n="Neff_half", frq="HRC_FRQ_A1"), "or", "O'Connell 2025 양극성장애 유럽인 (23andMe 제외)"),
 "LDL": ("LDL_GLGC2021_EUR.gz", "\t", dict(rs="rsID", ea="ALT", oa="REF", b="EFFECT_SIZE", se="SE", p="pvalue", n="N", frq="POOLED_ALT_AF"), "beta", "Graham 2021 LDL 유럽인 (GLGC, UKB 포함; GC 보정 전 pvalue)"),
}
def skip_hash(f):
    with gzip.open(f, "rt", encoding="utf-8", errors="replace") as fh:
        k = 0
        for line in fh:
            if not line.startswith("##"): return k
            k += 1
def munge(inp, out):
    subprocess.run(["python", "-W", "ignore", "munge_sumstats.py", "--sumstats", str(inp), "--signed-sumstats", "Z,0", "--out", str(out)], cwd=LDSC, check=True, env=env, capture_output=True)
    t = open(f"{out}.log", encoding="utf-8").read(); P(f"  munge {Path(out).name}: {t[t.find('Writing'):].splitlines()[0]}")
    P("  " + " | ".join(l.split(' - INFO - ')[-1].strip() for l in t.splitlines() if "Mean chi^2" in l or "Lambda GC" in l or "Max chi^2" in l))
for k, (fn, sep, c, kind, lab) in SRC.items():
    f = RAW/fn; sk = skip_hash(f); parts, n0 = [], 0
    for ch in pd.read_csv(f, sep=sep, usecols=list(c.values()), skiprows=sk, chunksize=4_000_000, low_memory=False):
        n0 += len(ch); ch = ch.rename(columns={v: kk for kk, v in c.items()}); parts.append(ch[ch.rs.isin(HM3)])
    M = pd.concat(parts, ignore_index=True)
    M = M[(M.ea.astype(str).str.len() == 1) & (M.oa.astype(str).str.len() == 1)].copy()
    for col in ["b", "se", "p", "n", "frq"]: M[col] = pd.to_numeric(M[col], errors="coerce")
    if kind == "or": M = M[M.b > 0]; M["b"] = np.log(M.b)
    if k == "BIP": M["n"] = 2 * M.n
    M = M[(M.se > 0) & M.b.notna() & M.n.notna() & (M.n > 0)].drop_duplicates("rs")   # 3f v04와 같이 첫 행 유지
    M["Z"] = M.b / M.se; M["ea"] = M.ea.str.upper(); M["oa"] = M.oa.str.upper()
    out = M.rename(columns={"rs": "SNP", "ea": "A1", "oa": "A2", "p": "P", "n": "N", "frq": "FRQ"})[["SNP", "A1", "A2", "Z", "P", "N", "FRQ"]]
    cov = len(out) / len(HM3); p0 = int((out.P <= 0).sum())
    P(f"[{k}] {lab}: 전체 {n0:,}행 → HapMap3 no-MHC SNV {len(out):,}개 ({cov * 100:.0f} %), N 중앙값 {out.N.median():,.0f} (범위 {out.N.min():,.0f}–{out.N.max():,.0f})"
      f" | 파일 P가 0인 SNP {p0:,}개(munge가 범위 밖으로 뺀다; |Z| 최대 {out.Z.abs().max():.1f})")
    assert cov >= .8, f"{k}: HapMap3 포함률 {cov:.2f} < 0.8 — 멈춘다"
    tmp = OUT/f"{k}_in_hm3.tsv"; out.to_csv(tmp, sep="\t", index=False); munge(tmp, OUT/k)
(EXP/"results"/"n01_munge.txt").write_text("\n".join(log), encoding="utf-8")
