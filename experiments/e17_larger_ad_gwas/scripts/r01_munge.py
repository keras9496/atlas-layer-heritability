"""실험 3r r01: 2026 AD 합의 메타분석 proxy 제외판(GCST90704647) 요약통계 정리 (PROTOCOL_exp3r_2026-10-25, 2026-10-25 고정).
 원래 AD 처리(code/00_prepare_sumstats.py)와 같은 규칙: SNV만, APOE 영역(chr19:43.9–45.9 Mb, GRCh38) 제외, HapMap3 no-MHC rsID만, rsID 중복은 첫 행.
 Z = beta / standard_error, P = p_value, FRQ = effect_allele_frequency, N = Neff_total(SNP별 유효 표본 수; 3n 환자–대조 형질과 같은 규칙).
 → munge_sumstats.py --signed-sumstats Z,0 (3n n01과 같다). 확인: HapMap3 no-MHC 목록의 80 % 이상이 남아야 한다(아니면 멈춘다).
출력: exp3r/sumstats/N_AD.sumstats.gz (+ .log), results/r01_munge.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import os, subprocess
import numpy as np, pandas as pd
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); R = Path(f"{DATA}/exp3r"); OUT = R/"sumstats"; OUT.mkdir(exist_ok=True)
LDSC = f"{DATA}/tools/ldsc"; env = dict(os.environ, PYTHONIOENCODING="utf-8")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
HM3 = set(pd.read_csv(X/"hm3_no_MHC.list.txt", header=None)[0])
COLS = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "p_value", "rs_id", "Neff_total"]
parts, n0, napoe = [], 0, 0
for ch in pd.read_csv(R/"raw"/"GCST90704647.tsv.gz", sep="\t", usecols=COLS, dtype={"chromosome": str}, chunksize=3_000_000, low_memory=False):
    n0 += len(ch)
    apoe = (ch.chromosome == "19") & (ch.base_pair_location >= 43_900_000) & (ch.base_pair_location <= 45_900_000); napoe += int(apoe.sum())
    ch = ch[~apoe & ch.rs_id.isin(HM3)]
    parts.append(ch[(ch.effect_allele.astype(str).str.len() == 1) & (ch.other_allele.astype(str).str.len() == 1)])
M = pd.concat(parts, ignore_index=True)
for c in ["beta", "standard_error", "p_value", "Neff_total", "effect_allele_frequency"]: M[c] = pd.to_numeric(M[c], errors="coerce")
M = M[(M.standard_error > 0) & M.beta.notna() & M.Neff_total.notna() & (M.Neff_total > 0)].drop_duplicates("rs_id")
M["Z"] = M.beta / M.standard_error
out = M.rename(columns={"rs_id": "SNP", "effect_allele": "A1", "other_allele": "A2", "p_value": "P", "Neff_total": "N", "effect_allele_frequency": "FRQ"})
out["A1"] = out.A1.str.upper(); out["A2"] = out.A2.str.upper(); out = out[["SNP", "A1", "A2", "Z", "P", "N", "FRQ"]]
cov = len(out) / len(HM3); p0 = int((out.P <= 0).sum())
P(f"[N_AD] 2026 AD 합의 메타분석 proxy 제외판 GCST90704647: 전체 {n0:,}행, APOE 영역 {napoe:,}행 제외 → HapMap3 no-MHC SNV {len(out):,}개 ({cov * 100:.0f} %),"
  f" Neff_total 중앙값 {out.N.median():,.0f} (범위 {out.N.min():,.0f}–{out.N.max():,.0f}) | 파일 P가 0인 SNP {p0:,}개")
assert cov >= .8, f"HapMap3 포함률 {cov:.2f} < 0.8 — 멈춘다"
tmp = OUT/"N_AD_in_hm3.tsv"; out.to_csv(tmp, sep="\t", index=False)
subprocess.run(["python", "-W", "ignore", "munge_sumstats.py", "--sumstats", str(tmp), "--signed-sumstats", "Z,0", "--out", str(OUT/"N_AD")], cwd=LDSC, check=True, env=env, capture_output=True)
t = open(OUT/"N_AD.log", encoding="utf-8").read(); P(f"  munge N_AD: {t[t.find('Writing'):].splitlines()[0]}")
P("  " + " | ".join(l.split(' - INFO - ')[-1].strip() for l in t.splitlines() if "Mean chi^2" in l or "Lambda GC" in l or "Max chi^2" in l))
(EXP/"results"/"r01_munge.txt").write_text("\n".join(log), encoding="utf-8")
