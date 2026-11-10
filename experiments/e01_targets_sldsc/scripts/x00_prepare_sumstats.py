"""실험 3 준비: Bellenguez stage I 요약통계 → S-LDSC 입력 (APOE 제외, HapMap3 no-MHC SNP) → munge_sumstats.
※ 2026-09-25에 대화형으로 실행한 명령을 그대로 스크립트로 옮긴 것이다(당시 별도 스크립트 파일은 없었다). 결과 파일의 로그
  (<DATA>/exp3/bellenguez_noAPOE.log)와 수치가 같다: 1,177,203 SNP, mean χ² 1.269.
입력: <DATA>/sumstats/GCST90027158_buildGRCh38.tsv.gz, exp3/hm3_no_MHC.list.txt (Zenodo 10515792)
출력: exp3/bellenguez_ldsc_in_hm3.tsv, exp3/bellenguez_noAPOE.sumstats.gz (+ .log)
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import os, subprocess
import pandas as pd

X = f"{DATA}/exp3"; LDSC = f"{DATA}/tools/ldsc"
S = pd.read_csv(f"{DATA}/sumstats/GCST90027158_buildGRCh38.tsv.gz", sep="\t",
                usecols=["variant_id", "p_value", "chromosome", "base_pair_location", "effect_allele", "other_allele",
                         "effect_allele_frequency", "beta", "standard_error", "n_cases", "n_controls"], dtype={"chromosome": str})
S = S[(S.effect_allele.str.len() == 1) & (S.other_allele.str.len() == 1)]
c, p = S.chromosome, S.base_pair_location
S = S[~((c == "19") & (p >= 43_900_000) & (p <= 45_900_000))]              # APOE 제외
S["N"] = S.n_cases + S.n_controls; S["Z"] = S.beta / S.standard_error
S = S[["variant_id", "effect_allele", "other_allele", "Z", "p_value", "N", "effect_allele_frequency"]].rename(
    columns={"variant_id": "SNP", "effect_allele": "A1", "other_allele": "A2", "p_value": "P", "effect_allele_frequency": "FRQ"})
h = set(pd.read_csv(f"{X}/hm3_no_MHC.list.txt", header=None)[0])
S[S.SNP.isin(h)].to_csv(f"{X}/bellenguez_ldsc_in_hm3.tsv", sep="\t", index=False)
# (첫 시도의 --merge-alleles 는 패치 전 munge 오류로 실패해, HapMap3 필터를 위에서 먼저 적용했다)
subprocess.run(["python", "-W", "ignore", "munge_sumstats.py", "--sumstats", f"{X}/bellenguez_ldsc_in_hm3.tsv",
                "--signed-sumstats", "Z,0", "--out", f"{X}/bellenguez_noAPOE"], cwd=LDSC, check=True,
               env=dict(os.environ, PYTHONIOENCODING="utf-8"))
