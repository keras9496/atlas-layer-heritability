"""실험 3b x06: 요약통계 준비 (PROTOCOL_exp3b_2026-09-28).
 1) APOE 포함 Bellenguez: 00_prepare_sumstats 와 같은 절차에서 APOE 제외만 뺀다 → exp3b/bellenguez_withAPOE.sumstats.gz
 2) 대형 좌위 제외: R_gws(stage I p < 5e-8 SNV ±1 Mb 병합) ∪ R_75(75 좌위 선도 ±1 Mb) → 회귀 SNP 에서 제거
    → exp3b/bellenguez_noAPOE_noGWS.sumstats.gz, 부: BIN1(rs6733839 ±1 Mb)만 제거 → bellenguez_noAPOE_noBIN1.sumstats.gz
 3) 단핵구 수(Vuckovic 2020, GCST90002393 조화 파일): rsid·effect/other allele·beta/se → Z, N 408,112, HapMap3 no-MHC → munge
    (다운로드가 끝나지 않았으면 건너뛴다: md5 확인)
 PASS_* 5개는 LDSC 형식 그대로 쓰므로 준비하지 않는다.
출력: <DATA>/exp3b/, results/x06_prepare_3b.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import os, subprocess, hashlib, sys
import numpy as np, pandas as pd
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); Y = Path(f"{DATA}/exp3b"); SS = Y/"sumstats"
LDSC = f"{DATA}/tools/ldsc"; UP = Path(__file__).resolve().parents[1]/"inputs"
env = dict(os.environ, PYTHONIOENCODING="utf-8")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
def munge(inp, out):
    subprocess.run(["python", "-W", "ignore", "munge_sumstats.py", "--sumstats", str(inp), "--signed-sumstats", "Z,0", "--out", str(out)], cwd=LDSC, check=True, env=env, capture_output=True)
    t = open(f"{out}.log", encoding="utf-8").read(); P(f"  munge {Path(out).name}: {t[t.find('Writing'):].splitlines()[0]}"); P("  " + " | ".join(l.split(' - INFO - ')[-1] for l in t.splitlines() if "Mean chi^2" in l or "Max chi^2" in l))
HM3 = set(pd.read_csv(X/"hm3_no_MHC.list.txt", header=None)[0])
ONLY_MONO = os.environ.get("X06_ONLY_MONO") == "1"   # 단핵구 단계만 (x07 실행 중 AD 파일 덮어쓰기 방지)

if not ONLY_MONO:
    # ---- 1) APOE 포함 Bellenguez
    S = pd.read_csv(f"{DATA}/sumstats/GCST90027158_buildGRCh38.tsv.gz", sep="\t",
                    usecols=["variant_id", "p_value", "chromosome", "base_pair_location", "effect_allele", "other_allele",
                             "effect_allele_frequency", "beta", "standard_error", "n_cases", "n_controls"], dtype={"chromosome": str})
    S = S[(S.effect_allele.str.len() == 1) & (S.other_allele.str.len() == 1)]
    isAPOE = (S.chromosome == "19") & (S.base_pair_location >= 43_900_000) & (S.base_pair_location <= 45_900_000)
    S["N"] = S.n_cases + S.n_controls; S["Z"] = S.beta / S.standard_error
    out = S[S.variant_id.isin(HM3)][["variant_id", "effect_allele", "other_allele", "Z", "p_value", "N", "effect_allele_frequency"]].rename(
        columns={"variant_id": "SNP", "effect_allele": "A1", "other_allele": "A2", "p_value": "P", "effect_allele_frequency": "FRQ"})
    P(f"[APOE 포함] HapMap3 SNP {len(out):,} (그중 APOE 영역 {int(S[S.variant_id.isin(HM3)].pipe(lambda d: ((d.chromosome == '19') & (d.base_pair_location.between(43_900_000, 45_900_000))).sum())):,})")
    out.to_csv(Y/"bellenguez_withAPOE_in_hm3.tsv", sep="\t", index=False); munge(Y/"bellenguez_withAPOE_in_hm3.tsv", Y/"bellenguez_withAPOE")

    # ---- 2) 대형 좌위 제외 영역
    g = S[(~isAPOE) & (S.p_value < 5e-8)][["chromosome", "base_pair_location"]].rename(columns={"chromosome": "chr", "base_pair_location": "pos"})
    L = pd.concat([pd.read_csv(UP/"loci.csv")[["chr", "pos"]], pd.read_csv(UP/"rep25_loci.csv")[["chr", "pos"]]]); L["chr"] = L.chr.astype(str)
    P(f"[제외 영역] stage I p < 5e-8 SNV {len(g):,}개(APOE 제외), 75 좌위 선도 {len(L)}개")
    def merge_iv(df):
        rows = []
        for c, d in df.groupby("chr"):
            iv = sorted(zip(d.pos - 1_000_000, d.pos + 1_000_000)); s0, e0 = iv[0]
            for s, e in iv[1:]:
                if s <= e0: e0 = max(e0, e)
                else: rows.append((c, max(s0, 0), e0)); s0, e0 = s, e
            rows.append((c, max(s0, 0), e0))
        return pd.DataFrame(rows, columns=["chr", "start", "end"])
    R = merge_iv(pd.concat([g, L])); R.to_csv(Y/"exclude_regions_gws_75.tsv", sep="\t", index=False)
    P(f"  R_gws ∪ R_75 병합 구간 {len(R)}개, 총 {(R.end - R.start).sum() / 1e6:.0f} Mb (R_gws 단독 {len(merge_iv(g))}개)")
    pos = S[S.variant_id.isin(HM3)].drop_duplicates("variant_id").set_index("variant_id")[["chromosome", "base_pair_location"]]
    def in_regions(snps, R):
        p = pos.reindex(snps); hit = np.zeros(len(p), bool)
        for c, s, e in R.itertuples(index=False): hit |= ((p.chromosome == c) & p.base_pair_location.between(s, e)).to_numpy()
        return hit
    base = pd.read_csv(X/"bellenguez_noAPOE.sumstats.gz", sep="\t")
    for name, RR in [("noGWS", R), ("noBIN1", pd.DataFrame([("2", 127135234 - 1_000_000, 127135234 + 1_000_000)], columns=["chr", "start", "end"]))]:
        hit = in_regions(base.SNP, RR); k = base[~hit]
        k.to_csv(Y/f"bellenguez_noAPOE_{name}.sumstats.gz", sep="\t", index=False, compression="gzip")
        P(f"  {name}: 회귀 SNP {len(base):,} → {len(k):,} (제거 {int(hit.sum()):,}); mean χ² {(base.Z ** 2).mean():.3f} → {(k.Z ** 2).mean():.3f}, max χ² {(k.Z ** 2).max():.1f}")

# ---- 3) 단핵구 수
f = SS/"GCST90002393.h.tsv.gz"; md5 = SS/"GCST90002393_md5.txt"
def md5sum(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""): h.update(b)
    return h.hexdigest()
ok = f.exists() and md5.exists() and any(md5sum(f) in l for l in open(md5) if "h.tsv.gz" in l and "tbi" not in l)
if not ok: P("[단핵구] 다운로드 미완료 또는 md5 불일치 → 건너뜀");
else:
    M = pd.read_csv(f, sep="\t", usecols=["rsid", "effect_allele", "other_allele", "beta", "standard_error", "p_value", "effect_allele_frequency", "INFO"], low_memory=False)
    M = M[M.rsid.isin(HM3) & (M.effect_allele.str.len() == 1) & (M.other_allele.str.len() == 1)].copy()
    M["Z"] = M.beta / M.standard_error; M["N"] = 408112
    M = M.rename(columns={"rsid": "SNP", "effect_allele": "A1", "other_allele": "A2", "p_value": "P", "effect_allele_frequency": "FRQ"})[["SNP", "A1", "A2", "Z", "P", "N", "FRQ", "INFO"]]
    P(f"[단핵구] HapMap3 SNV {len(M):,}"); M.to_csv(Y/"monocyte_in_hm3.tsv", sep="\t", index=False); munge(Y/"monocyte_in_hm3.tsv", Y/"monocyte_count")
(EXP/"results"/"x06_prepare_3b.txt").write_text("\n".join(log), encoding="utf-8")
