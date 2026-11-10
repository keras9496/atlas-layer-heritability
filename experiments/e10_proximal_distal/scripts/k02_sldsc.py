"""실험 3k k02: 근위·원위 6주석 S-LDSC (PROTOCOL_exp3k_2026-10-15). 순차 실행, 한 프로세스.
 PD_{g} = baselineLD v2.2 + P_{g}_{CAGE,DNASE,H3K27ac} + D_{g}_{CAGE,DNASE,H3K27ac}. 형질 14개 × 자기 계통 = 14회.
 jackknife는 3j의 정렬 블록(exp3j/blocks.npz)을 쓴다. 호출과 패치는 3j j02와 같다(원래 ldsc main(), LstsqJackknifeFast의 separators만 바꿈).
 단계: dry = 가짜 형질(3j FAKE1, Z 섞음, 결과 아님) × PD_myeloid / main = 14회.  재실행 안전: npz가 있으면 건너뛴다.
출력: exp3k/runs[/test]/<형질>__PD_<계통>.npz (names, coef, part_delete, tot, h2_log, n_snp, sep), .log; results/k02_<단계>.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import sys, re, io, time, contextlib
import numpy as np, pandas as pd
from pathlib import Path
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); B = X/"GRCh38"; F = Path(f"{DATA}/exp3f")
J = Path(f"{DATA}/exp3j"); K = Path(f"{DATA}/exp3k"); BV = Path(f"{DATA}/tools/ldsc_bv"); SS = Path(f"{DATA}/exp3b/sumstats")
RUN = K/"test"/"runs" if STAGE == "dry" else K/"runs"; RUN.mkdir(parents=True, exist_ok=True)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
TRAITS = {"AD": X/"bellenguez_noAPOE.sumstats.gz", "MONO": Path(f"{DATA}/exp3b/monocyte_count.sumstats.gz"), "EA": SS/"PASS_Years_of_Education2.sumstats",
          "SCZ": SS/"PASS_Schizophrenia.sumstats", "BIP": SS/"PASS_Bipolar_Disorder.sumstats", "NEURO": SS/"PASS_Neuroticism.sumstats", "BMI": SS/"PASS_BMI1.sumstats",
          "IBD": SS/"PASS_IBD.sumstats", "RA": SS/"PASS_Rheumatoid_Arthritis.sumstats", "LDL": SS/"PASS_LDL.sumstats", "HEIGHT": SS/"PASS_Height1.sumstats",
          "BMD": F/"sumstats"/"BMD.sumstats.gz", "NEUT": F/"sumstats"/"NEUT.sumstats.gz", "TC": F/"sumstats"/"TC.sumstats.gz"}   # 3f·3j와 같다
OWN = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
       "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal", "BMD": "stromal", "NEUT": "myeloid", "TC": "liver"}
if STAGE == "dry": TRAITS, OWN = {"FAKE1": J/"test"/"FAKE1.sumstats.gz"}, {"FAKE1": "myeloid"}
# --- ldsc와 패치 (3j j02와 같다)
sys.path.insert(0, str(BV))
import ldscore.ldsc as LD
from ldscore import sumstats as S, jackknife as jk, parse as ps
BOUNDS = np.load(J/"blocks.npz")["bounds"]; NB = len(BOUNDS) + 1
KEY = pd.read_parquet(J/"snp_key.parquet").drop_duplicates("SNP").set_index("SNP").key
ST = dict(sep=None, names=None, hsq=None)
_ld = ps.ldscore; CACHE = {}
def ld_cached(fh, num=None):
    if "baselineLD_v2.2" in fh or "weights.hm3" in fh:
        if (fh, num) not in CACHE: CACHE[(fh, num)] = _ld(fh, num)
        return CACHE[(fh, num)].copy()
    return _ld(fh, num)
ps.ldscore = ld_cached
def separators(snps):
    k = KEY.reindex(snps).to_numpy(); assert not np.isnan(k).any(); assert (np.diff(k) >= 0).all()
    c = np.bincount(np.searchsorted(BOUNDS, k.astype(np.int64), side="right"), minlength=NB); assert len(c) == NB and (c > 0).all(), "빈 블록"
    return np.concatenate([[0], np.cumsum(c)]).astype(int)
_rls = S._read_ld_sumstats
def rls(args, lg, fh, alleles=False, dropna=True):
    out = _rls(args, lg, fh, alleles=alleles, dropna=dropna); ss = out[3]; ST["names"] = list(out[2])
    cm = max(0.001 * ss.N.max(), 80); ST["sep"] = separators(ss.SNP.to_numpy()[np.ravel(ss.Z ** 2 < cm)]); return out
S._read_ld_sumstats = rls
_LJF = jk.LstsqJackknifeFast
class LJF(_LJF):
    def __init__(self, x, y, n_blocks=None, separators=None):
        separators = list(ST["sep"]); assert separators[-1] == x.shape[0]; _LJF.__init__(self, x, y, None, separators)
jk.LstsqJackknifeFast = LJF
_eh = S.estimate_h2
def eh(args, lg): r = _eh(args, lg); ST["hsq"] = r; return r
S.estimate_h2 = eh
LOGS = []; _li = LD.Logger.__init__
def li(self, fh): _li(self, fh); LOGS.append(self)
LD.Logger.__init__ = li
t0 = time.time()
for i, t in enumerate(TRAITS):
    g = OWN[t]; out = RUN/f"{t}__PD_{g}"
    if not Path(f"{out}.npz").exists():
        ST.update(sep=None, names=None, hsq=None)
        sys.argv = ["ldsc", "--h2", str(TRAITS[t]), "--ref-ld-chr", ",".join([str(B/"baselineLD_v2.2"/"baselineLD."), str(K/"annot"/f"PD_{g}"/f"PD_{g}.")]),
                    "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."), "--out", str(out)]
        with contextlib.redirect_stdout(io.StringIO()): LD.main()
        for l in LOGS: l.log_fh.close()
        LOGS.clear(); h = ST["hsq"]
        h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log", encoding="utf-8").read()).group(1))
        r = dict(names=np.array(ST["names"]), coef=np.asarray(h.coef).ravel(), part_delete=np.asarray(h.part_delete_values), tot=float(np.asarray(h.tot).ravel()[0]),
                 h2_log=h2, n_snp=int(ST["sep"][-1]), sep=np.asarray(ST["sep"]))
        assert r["part_delete"].shape == (200, len(r["coef"])); np.savez(f"{out}.npz", **r)
    msg = f"[{i + 1}/{len(TRAITS)}] {t} PD_{g} 완료 ({(time.time() - t0) / 60:.1f}분)"; P("  " + msg)
    if STAGE == "main": (K/"PIPELINE_STATUS.txt").write_text(msg, encoding="utf-8")
(EXP/"results"/f"k02_{STAGE}.txt").write_text("\n".join(log), encoding="utf-8")
