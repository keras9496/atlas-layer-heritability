"""실험 3o o02: 연속(CONT)·구간(BIN) 모형 S-LDSC (PROTOCOL_exp3o_2026-10-21, 2026-10-21 고정). 순차 실행, 한 프로세스.
 실행기는 3j j02 / 3n n02와 같다(원래 ldsc의 main(), 같은 회귀 인자, jackknife 블록만 3j 정렬 경계).
   CONT = baselineLD + C3_{g} (연속 3층위)   BIN = baselineLD + B12_{g} (구간 4 × 3층위)
 형질: 원래 GWAS 14개(3j와 같은 파일) + 3n 새 GWAS 다섯 개(N_ 접두어). 계통은 주 지도와 같다.
 단계 (python o02_sldsc.py dry | main):
   dry  = 가짜 형질(3j FAKE1, AD Z 섞음, 결과 아님) × 골수계 CONT·BIN — 주석 수(3, 12)와 형식 확인
   main = 19형질 × 2 = 38회. 재실행 안전: npz가 있으면 건너뛴다.
출력: exp3o/runs[/test]/<형질>__<모형>.npz, .log; results/o02_<단계>.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import sys, re, io, time, contextlib
import numpy as np, pandas as pd
from pathlib import Path
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); B = X/"GRCh38"; F = Path(f"{DATA}/exp3f"); J = Path(f"{DATA}/exp3j")
N3 = Path(f"{DATA}/exp3n"); O = Path(f"{DATA}/exp3o"); BV = Path(f"{DATA}/tools/ldsc_bv"); SS = Path(f"{DATA}/exp3b/sumstats")
RUN = O/"test" if STAGE == "dry" else O/"runs"; RUN.mkdir(parents=True, exist_ok=True)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
TRAITS = {"AD": X/"bellenguez_noAPOE.sumstats.gz", "MONO": Path(f"{DATA}/exp3b/monocyte_count.sumstats.gz"), "EA": SS/"PASS_Years_of_Education2.sumstats",
          "SCZ": SS/"PASS_Schizophrenia.sumstats", "BIP": SS/"PASS_Bipolar_Disorder.sumstats", "NEURO": SS/"PASS_Neuroticism.sumstats", "BMI": SS/"PASS_BMI1.sumstats",
          "IBD": SS/"PASS_IBD.sumstats", "RA": SS/"PASS_Rheumatoid_Arthritis.sumstats", "LDL": SS/"PASS_LDL.sumstats", "HEIGHT": SS/"PASS_Height1.sumstats",
          "BMD": F/"sumstats"/"BMD.sumstats.gz", "NEUT": F/"sumstats"/"NEUT.sumstats.gz", "TC": F/"sumstats"/"TC.sumstats.gz"}   # 3j j02와 같다
TRAITS.update({f"N_{t}": N3/"sumstats"/f"{t}.sumstats.gz" for t in ["HEIGHT", "BMI", "SCZ", "BIP", "LDL"]})
OWN = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
       "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal", "BMD": "stromal", "NEUT": "myeloid", "TC": "liver", "FAKE1": "myeloid"}
OWN.update({f"N_{t}": OWN[t] for t in ["HEIGHT", "BMI", "SCZ", "BIP", "LDL"]})
def refs(t, model):
    g = OWN[t]; grp = {"CONT": f"C3_{g}", "BIN": f"B12_{g}"}[model]; return [str(O/"annot"/grp/f"{grp}.")]
NADD = {"CONT": 3, "BIN": 12}
MAIN = [(t, m) for m in ["CONT", "BIN"] for t in TRAITS]; assert len(MAIN) == 38
# --- ldsc 불러오기와 패치 (3j j02와 같다) ---
sys.path.insert(0, str(BV))
import ldscore.ldsc as LD
from ldscore import sumstats as S, jackknife as jk, parse as ps
BL = np.load(J/"blocks.npz"); BOUNDS = BL["bounds"]; NB = len(BOUNDS) + 1
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
    k = KEY.reindex(snps).to_numpy(); assert not np.isnan(k).any(), "위치 없는 SNP"
    assert (np.diff(k) >= 0).all(), "회귀 SNP가 위치 순서가 아니다"
    c = np.bincount(np.searchsorted(BOUNDS, k.astype(np.int64), side="right"), minlength=NB); assert len(c) == NB and (c > 0).all(), "빈 블록"
    return np.concatenate([[0], np.cumsum(c)]).astype(int)
_rls = S._read_ld_sumstats
def rls(args, lg, fh, alleles=False, dropna=True):
    out = _rls(args, lg, fh, alleles=alleles, dropna=dropna); M_annot, w_ld_cname, ref_ld_cnames, ss, novar = out
    ST["names"] = list(ref_ld_cnames)
    cm = max(0.001 * ss.N.max(), 80) if len(ref_ld_cnames) > 1 else None
    snps = ss.SNP.to_numpy()[np.ravel(ss.Z ** 2 < cm)] if cm is not None else ss.SNP.to_numpy()
    ST["sep"] = separators(snps); return out
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
def run(t, model, ss):
    out = RUN/f"{t}__{model}"
    if Path(f"{out}.npz").exists(): return dict(np.load(f"{out}.npz", allow_pickle=False))
    ST.update(sep=None, names=None, hsq=None)
    sys.argv = ["ldsc", "--h2", str(ss), "--ref-ld-chr", ",".join([str(B/"baselineLD_v2.2"/"baselineLD.")] + refs(t, model)), "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."), "--out", str(out)]
    with contextlib.redirect_stdout(io.StringIO()): LD.main()
    for l in LOGS: l.log_fh.close()
    LOGS.clear(); h = ST["hsq"]
    h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log", encoding="utf-8").read()).group(1))
    r = dict(names=np.array(ST["names"]), coef=np.asarray(h.coef).ravel(), part_delete=np.asarray(h.part_delete_values), tot=float(np.asarray(h.tot).ravel()[0]),
             h2_log=h2, n_snp=int(ST["sep"][-1]), sep=np.asarray(ST["sep"]))
    assert r["part_delete"].shape == (200, len(r["coef"])), r["part_delete"].shape
    assert len(r["coef"]) == 97 + NADD[model], (model, len(r["coef"]))
    np.savez(f"{out}.npz", **r); return r
t0 = time.time()
if STAGE == "dry":
    for model in ["CONT", "BIN"]:
        r = run("FAKE1", model, J/"test"/"FAKE1.sumstats.gz"); P(f"  FAKE1 {model}: 주석 {len(r['coef'])}개 (97 + {NADD[model]}), 회귀 SNP {r['n_snp']:,} — 형식 확인 (결과 아님)")
    P("형식 확인 통과")
else:
    for i, (t, model) in enumerate(MAIN):
        run(t, model, TRAITS[t]); el = (time.time() - t0) / 60
        msg = f"[{i + 1}/{len(MAIN)}] {t} {model} 완료 ({el:.0f}분, 남은 약 {el / (i + 1) * (len(MAIN) - i - 1):.0f}분)"
        P("  " + msg); (O/"PIPELINE_STATUS.txt").write_text(msg, encoding="utf-8")
(EXP/"results"/f"o02_{STAGE}.txt").write_text("\n".join(log), encoding="utf-8")
