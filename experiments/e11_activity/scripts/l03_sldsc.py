"""실험 3l l03: 활성도 조건화 S-LDSC (PROTOCOL_exp3l_2026-10-17, 2026-10-18 고정). 3j j02·3m m04b와 같은 ldsc 호출과 정렬 블록 패치.
 모형(형질 t, 자기 계통 g; baselineLD v2.2에 더함):
   EA  = 효과 L5_{g}_{CAGE,DNASE,H3K27ac} (3d) + 활성도 A5_{g} (l02)        ← 주
   A   = 활성도 A5_{g}
   EAM = 효과 3 + 활성도 3 + 실측 BM_{g} (3e)
   E   = 효과 3 = 3j의 T3all (이미 실행, 다시 돌리지 않는다)
 단계: python l03_sldsc.py pretest | dry | main
   pretest = 코드 경로 시험(결과 아님): 활성도 자리에 3i 세포형 일치 주석 M3_{g}를 대역으로 넣고 가짜 형질 1개(AD Z 섞음)로 EA·A·EAM 한 번씩, exp3l/test
   dry     = 가짜 형질 1개로 실제 A5_{g}를 넣은 EA 한 번 (형식 확인, 결과 아님), exp3l/test
   main    = 14형질 × (EA, A, EAM) = 42회, exp3l/runs
출력: exp3l/runs[/test]/<형질>__<모형>.npz (names, coef, part_delete, tot, h2_log, n_snp, sep), .log; results/l03_<단계>.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import sys, re, io, time, contextlib
import numpy as np, pandas as pd
from pathlib import Path
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); B = X/"GRCh38"; Y = Path(f"{DATA}/exp3d"); F = Path(f"{DATA}/exp3f")
I = Path(f"{DATA}/exp3i"); W = Path(f"{DATA}/exp3e"); J = Path(f"{DATA}/exp3j"); L = Path(f"{DATA}/exp3l")
BV = Path(f"{DATA}/tools/ldsc_bv"); SS = Path(f"{DATA}/exp3b/sumstats")
RUN = L/"test" if STAGE in ("pretest", "dry") else L/"runs"; RUN.mkdir(parents=True, exist_ok=True)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
TRAITS = {"AD": X/"bellenguez_noAPOE.sumstats.gz", "MONO": Path(f"{DATA}/exp3b/monocyte_count.sumstats.gz"), "EA": SS/"PASS_Years_of_Education2.sumstats",
          "SCZ": SS/"PASS_Schizophrenia.sumstats", "BIP": SS/"PASS_Bipolar_Disorder.sumstats", "NEURO": SS/"PASS_Neuroticism.sumstats", "BMI": SS/"PASS_BMI1.sumstats",
          "IBD": SS/"PASS_IBD.sumstats", "RA": SS/"PASS_Rheumatoid_Arthritis.sumstats", "LDL": SS/"PASS_LDL.sumstats", "HEIGHT": SS/"PASS_Height1.sumstats",
          "BMD": F/"sumstats"/"BMD.sumstats.gz", "NEUT": F/"sumstats"/"NEUT.sumstats.gz", "TC": F/"sumstats"/"TC.sumstats.gz"}   # 3j j02와 같다
OWN = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
       "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal", "BMD": "stromal", "NEUT": "myeloid", "TC": "liver"}
LAY = ["CAGE", "DNASE", "H3K27ac"]
def act(g): return str(I/"annot"/f"M3_{g}"/f"M3_{g}.") if STAGE == "pretest" else str(L/"annot"/f"A5_{g}"/f"A5_{g}.")
def refs(t, model):
    g = OWN.get(t, "myeloid"); eff = [str(Y/"annot"/f"L5_{g}_{m}"/f"L5_{g}_{m}.") for m in LAY]
    return {"EA": eff + [act(g)], "A": [act(g)], "EAM": eff + [act(g), str(W/"annot"/f"BM_{g}"/f"BM_{g}.")]}[model]
NADD = {"EA": 6, "A": 3, "EAM": 9}
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
    ST["names"] = list(ref_ld_cnames); cm = max(0.001 * ss.N.max(), 80)
    ST["sep"] = separators(ss.SNP.to_numpy()[np.ravel(ss.Z ** 2 < cm)]); return out
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
def run(t, model, ss=None):
    out = RUN/f"{t}__{model}"
    if Path(f"{out}.npz").exists(): return dict(np.load(f"{out}.npz", allow_pickle=False))
    ST.update(sep=None, names=None, hsq=None)
    sys.argv = ["ldsc", "--h2", str(ss or TRAITS[t]), "--ref-ld-chr", ",".join([str(B/"baselineLD_v2.2"/"baselineLD.")] + refs(t, model)),
                "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."), "--out", str(out)]
    with contextlib.redirect_stdout(io.StringIO()): LD.main()
    for l in LOGS: l.log_fh.close()
    LOGS.clear(); h = ST["hsq"]
    h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log", encoding="utf-8").read()).group(1))
    r = dict(names=np.array(ST["names"]), coef=np.asarray(h.coef).ravel(), part_delete=np.asarray(h.part_delete_values), tot=float(np.asarray(h.tot).ravel()[0]),
             h2_log=h2, n_snp=int(ST["sep"][-1]), sep=np.asarray(ST["sep"]))
    assert r["part_delete"].shape == (NB, len(r["coef"])), r["part_delete"].shape
    assert len(r["coef"]) == 97 + NADD[model], (model, len(r["coef"]))
    np.savez(f"{out}.npz", **r); return r
t0 = time.time()
if STAGE in ("pretest", "dry"):
    fk = L/"test"/"FAKE1.sumstats.gz"
    if not fk.exists():
        d = pd.read_csv(TRAITS["AD"], sep="\t"); d["Z"] = np.random.default_rng(20261018).permutation(d.Z.to_numpy()); d.to_csv(fk, sep="\t", index=False, compression="gzip")
    tag = "FAKEpre" if STAGE == "pretest" else "FAKE1"      # 대역 시험과 실제 형식 시험의 출력 이름을 나눈다
    for model in (["EA", "A", "EAM"] if STAGE == "pretest" else ["EA"]):
        r = run(tag, model, ss=fk); P(f"  {tag} {model} 완료 (주석 {len(r['coef'])}개, h2 {r['h2_log']}, SNP {r['n_snp']:,}; {(time.time() - t0) / 60:.1f}분)")
else:
    MAIN = [(t, m) for m in ["EA", "A", "EAM"] for t in TRAITS]; assert len(MAIN) == 42
    for i, (t, model) in enumerate(MAIN):
        run(t, model); el = (time.time() - t0) / 60
        msg = f"[{i + 1}/{len(MAIN)}] {t} {model} 완료 ({el:.0f}분, 남은 약 {el / (i + 1) * (len(MAIN) - i - 1):.0f}분)"; P("  " + msg)
        (L/"PIPELINE_STATUS_l03.txt").write_text(msg, encoding="utf-8")
(EXP/"results"/f"l03_{STAGE}.txt").write_text("\n".join(log), encoding="utf-8")
