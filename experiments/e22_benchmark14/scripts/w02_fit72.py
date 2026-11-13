"""실험 3w w02: 벤치마크 14형질 가운데 72칸 지도(조직 묶음 12 × 층위 6, 주석 하나씩)가 없는 아홉 GWAS — 더 큰 GWAS 여섯과 검증 형질 셋 (PROTOCOL_exp3w 2절).
 모형과 실행기는 3v v01_fit72.py와 같다(baselineLD v2.2 + L5_<묶음>_<층위> 하나; jackknife는 3d·3v처럼 ldsc 기본 200블록; 요약통계와 baselineLD를 한 번 읽고 마지막 열과 M만 바꿔 Hsq를 72번).
 주석 72개의 LD 점수 배열은 3v가 모아 둔 것(exp3v/annld.npy)을 읽기만 한다. 나머지 다섯 형질(EA, NEURO, MONO, IBD, RA)의 72칸은 3d y04_B.csv에 있다.
 단계 (python -I -X utf8 w02_fit72.py check | fit <k> <K>):
   check = 3d에 있는 형질 하나(EA)를 다시 돌려 y04_B.csv의 τ*, 표준오차와 같아야 한다(상대 차이 < 1e-6).
 출력: exp3w/runs72/<형질 열쇠>.npz (ann, coef, coef_se, delete[72 × 200], tot, tot_se, intercept, mean_chisq, n_snp), .log
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import os
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_k, "2")
import sys, io, time, copy, contextlib
import numpy as np, pandas as pd
from pathlib import Path
STAGE = sys.argv[1]
DATA = Path(DATA); EXP = Path(__file__).resolve().parents[1]; RES = EXP/"results"; RES.mkdir(exist_ok=True)
X = DATA/"exp3"; B = X/"GRCh38"; Y = DATA/"exp3d"; AN = Y/"annot"; F = DATA/"exp3f"; N3 = DATA/"exp3n"; R3 = DATA/"exp3r"; V = DATA/"exp3v"; W = DATA/"exp3w"; BV = DATA/"tools"/"ldsc_bv"
RUN = W/"runs72"; RUN.mkdir(parents=True, exist_ok=True)
LIN = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
MODS = ["CAGE", "DNASE", "H3K27ac", "H3K4me1", "H3K4me3", "RNA"]
ANN = [f"L5_{g}_{m}" for g in LIN for m in MODS]; assert len(ANN) == 72
SS = {"N_HEIGHT": N3/"sumstats"/"HEIGHT.sumstats.gz", "N_BMI": N3/"sumstats"/"BMI.sumstats.gz", "N_SCZ": N3/"sumstats"/"SCZ.sumstats.gz", "N_BIP": N3/"sumstats"/"BIP.sumstats.gz",
      "N_LDL": N3/"sumstats"/"LDL.sumstats.gz", "N_AD": R3/"sumstats"/"N_AD.sumstats.gz", "BMD": F/"sumstats"/"BMD.sumstats.gz", "NEUT": F/"sumstats"/"NEUT.sumstats.gz", "TC": F/"sumstats"/"TC.sumstats.gz"}
BASE = str(B/"baselineLD_v2.2"/"baselineLD."); WLD = str(B/"weights"/"weights.hm3_noMHC.")
def pref(a): return str(AN/a/f"{a}.")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
sys.path.insert(0, str(BV))
import ldscore.ldsc as LD
from ldscore import sumstats as S, regressions as reg, parse as ps
ASNP = pd.Index(open(V/"annld_snps.txt", encoding="utf-8").read().split("\n")); A = np.load(V/"annld.npy", mmap_mode="r")
AM = pd.read_csv(V/"annld_M.csv").set_index("annot").M_5_50
assert A.shape == (72, len(ASNP)) and list(AM.index) == ANN and ASNP.is_unique
_ld = ps.ldscore; CACHE = {}
def ld_cached(fh, num=None):
    if "baselineLD_v2.2" in fh or "weights.hm3" in fh:
        if (fh, num) not in CACHE: CACHE[(fh, num)] = _ld(fh, num)
        return CACHE[(fh, num)].copy()
    return _ld(fh, num)
ps.ldscore = ld_cached
OUT = {}
def eh72(args, lg):
    """ldsc sumstats.estimate_h2의 앞부분을 그대로 옮기고, Hsq를 주석 72개에 되풀이한다(3v v01과 같다)."""
    args = copy.deepcopy(args)
    assert args.samp_prev is None and args.intercept_h2 is None and not args.no_intercept and args.chisq_max is None and args.two_step is None
    M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols = S._read_ld_sumstats(args, lg, args.h2)
    assert not np.asarray(novar_cols).any()
    n_snp = len(sumstats); n_blocks = min(n_snp, args.n_blocks); n_annot = len(ref_ld_cnames); assert n_annot > 1 and n_blocks == 200
    chisq_max = max(0.001 * sumstats.N.max(), 80)
    chisq = np.array(sumstats.Z ** 2).reshape((n_snp, 1)); ii = np.ravel(chisq < chisq_max)
    sumstats = sumstats.iloc[ii, :]; n_snp = int(np.sum(ii))
    s = lambda x: np.array(x).reshape((n_snp, 1))
    ref_ld = np.array(sumstats[ref_ld_cnames]); chisq = chisq[ii].reshape((n_snp, 1)); w = s(sumstats[w_ld_cname]); N = s(sumstats.N)
    pos = ASNP.get_indexer(sumstats.SNP.to_numpy()); assert (pos >= 0).all()
    assert np.array_equal(ref_ld[:, -1], np.asarray(A[0])[pos]) and M_annot[0, -1] == AM[ANN[0]]
    r = dict(coef=[], coef_se=[], delete=[], tot=[], tot_se=[], intercept=[], mean_chisq=[]); first = None
    for j, a in enumerate(ANN):
        x = ref_ld.copy(); x[:, -1] = np.asarray(A[j])[pos]; M = M_annot.copy(); M[0, -1] = AM[a]
        h = reg.Hsq(chisq, x, w, N, M, n_blocks=n_blocks, intercept=None, twostep=None, old_weights=True)
        pd_ = np.asarray(h.part_delete_values); assert pd_.shape == (200, n_annot)
        r["coef"].append(float(np.asarray(h.coef).ravel()[-1])); r["coef_se"].append(float(np.asarray(h.coef_se).ravel()[-1])); r["delete"].append(pd_[:, -1].copy())
        r["tot"].append(float(np.asarray(h.tot).ravel()[0])); r["tot_se"].append(float(np.asarray(h.tot_se).ravel()[0]))
        r["intercept"].append(float(np.asarray(h.intercept).ravel()[0])); r["mean_chisq"].append(float(np.asarray(h.mean_chisq).ravel()[0]))
        if j == 0: first = h
    OUT.update({k: np.asarray(v) for k, v in r.items()}); OUT.update(ann=np.array(ANN), n_snp=n_snp)
    return first
S.estimate_h2 = eh72
LOGS = []; _li = LD.Logger.__init__
def li(self, fh): _li(self, fh); LOGS.append(self)
LD.Logger.__init__ = li
def run(sumstats, out):
    if Path(f"{out}.npz").exists(): return dict(np.load(f"{out}.npz", allow_pickle=False))
    OUT.clear()
    sys.argv = ["ldsc", "--h2", str(sumstats), "--ref-ld-chr", BASE + "," + pref(ANN[0]), "--w-ld-chr", WLD, "--out", str(out)]
    try:
        with contextlib.redirect_stdout(io.StringIO()): LD.main()
    finally:
        for l in LOGS: l.log_fh.close()
        LOGS.clear()
    assert len(OUT.get("coef", [])) == 72, "적합 결과가 없다"
    r = dict(OUT); np.savez(f"{out}.npz", **r); return r
SD = pd.read_csv(Y/"annot_sd_common.csv", index_col=0).sd_common
M_TOT = sum(float(open(B/"baselineLD_v2.2"/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
def taustar(r):
    h2 = np.round(r["tot"], 4); f = np.array([SD[a] for a in ANN]) * M_TOT / h2
    return r["coef"] * f, r["coef_se"] * f, h2
if STAGE == "check":
    T = W/"test"; T.mkdir(parents=True, exist_ok=True); REF = pd.read_csv(EXP.parent/"e03_single_layer_map"/"results"/"y04_B.csv")
    t0 = time.time(); r = run(DATA/"exp3b"/"sumstats"/"PASS_Years_of_Education2.sumstats", T/"CHECK72_EA"); ts, se, h2 = taustar(r); ref = REF[REF.trait == "EA"].set_index("annot").loc[ANN]
    d1 = float(np.max(np.abs(ts - ref.tau_star.to_numpy()) / np.abs(ref.tau_star.to_numpy()))); d2 = float(np.max(np.abs(se - ref.tau_star_se.to_numpy()) / ref.tau_star_se.to_numpy()))
    P(f"확인 EA: 72개 모형 — τ* 상대 차이 최대 {d1:.1e}, 표준오차 {d2:.1e} ({(time.time() - t0) / 60:.1f}분)"); assert max(d1, d2) < 1e-6, "재현 실패 — 멈춘다"
    P("통과"); (RES/"w02_check.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit(0)
if STAGE == "fit":
    k, K = int(sys.argv[2]), int(sys.argv[3]); keys = list(SS)[k::K]; t0 = time.time()
    for t in keys: run(SS[t], RUN/t); P(f"{t}: 72개 모형 끝 ({(time.time() - t0) / 60:.1f}분)")
    (RES/f"w02_fit_{k}of{K}.txt").write_text("\n".join(log), encoding="utf-8")
