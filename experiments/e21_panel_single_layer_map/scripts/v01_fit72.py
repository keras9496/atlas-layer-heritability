"""실험 3v v01: 형질 유형 — 106개 GWAS × 주석 72개(조직 묶음 12 × 층위 6)를 하나씩 넣은 S-LDSC (PROTOCOL_exp3v).
 모형은 주 지도의 행렬(3d y04 단계 B, model D1)과 같다: baselineLD v2.2 + L5_<묶음>_<층위> 하나. jackknife 블록도 3d처럼 ldsc 기본값(SNP 수로 나눈 200블록)이다.
 빠른 실행기: 형질마다 요약통계와 baselineLD LD 점수를 한 번만 읽어 합치고(ldsc의 _read_ld_sumstats 그대로), 마지막 열(주석의 LD 점수)과 M만 바꿔 ldsc의 Hsq를 72번 부른다.
   회귀 코드는 ldsc(ldsc_bv)의 것이고, estimate_h2의 앞부분(χ² 상한, 배열 만들기)을 그대로 옮겼다.
 단계 (python -I -X utf8 v01_fit72.py prepld | check | fit <k> <K>):
   prepld = 주석 72개의 LD 점수를 baselineLD의 SNP 순서로 한 배열에 모은다(exp3v/annld.npy, annld_snps.txt, annld_M.csv). SNP 열이 baselineLD와 같아야 한다.
   check  = 주 지도의 두 형질(HEIGHT, AD)에 72개 모형을 돌려 3d y04_B.csv의 τ*, 표준오차, h²와 같아야 한다(상대 차이 < 1e-6). 다르면 멈춘다.
   fit    = 형질 목록(3u traits_exp3u.csv의 106개)의 k번째 조각(k = 0..K−1). 재실행 안전(npz가 있으면 건너뛴다).
 출력: exp3v/runs/<id>.npz (ann, coef, coef_se, delete[72 × 200], tot, tot_se, intercept, mean_chisq, n_snp), <id>.log
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
X = DATA/"exp3"; B = X/"GRCh38"; Y = DATA/"exp3d"; AN = Y/"annot"; U = DATA/"exp3u"; V = DATA/"exp3v"; BV = DATA/"tools"/"ldsc_bv"
V.mkdir(exist_ok=True); RUN = V/"runs"; RUN.mkdir(exist_ok=True)
LIN = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
MODS = ["CAGE", "DNASE", "H3K27ac", "H3K4me1", "H3K4me3", "RNA"]
ANN = [f"L5_{g}_{m}" for g in LIN for m in MODS]; assert len(ANN) == 72
BASE = str(B/"baselineLD_v2.2"/"baselineLD."); WLD = str(B/"weights"/"weights.hm3_noMHC.")
def pref(a): return str(AN/a/f"{a}.")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
sys.path.insert(0, str(BV))
import ldscore.ldsc as LD
from ldscore import sumstats as S, regressions as reg, parse as ps

if STAGE == "prepld":
    t0 = time.time(); b = ps.ldscore(BASE, 22); snps = b.SNP.to_numpy(); n = len(snps); del b
    A = np.empty((len(ANN), n)); Ms = []
    for j, a in enumerate(ANN):
        d = ps.ldscore(pref(a), 22); assert d.shape[1] == 2 and np.array_equal(d.SNP.to_numpy(), snps), a
        A[j] = d.iloc[:, 1].to_numpy(float); assert np.isfinite(A[j]).all() and A[j].var() > 0, a
        m = ps.M(pref(a), 22, common=True); assert m.shape == (1, 1); Ms.append(float(m[0, 0]))
    np.save(V/"annld.npy", A); (V/"annld_snps.txt").write_text("\n".join(snps), encoding="utf-8")
    pd.DataFrame(dict(annot=ANN, M_5_50=Ms)).to_csv(V/"annld_M.csv", index=False)
    P(f"주석 {len(ANN)}개 × SNP {n}개의 LD 점수를 모았다 ({(time.time() - t0) / 60:.1f}분)"); (RES/"v01_prepld.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit(0)

# --- 모아 둔 주석 LD 점수 ---
ASNP = pd.Index(open(V/"annld_snps.txt", encoding="utf-8").read().split("\n")); A = np.load(V/"annld.npy", mmap_mode="r")
AM = pd.read_csv(V/"annld_M.csv").set_index("annot").M_5_50
assert A.shape == (72, len(ASNP)) and list(AM.index) == ANN and ASNP.is_unique
# --- ldsc 패치: baselineLD와 가중치 LD 점수는 한 번만 읽는다(3u u01과 같다). 블록 경계는 바꾸지 않는다 ---
_ld = ps.ldscore; CACHE = {}
def ld_cached(fh, num=None):
    if "baselineLD_v2.2" in fh or "weights.hm3" in fh:
        if (fh, num) not in CACHE: CACHE[(fh, num)] = _ld(fh, num)
        return CACHE[(fh, num)].copy()
    return _ld(fh, num)
ps.ldscore = ld_cached
OUT = {}
def eh72(args, lg):
    """ldsc sumstats.estimate_h2의 앞부분을 그대로 옮기고, Hsq를 주석 72개에 되풀이한다."""
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
    assert np.array_equal(ref_ld[:, -1], np.asarray(A[0])[pos]) and M_annot[0, -1] == AM[ANN[0]]      # 파일에서 읽은 첫 주석의 열 = 모아 둔 배열
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
    T = V/"test"; T.mkdir(exist_ok=True); REF = pd.read_csv(Path(__file__).resolve().parents[2]/"e03_single_layer_map"/"results"/"y04_B.csv")
    SRC = {"HEIGHT": DATA/"exp3b"/"sumstats"/"PASS_Height1.sumstats", "AD": X/"bellenguez_noAPOE.sumstats.gz"}; worst = 0
    for t, f in SRC.items():
        t0 = time.time(); r = run(f, T/f"CHECK_{t}"); ts, se, h2 = taustar(r); ref = REF[REF.trait == t].set_index("annot").loc[ANN]
        d1 = float(np.max(np.abs(ts - ref.tau_star.to_numpy()) / np.abs(ref.tau_star.to_numpy()))); d2 = float(np.max(np.abs(se - ref.tau_star_se.to_numpy()) / ref.tau_star_se.to_numpy()))
        d3 = float(np.max(np.abs(h2 - ref.h2_obs.to_numpy()))); worst = max(worst, d1, d2, d3)
        P(f"확인 {t}: 72개 모형 — τ* 상대 차이 최대 {d1:.1e}, 표준오차 {d2:.1e}, h² 차이 {d3:.1e} ({(time.time() - t0) / 60:.1f}분)")
    assert worst < 1e-6, "재현 실패 — 멈춘다"; P("통과"); (RES/"v01_check.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit(0)

if STAGE == "fit":
    k, K = int(sys.argv[2]), int(sys.argv[3]); TR = pd.read_csv(Path(__file__).resolve().parents[2]/"e20_gwas_panel"/"traits_exp3u.csv"); TR = TR[TR.include].reset_index(drop=True); assert len(TR) == 106
    ids = TR.id.tolist()[k::K]; t0 = time.time()
    for i, tid in enumerate(ids):
        run(U/"sumstats"/f"{tid}.sumstats.gz", RUN/tid)
        (V/f"STATUS_fit_{k}of{K}.txt").write_text(f"[{i + 1}/{len(ids)}] {tid} ({(time.time() - t0) / 60:.1f}분)", encoding="utf-8")
    P(f"fit 조각 {k}/{K}: 형질 {len(ids)}개 × 72모형 끝, {(time.time() - t0) / 60:.1f}분"); (RES/f"v01_fit_{k}of{K}.txt").write_text("\n".join(log), encoding="utf-8")
