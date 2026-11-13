"""실험 3u u01: 형질 확대 S-LDSC 적합 (PROTOCOL_exp3u). 주석과 LD 점수는 모두 기존 것(3c 계통 주석, 3d 층위 주석)이고 새 주석은 없다.
 실행기는 3t t08(= 3j j02)과 같다: 원래 ldsc(ldsc_bv)의 main()을 부르고 jackknife 블록 경계만 3j 정렬 경계로 바꾼다.
 모형 (형질 × 12계통):
   M1_<계통>    = baselineLD v2.2 + L5_<계통>                          (계통 지도)
   T3all_<계통> = baselineLD v2.2 + L5_<계통>_{CAGE, DNASE, H3K27ac}    (세 층위 결합, 전체 트랙)
 단계 (python -I -X utf8 u01_fit.py prep | check | fit <k> <K> | fitR <k> <K> | fitS <k> <K>):
   prep  = 요약통계 정리본(SNP A1 A2 N Z; 결측 제거, 같은 rsID는 첫 행만)을 exp3u/sumstats/<id>.sumstats.gz에 쓴다. 회귀 없음.
   check = 기존 적합 둘을 다시 돌려 계수와 delete 값이 같아야 한다(상대 차이 < 1e-9): 3j HEIGHT__T3all_adipose, 3n HEIGHT__M1_stromal. 다르면 멈춘다.
   fit   = 일 목록(형질 × 24모형)의 k번째 조각(k = 0..K−1)을 적합한다. 재실행 안전(npz가 있으면 건너뛴다).
   fitR  = 재추출 모형 T3r<01..20>_<분석 계통>(3f 주석 T3_<s>_<계통>)을 u02가 정한 형질에 적합한다.
   fitS  = 새 형질 가운데 stromal 계통 단측 P < 0.05인 형질에 stromal 재추출 모형 T3r<01..20>_stromal을 적합한다(사후).
 출력: exp3u/runs/<id>__<모형>.npz (names, coef, part_delete, tot, intercept, h2_log, n_snp, sep), .log; 실패는 <id>__<모형>.FAILED.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import sys, re, io, time, contextlib
import numpy as np, pandas as pd
from pathlib import Path
STAGE = sys.argv[1]
DATA = Path(DATA); EXP = Path(__file__).resolve().parents[1]; RES = EXP/"results"; RES.mkdir(exist_ok=True)
X = DATA/"exp3"; B = X/"GRCh38"; C = DATA/"exp3c"; Y = DATA/"exp3d"; F = DATA/"exp3f"; J = DATA/"exp3j"; N3 = DATA/"exp3n"; U = DATA/"exp3u"
RAW = U/"raw"/"sumstats_107"; SS = U/"sumstats"; BV = DATA/"tools"/"ldsc_bv"
LIN = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
LAY = ["CAGE", "DNASE", "H3K27ac"]; MODELS = [f"M1_{g}" for g in LIN] + [f"T3all_{g}" for g in LIN]
TR = pd.read_csv(EXP/"traits_exp3u.csv"); TR = TR[TR.include].reset_index(drop=True)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
if STAGE == "prep":
    SS.mkdir(parents=True, exist_ok=True); rows = []
    for r in TR.itertuples():
        d = pd.read_csv(RAW/f"{r.file}.sumstats.gz", sep=r"\s+")[["SNP", "A1", "A2", "N", "Z"]]; n0 = len(d)
        d = d.dropna(); n1 = len(d); d = d.drop_duplicates("SNP", keep="first"); n2 = len(d)
        if not (SS/f"{r.id}.sumstats.gz").exists(): d.to_csv(SS/f"{r.id}.sumstats.gz", sep="\t", index=False, compression="gzip")
        rows.append(dict(id=r.id, file=r.file, n_raw=n0, n_na=n0 - n1, n_dup=n1 - n2, n_out=n2))
    Q = pd.DataFrame(rows); Q.to_csv(RES/"u01_prep.csv", index=False); P(f"정리본 {len(Q)}개. 결측 행 합계 {int(Q.n_na.sum())}, 중복 행 합계 {int(Q.n_dup.sum())} (파일 {int((Q.n_dup > 0).sum())}개)")
    (RES/"u01_prep.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit(0)
# --- ldsc 불러오기와 패치 (3j j02, 3t t08과 같다) ---
sys.path.insert(0, str(BV))
import ldscore.ldsc as LD
from ldscore import sumstats as S, jackknife as jk, parse as ps
BOUNDS = np.load(J/"blocks.npz")["bounds"]; NB = len(BOUNDS) + 1
KEY = pd.read_parquet(J/"snp_key.parquet").drop_duplicates("SNP").set_index("SNP").key
ST = dict(sep=None, names=None, hsq=None); _ld = ps.ldscore; CACHE = {}
def ld_cached(fh, num=None):
    if "baselineLD_v2.2" in fh or "weights.hm3" in fh:
        if (fh, num) not in CACHE: CACHE[(fh, num)] = _ld(fh, num)
        return CACHE[(fh, num)].copy()
    return _ld(fh, num)
ps.ldscore = ld_cached
class EmptyBlock(Exception): pass
def separators(snps):
    k = KEY.reindex(snps).to_numpy(); assert not np.isnan(k).any() and (np.diff(k) >= 0).all()
    c = np.bincount(np.searchsorted(BOUNDS, k.astype(np.int64), side="right"), minlength=NB); assert len(c) == NB
    if not (c > 0).all(): raise EmptyBlock(f"빈 블록 {int((c == 0).sum())}개")
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
def refs(model):
    kind, g = model.split("_", 1)
    if kind == "M1": return [str(C/"annot"/f"L5_{g}"/f"L5_{g}.")]
    if kind == "T3all": return [str(Y/"annot"/f"L5_{g}_{m}"/f"L5_{g}_{m}.") for m in LAY]
    if kind.startswith("T3r"): s = kind[3:]; return [str(F/"annot"/f"T3_{s}_{g}"/f"T3_{s}_{g}.")]   # 3f 재추출 주석(세 층위가 한 파일)
    raise ValueError(model)
def run(sumstats, model, out):
    if Path(f"{out}.npz").exists(): return dict(np.load(f"{out}.npz", allow_pickle=False))
    ST.update(sep=None, names=None, hsq=None)
    sys.argv = ["ldsc", "--h2", str(sumstats), "--ref-ld-chr", ",".join([str(B/"baselineLD_v2.2"/"baselineLD.")] + refs(model)), "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."), "--out", str(out)]
    try:
        with contextlib.redirect_stdout(io.StringIO()): LD.main()
    finally:
        for l in LOGS: l.log_fh.close()
        LOGS.clear()
    h = ST["hsq"]; txt = open(f"{out}.log", encoding="utf-8").read(); h2 = float(re.search(r"Total Observed scale h2: ([\d.eE+-]+)", txt).group(1))
    r = dict(names=np.array(ST["names"]), coef=np.asarray(h.coef).ravel(), part_delete=np.asarray(h.part_delete_values), tot=float(np.asarray(h.tot).ravel()[0]),
             intercept=float(np.asarray(h.intercept).ravel()[0]), h2_log=h2, n_snp=int(ST["sep"][-1]), sep=np.asarray(ST["sep"]))
    assert r["part_delete"].shape == (NB, len(r["coef"])); np.savez(f"{out}.npz", **r); return r
def rel(a, b): return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))
if STAGE == "check":
    T = U/"test"; T.mkdir(parents=True, exist_ok=True)
    a = run(DATA/"exp3b"/"sumstats"/"PASS_Height1.sumstats", "T3all_adipose", T/"CHECK_HEIGHT__T3all_adipose"); b = np.load(J/"runs"/"HEIGHT__T3all_adipose.npz")
    d1, d2 = rel(a["coef"], b["coef"]), rel(a["part_delete"], b["part_delete"]); P(f"확인 1: 3j HEIGHT T3all_adipose — 계수 상대 차이 {d1:.1e}, delete 값 {d2:.1e}")
    a = run(N3/"sumstats"/"HEIGHT.sumstats.gz", "M1_stromal", T/"CHECK_HEIGHT3n__M1_stromal"); b = np.load(N3/"runs"/"HEIGHT__M1_stromal.npz")
    d3, d4 = rel(a["coef"], b["coef"]), rel(a["part_delete"], b["part_delete"]); P(f"확인 2: 3n HEIGHT M1_stromal — 계수 상대 차이 {d3:.1e}, delete 값 {d4:.1e}")
    assert max(d1, d2, d3, d4) < 1e-9, "재현 실패 — 멈춘다"; P("통과"); (RES/"u01_check.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit(0)
if STAGE in ("fit", "fitR", "fitS"):
    k, K = int(sys.argv[2]), int(sys.argv[3]); RUN = U/"runs"; RUN.mkdir(parents=True, exist_ok=True)
    if STAGE == "fit": jobs = [(r.id, m) for r in TR.itertuples() for m in MODELS]
    elif STAGE == "fitS":   # 사후: 새 형질 가운데 stromal 계통 단측 P < 0.05인 형질에 stromal 재추출 모형 20회
        A = pd.read_csv(RES/"u02_analysis_lineage.csv"); A = A[(A.status == "new") & (A.stromal_p1 < .05)]
        jobs = [(r.id, f"T3r{s:02d}_stromal") for r in A.itertuples() for s in range(1, 21)]
    else:   # 재추출 모형: 분석 계통이 재추출 주석이 있는 다섯 계통인 형질 × 20회 (u02의 분석 계통 표)
        A = pd.read_csv(RES/"u02_analysis_lineage.csv"); A = A[A.eligible & A.analysis_lineage.isin(["CNS_tissue", "liver", "lymphoid", "myeloid", "stromal"])]
        jobs = [(r.id, f"T3r{s:02d}_{r.analysis_lineage}") for r in A.itertuples() for s in range(1, 21)]
    jobs = jobs[k::K]; t0 = time.time(); nf = 0
    for i, (tid, m) in enumerate(jobs):
        out = RUN/f"{tid}__{m}"; fail = Path(f"{out}.FAILED.txt")
        if fail.exists(): nf += 1; continue
        try: run(SS/f"{tid}.sumstats.gz", m, out)
        except EmptyBlock as e: fail.write_text(str(e), encoding="utf-8"); nf += 1
        if (i + 1) % 25 == 0 or i + 1 == len(jobs):
            (U/f"PIPELINE_STATUS_u01_{STAGE}_{k}of{K}.txt").write_text(f"[{i + 1}/{len(jobs)}] {tid} {m} ({(time.time() - t0) / 60:.1f}분, 실패 {nf})", encoding="utf-8")
    P(f"{STAGE} 조각 {k}/{K}: {len(jobs)}개 끝, 실패 {nf}개, {(time.time() - t0) / 60:.1f}분"); (RES/f"u01_{STAGE}_{k}of{K}.txt").write_text("\n".join(log), encoding="utf-8")
