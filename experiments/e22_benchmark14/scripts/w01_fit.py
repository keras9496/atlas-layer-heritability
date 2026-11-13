"""실험 3w w01: 벤치마크 14형질과 106개 GWAS의 빠진 S-LDSC 적합 (PROTOCOL_exp3w v1.0). 주석과 LD 점수는 모두 기존 것이다(활성도 구간 주석만 w03이 새로 만든다).
 실행기는 3u u01(= 3j j02, 3t t08)과 같다: 원래 ldsc(ldsc_bv)의 main()을 부르고 jackknife 블록 경계만 3j 정렬 경계로 바꾼다.
 모형 (baselineLD v2.2에 더하는 주석; <g> = 조직 묶음):
   M1_<g>     = L5_<g>                                  (조직 묶음 지도; 3c)
   M2_<g>     = L5_all + L5_<g>                          (모든 묶음의 큰 효과 SNP를 보정; 3c의 M2)
   T3all_<g>  = L5_<g>_{CAGE, DNASE, H3K27ac}            (세 층위 함께, 전체 트랙; 3d 주석)
   T3r<ss>_<g>= T3_<ss>_<g>                              (세 층위 함께, 층위당 세 단위 재추출; 3f 주석)
   M3_<g>     = M3_<g>                                   (세포형 일치 세 층위; 3i)
   PD_<g>     = PD_<g>                                   (근위·원위로 나눈 세 층위; 3k)
   BM_<g>     = BM_<g>                                   (실측 세 층위; 3e)
   J1C_<g>    = BM_<g> + L5_<g>_CAGE                     (실측 + 예측 CAGE; 3e의 J1)
   J1D_<g>    = BM_<g> + L5_<g>_DNASE                    (실측 + 예측 DNase; 3e의 J1own_DNASE)
   CONT_<g>   = C3_<g>,  BIN_<g> = B12_<g>               (연속 백분위, 겹치지 않는 구간; 3o)
   EAb_<g>    = L5_<g>_{CAGE, DNASE, H3K27ac} + AB_<g>   (효과 + 활성도 네 구간 × 세 층위; AB는 w03이 만든다)
 단계 (python -I -X utf8 w01_fit.py check | fit <bench|panel|actbin|panelpd> <k> <K> | jobs):
   check = 기존 적합 셋을 다시 돌려 계수와 delete 값이 같아야 한다(상대 차이 < 1e-9). 다르면 멈춘다.
   fit   = 일 목록의 k번째 조각. 재실행 안전(npz가 있으면 건너뛴다).
 출력: exp3w/runs/<형질 열쇠>__<모형>.npz (names, coef, part_delete, tot, intercept, h2_log, n_snp, sep), .log
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import sys, re, io, time, contextlib
import numpy as np, pandas as pd
from pathlib import Path
STAGE = sys.argv[1]
DATA = Path(DATA); EXP = Path(__file__).resolve().parents[1]; RES = EXP/"results"; RES.mkdir(exist_ok=True)
X = DATA/"exp3"; B = X/"GRCh38"; C = DATA/"exp3c"; Y = DATA/"exp3d"; E = DATA/"exp3e"; F = DATA/"exp3f"; I = DATA/"exp3i"; J = DATA/"exp3j"; K3 = DATA/"exp3k"
N3 = DATA/"exp3n"; O3 = DATA/"exp3o"; R3 = DATA/"exp3r"; U = DATA/"exp3u"; W = DATA/"exp3w"; BV = DATA/"tools"/"ldsc_bv"
RUN = W/"runs"; RUN.mkdir(parents=True, exist_ok=True)
LIN = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
FIVE = ["CNS_tissue", "liver", "lymphoid", "myeloid", "stromal"]; LAY = ["CAGE", "DNASE", "H3K27ac"]
PASS = DATA/"exp3b"/"sumstats"
SS = {"EA": PASS/"PASS_Years_of_Education2.sumstats", "NEURO": PASS/"PASS_Neuroticism.sumstats", "IBD": PASS/"PASS_IBD.sumstats", "RA": PASS/"PASS_Rheumatoid_Arthritis.sumstats",
      "MONO": DATA/"exp3b"/"monocyte_count.sumstats.gz", "BMD": F/"sumstats"/"BMD.sumstats.gz", "NEUT": F/"sumstats"/"NEUT.sumstats.gz", "TC": F/"sumstats"/"TC.sumstats.gz",
      "N_HEIGHT": N3/"sumstats"/"HEIGHT.sumstats.gz", "N_BMI": N3/"sumstats"/"BMI.sumstats.gz", "N_SCZ": N3/"sumstats"/"SCZ.sumstats.gz", "N_BIP": N3/"sumstats"/"BIP.sumstats.gz",
      "N_LDL": N3/"sumstats"/"LDL.sumstats.gz", "N_AD": R3/"sumstats"/"N_AD.sumstats.gz"}
OWN = {"EA": "CNS_tissue", "N_SCZ": "CNS_tissue", "N_BIP": "CNS_tissue", "NEURO": "CNS_tissue", "N_BMI": "CNS_tissue", "N_AD": "myeloid", "MONO": "myeloid", "NEUT": "myeloid",
       "IBD": "lymphoid", "RA": "lymphoid", "N_LDL": "liver", "TC": "liver", "N_HEIGHT": "stromal", "BMD": "stromal"}
B14 = list(OWN); NEW6 = [t for t in B14 if t.startswith("N_")]; KEEP8 = [t for t in B14 if not t.startswith("N_")]; assert len(B14) == 14 and len(NEW6) == 6 and set(B14) == set(SS)
DNASE_TRAITS = ["N_HEIGHT", "BMD"]                                # 층위 가운데 DNase의 τ*가 가장 큰 형질(재추출 지도)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
def refs(model):
    kind, g = model.split("_", 1); d3 = [str(Y/"annot"/f"L5_{g}_{m}"/f"L5_{g}_{m}.") for m in LAY]
    if kind == "M1": return [str(C/"annot"/f"L5_{g}"/f"L5_{g}.")]
    if kind == "M2": return [str(C/"annot"/"L5_all"/"L5_all."), str(C/"annot"/f"L5_{g}"/f"L5_{g}.")]
    if kind == "T3all": return d3
    if kind.startswith("T3r"): s = kind[3:]; return [str(F/"annot"/f"T3_{s}_{g}"/f"T3_{s}_{g}.")]
    if kind == "M3": return [str(I/"annot"/f"M3_{g}"/f"M3_{g}.")]
    if kind == "PD": return [str(K3/"annot"/f"PD_{g}"/f"PD_{g}.")]
    if kind == "BM": return [str(E/"annot"/f"BM_{g}"/f"BM_{g}.")]
    if kind == "J1C": return [str(E/"annot"/f"BM_{g}"/f"BM_{g}."), d3[0]]
    if kind == "J1D": return [str(E/"annot"/f"BM_{g}"/f"BM_{g}."), d3[1]]
    if kind == "CONT": return [str(O3/"annot"/f"C3_{g}"/f"C3_{g}.")]
    if kind == "BIN": return [str(O3/"annot"/f"B12_{g}"/f"B12_{g}.")]
    if kind == "EAb": return d3 + [str(W/"annot"/f"AB_{g}"/f"AB_{g}.")]
    raise ValueError(model)
def jobs_bench():
    j = [(t, f"M1_{g}") for t in KEEP8 for g in LIN] + [(t, f"M2_{g}") for t in B14 for g in LIN] + [(t, f"T3all_{g}") for t in NEW6 for g in LIN]
    j += [(t, f"{k}_{OWN[t]}") for t in NEW6 for k in ("M3", "PD")] + [(t, f"{k}_{OWN[t]}") for t in B14 for k in ("BM", "J1C")] + [(t, f"J1D_{OWN[t]}") for t in DNASE_TRAITS]
    j += [("N_AD", "CONT_myeloid"), ("N_AD", "BIN_myeloid")]
    return [(t, m, SS[t]) for t, m in j]
def jobs_panel():
    """106개 GWAS: 농축된 GWAS–조직 묶음 쌍(조직 묶음 모형 단측 P < 0.05/12) 가운데 다섯 묶음의 쌍. 재추출 20회(없는 것만)와 세포형 일치 모형."""
    u = pd.read_csv(EXP.parent/"e20_gwas_panel"/"results"/"u02_lineage.csv", encoding="utf-8-sig"); e = u[(u.p1 < .05 / 12) & (u.tau > 0) & u.lineage.isin(FIVE)]
    j = []
    for r in e.itertuples():
        j += [(r.id, f"T3r{s:02d}_{r.lineage}") for s in range(1, 21) if not (U/"runs"/f"{r.id}__T3r{s:02d}_{r.lineage}.npz").exists()]
        j.append((r.id, f"M3_{r.lineage}"))
    return [(t, m, U/"sumstats"/f"{t}.sumstats.gz") for t, m in j]
def jobs_actbin(): return [(t, f"EAb_{OWN[t]}", SS[t]) for t in B14]
def jobs_panelpd():
    """106개 GWAS: 다섯 묶음의 농축된 쌍에 근위·원위 모형 (PROTOCOL_exp3w_추가)."""
    u = pd.read_csv(EXP.parent/"e20_gwas_panel"/"results"/"u02_lineage.csv", encoding="utf-8-sig"); e = u[(u.p1 < .05 / 12) & (u.tau > 0) & u.lineage.isin(FIVE)]
    return [(r.id, f"PD_{r.lineage}", U/"sumstats"/f"{r.id}.sumstats.gz") for r in e.itertuples()]
JOBS = {"bench": jobs_bench, "panel": jobs_panel, "actbin": jobs_actbin, "panelpd": jobs_panelpd}
if STAGE == "jobs":
    for k, f in JOBS.items():
        j = f(); P(k, len(j), "개 | 모형 종류:", pd.Series([m.split("_")[0][:3] for _, m, _ in j]).value_counts().to_dict(), "| 이미 있는 것:", sum((RUN/f"{t}__{m}.npz").exists() for t, m, _ in j))
    sys.exit(0)
# --- ldsc 불러오기와 패치 (3u u01과 같다) ---
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
    T = W/"test"; T.mkdir(parents=True, exist_ok=True); worst = 0
    for name, ss, model, old in [("3n 키(더 큰 GWAS) M1_stromal", SS["N_HEIGHT"], "M1_stromal", N3/"runs"/"HEIGHT__M1_stromal.npz"),
                                 ("3j 골밀도 M3_stromal", SS["BMD"], "M3_stromal", J/"runs"/"BMD__M3_stromal.npz"),
                                 ("3k 교육 연수 PD_CNS_tissue", SS["EA"], "PD_CNS_tissue", K3/"runs"/"EA__PD_CNS_tissue.npz")]:
        a = run(ss, model, T/("CHECK_" + re.sub(r"\W+", "_", name))); b = np.load(old); d1, d2 = rel(a["coef"], b["coef"]), rel(a["part_delete"], b["part_delete"]); worst = max(worst, d1, d2)
        P(f"확인: {name} — 계수 상대 차이 {d1:.1e}, delete 값 {d2:.1e}")
    assert worst < 1e-9, "재현 실패 — 멈춘다"; P("통과"); (RES/"w01_check.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit(0)
if STAGE == "fit":
    which, k, K = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]); jobs = JOBS[which]()[k::K]; t0 = time.time(); nf = 0
    for i, (t, m, ss) in enumerate(jobs):
        out = RUN/f"{t}__{m}"; fail = Path(f"{out}.FAILED.txt")
        if fail.exists(): nf += 1; continue
        try: run(ss, m, out)
        except EmptyBlock as e: fail.write_text(str(e), encoding="utf-8"); nf += 1
        if (i + 1) % 10 == 0 or i + 1 == len(jobs):
            (W/f"STATUS_w01_{which}_{k}of{K}.txt").write_text(f"[{i + 1}/{len(jobs)}] {t} {m} ({(time.time() - t0) / 60:.1f}분, 실패 {nf})", encoding="utf-8")
    P(f"fit {which} 조각 {k}/{K}: {len(jobs)}개 끝, 실패 {nf}개, {(time.time() - t0) / 60:.1f}분"); (RES/f"w01_fit_{which}_{k}of{K}.txt").write_text("\n".join(log), encoding="utf-8")
