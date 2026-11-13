"""실험 3t t08: 형질 × 계통 전체 행렬 (PROTOCOL_exp3t_2026-11-10 T8).
 전체 트랙 세 층위 결합 모형(L5_<계통>_{CAGE, DNASE, H3K27ac})을 14형질 × 12계통에 적합한다. 정렬 블록, 실행기는 3j j02와 같다.
 기존 47개(exp3j/runs의 T3all, T3all_<계통>)는 읽기만 하고, 나머지 121개를 exp3t/runs에 쓴다.
 단계 (python t08_matrix.py fit | test):
   fit  = 확인(HEIGHT adipose를 다시 돌려 exp3j 값과 계수 상대 차이 < 1e-9) 뒤 121회. 재실행 안전.
   test = 14 × 12 행렬(CAGE − DNase, φD), 계통별 Q, 계통별 키 대 다른 형질 비교.
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import sys, re, io, time, contextlib
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
STAGE = sys.argv[1]
DATA = Path(DATA); EXP = Path(__file__).resolve().parents[1]; RES = EXP/"results"; RES.mkdir(exist_ok=True)
X = DATA/"exp3"; B = X/"GRCh38"; Y = DATA/"exp3d"; F = DATA/"exp3f"; J = DATA/"exp3j"; T = DATA/"exp3t"
RUN = T/"runs"; RUN.mkdir(parents=True, exist_ok=True); BV = DATA/"tools"/"ldsc_bv"; SS = DATA/"exp3b"/"sumstats"
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
TRAITS = {"AD": X/"bellenguez_noAPOE.sumstats.gz", "MONO": DATA/"exp3b"/"monocyte_count.sumstats.gz", "EA": SS/"PASS_Years_of_Education2.sumstats",
          "SCZ": SS/"PASS_Schizophrenia.sumstats", "BIP": SS/"PASS_Bipolar_Disorder.sumstats", "NEURO": SS/"PASS_Neuroticism.sumstats", "BMI": SS/"PASS_BMI1.sumstats",
          "IBD": SS/"PASS_IBD.sumstats", "RA": SS/"PASS_Rheumatoid_Arthritis.sumstats", "LDL": SS/"PASS_LDL.sumstats", "HEIGHT": SS/"PASS_Height1.sumstats",
          "BMD": F/"sumstats"/"BMD.sumstats.gz", "NEUT": F/"sumstats"/"NEUT.sumstats.gz", "TC": F/"sumstats"/"TC.sumstats.gz"}
OWN = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
       "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal", "BMD": "stromal", "NEUT": "myeloid", "TC": "liver"}
MAP = ["AD", "MONO", "EA", "SCZ", "BIP", "NEURO", "BMI", "IBD", "RA", "LDL", "HEIGHT"]; VAL = ["BMD", "NEUT", "TC"]; ALL = MAP + VAL
LIN = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
LAY = ["CAGE", "DNASE", "H3K27ac"]
M_TOT = sum(float(open(B/"baselineLD_v2.2"/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
strip = lambda n: re.sub(r"L2_\d+$", "", n)
def old_path(t, g):
    p = J/"runs"/(f"{t}__T3all.npz" if OWN[t] == g else f"{t}__T3all_{g}.npz"); return p if p.exists() else None
def path(t, g): return old_path(t, g) or RUN/f"{t}__T3all_{g}.npz"
def jcov(Z): Bn = Z.shape[0]; z = Z - Z.mean(0); return (Bn - 1) / Bn * z.T @ z
def jse(v): return float(np.sqrt(jcov(np.asarray(v).reshape(-1, 1))[0, 0]))
if STAGE == "fit":
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
    def separators(snps):
        k = KEY.reindex(snps).to_numpy(); assert not np.isnan(k).any() and (np.diff(k) >= 0).all()
        c = np.bincount(np.searchsorted(BOUNDS, k.astype(np.int64), side="right"), minlength=NB); assert len(c) == NB and (c > 0).all()
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
    def run(t, g, out):
        if Path(f"{out}.npz").exists(): return dict(np.load(f"{out}.npz", allow_pickle=False))
        ST.update(sep=None, names=None, hsq=None)
        sys.argv = ["ldsc", "--h2", str(TRAITS[t]), "--ref-ld-chr", ",".join([str(B/"baselineLD_v2.2"/"baselineLD.")] + [str(Y/"annot"/f"L5_{g}_{m}"/f"L5_{g}_{m}.") for m in LAY]),
                    "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."), "--out", str(out)]
        with contextlib.redirect_stdout(io.StringIO()): LD.main()
        for l in LOGS: l.log_fh.close()
        LOGS.clear(); h = ST["hsq"]; h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log", encoding="utf-8").read()).group(1))
        r = dict(names=np.array(ST["names"]), coef=np.asarray(h.coef).ravel(), part_delete=np.asarray(h.part_delete_values), tot=float(np.asarray(h.tot).ravel()[0]), h2_log=h2, n_snp=int(ST["sep"][-1]), sep=np.asarray(ST["sep"]))
        assert r["part_delete"].shape == (200, len(r["coef"])); np.savez(f"{out}.npz", **r); return r
    t0 = time.time(); a = run("HEIGHT", "adipose", RUN/"CHECK_HEIGHT__T3all_adipose"); b = np.load(J/"runs"/"HEIGHT__T3all_adipose.npz")
    d = float(np.max(np.abs(a["coef"] - b["coef"])) / np.max(np.abs(b["coef"]))); d2 = float(np.max(np.abs(a["part_delete"] - b["part_delete"])) / np.max(np.abs(b["part_delete"])))
    P(f"확인: HEIGHT adipose 다시 적합 — 계수 상대 차이 {d:.1e}, delete 값 {d2:.1e}"); assert d < 1e-9 and d2 < 1e-9, "재현 실패 — 멈춘다"
    jobs = [(t, g) for g in LIN for t in ALL if old_path(t, g) is None]; P(f"새 적합 {len(jobs)}회 (기존 {14 * 12 - len(jobs)}개는 exp3j에서 읽음)")
    for i, (t, g) in enumerate(jobs):
        run(t, g, RUN/f"{t}__T3all_{g}"); msg = f"[{i + 1}/{len(jobs)}] {t} {g} ({(time.time() - t0) / 60:.1f}분)"; (T/"PIPELINE_STATUS_t08.txt").write_text(msg, encoding="utf-8")
        if (i + 1) % 20 == 0 or i + 1 == len(jobs): P("  " + msg)
    (RES/"t08_fit.txt").write_text("\n".join(log), encoding="utf-8")
else:
    SD = pd.read_csv(Y/"annot_sd_common.csv", index_col=0).sd_common; ang = lambda y, x: np.degrees(np.arctan2(y, x))
    def load(t, g):
        z = np.load(path(t, g)); h2 = float(z["h2_log"]); names, coef, D = [strip(n) for n in z["names"]], z["coef"], z["part_delete"]; res = {}
        for k, n in enumerate(names):
            if n in SD.index: f = SD[n] * M_TOT / h2; res[[m for m in LAY if n.endswith("_" + m)][0]] = (coef[k] * f, D[:, k] * f)
        assert sorted(res) == sorted(LAY); return res
    L = {(t, g): load(t, g) for t in ALL for g in LIN}; rows = []
    for (t, g), l in L.items():
        C, D, H = l["CAGE"], l["DNASE"], l["H3K27ac"]
        rows.append(dict(trait=t, lineage=g, own=OWN[t] == g, cage=C[0], cage_se=jse(C[1]), dnase=D[0], dnase_se=jse(D[1]), h3k27ac=H[0], h3k27ac_se=jse(H[1]),
                         c1=C[0] - D[0], c1_se=jse(C[1] - D[1]), phiD=float(ang(D[0], C[0])), phiD_se=jse(ang(D[1], C[1]))))
    M = pd.DataFrame(rows); M.to_csv(RES/"t08_matrix.csv", index=False)
    def Qtest(g, ts):
        d = np.concatenate([[L[t, g]["CAGE"][0] - L[t, g]["DNASE"][0], L[t, g]["CAGE"][0] - L[t, g]["H3K27ac"][0]] for t in ts])
        Z = np.hstack([np.column_stack([L[t, g]["CAGE"][1] - L[t, g]["DNASE"][1], L[t, g]["CAGE"][1] - L[t, g]["H3K27ac"][1]]) for t in ts]); Si = np.linalg.inv(jcov(Z)); K = len(ts)
        Xm = np.tile(np.eye(2), (K, 1)); mu = np.linalg.solve(Xm.T @ Si @ Xm, Xm.T @ Si @ d); r = d - Xm @ mu; Q = float(r @ Si @ r); df = 2 * (K - 1); nu = 199
        return Q, df, float(stats.chi2.sf(Q, df)), float(stats.f.sf((nu - df + 1) / (nu * df) * Q, df, nu - df + 1))
    P("=== T8 형질 × 계통 행렬 (전체 트랙, 정렬 블록) ===")
    P("\nCAGE − DNase (형질 행 × 계통 열; * = 원고에서 그 형질을 분석한 계통)")
    piv = M.pivot(index="trait", columns="lineage", values="c1").loc[ALL, LIN]; P(piv.round(2).to_string())
    P("\n계통마다 지도 형질 11개의 Q (차이 두 개, 자유도 20)")
    for g in LIN:
        Q, df, p, pF = Qtest(g, MAP); P(f"  {g:22s} Q = {Q:5.1f}, P = {p:.2g}, F 근사 P = {pF:.2g}")
    P("\nstromal 계통의 주석에서 (키·BMD와 같은 주석): τ* CAGE, DNase, H3K27ac (s.e.) | CAGE − DNase | φD")
    for t in ALL:
        r = M[(M.trait == t) & (M.lineage == "stromal")].iloc[0]
        P(f"  {t:7s} {r.cage:+.2f} ({r.cage_se:.2f}), {r.dnase:+.2f} ({r.dnase_se:.2f}), {r.h3k27ac:+.2f} ({r.h3k27ac_se:.2f}) | {r.c1:+.2f} ({r.c1_se:.2f}) | {r.phiD:+.0f}° ({r.phiD_se:.0f})")
    P("\n계통마다: 골격 형질(키, BMD) 대 나머지 지도 형질 열 개의 CAGE − DNase 차이 (나머지 − 골격; IVW, 공동 jackknife, 양측 P)")
    ests, dvs = [], []
    for g in LIN:
        others = [t for t in MAP if t != "HEIGHT"]; sk = ["HEIGHT", "BMD"]; c1 = {t: (L[t, g]["CAGE"][0] - L[t, g]["DNASE"][0], L[t, g]["CAGE"][1] - L[t, g]["DNASE"][1]) for t in others + sk}
        def ivw(ts): w = np.array([1 / jse(c1[t][1]) ** 2 for t in ts]); w /= w.sum(); return float(sum(wi * c1[t][0] for wi, t in zip(w, ts))), sum(wi * c1[t][1] for wi, t in zip(w, ts))
        eo, do = ivw(others); es, ds = ivw(sk); e = eo - es; dv = do - ds; s = jse(dv); ests.append(e); dvs.append(dv)
        P(f"  {g:22s} 나머지 {eo:+.2f} | 골격 {es:+.2f} | 차이 {e:+.2f} ± {s:.2f}, P = {2 * stats.norm.sf(abs(e / s)):.2g}")
    w = np.array([1 / jse(d) ** 2 for d in dvs]); w /= w.sum(); e = float(w @ np.array(ests)); s = jse(sum(wi * d for wi, d in zip(w, dvs)))
    P(f"  12계통 합침: {e:+.2f} ± {s:.2f}, P = {2 * stats.norm.sf(abs(e / s)):.2g}; {sum(x > 0 for x in ests)}/12 계통에서 양수")
    (RES/"t08_tests.txt").write_text("\n".join(log), encoding="utf-8")
