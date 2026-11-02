"""실험 3n n02: 최신 GWAS 다섯 개의 S-LDSC (PROTOCOL_exp3n_2026-10-21, 2026-10-21 고정). 순차 실행, 한 프로세스.
 실행기는 3j j02와 같다: 원래 ldsc(ldsc_bv)의 main()을 같은 회귀 인자로 부르고 jackknife 블록 경계만 3j 정렬 경계로 바꾼다.
 모형(주석은 모두 기존 것, 새 주석 없음; 계통은 원래 지도와 같다):
   T3_01–T3_20 (3f 재추출) + T3all (3d 전체 트랙)  — 층위 지도
   M1_<계통> 12개 (3c 계통 주석 L5_<계통>)            — 계통 지도 (기술)
   EA (3d 효과 3층위 + 3l 활성도 3층위)               — 활성도 조건화 (기술)
 단계 (python n02_sldsc.py check | main):
   check = 원래 키 GWAS(PASS_Height1)로 이 실행기를 돌려 기존 출력과 같아야 한다: 3j HEIGHT__T3all, 3l HEIGHT__EA (계수·delete 값 상대 차이 < 1e-9),
           3c HEIGHT M1 stromal (계수, .results). 어긋나면 멈춘다.
   main  = 5형질 × 34 = 170회. 재실행 안전: npz가 있으면 건너뛴다.
출력: exp3n/runs[/test]/<형질>__<모형>.npz (names, coef, part_delete, tot, h2_log, n_snp, sep), .log; results/n02_<단계>.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import sys, re, io, time, contextlib
import numpy as np, pandas as pd
from pathlib import Path
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); B = X/"GRCh38"; C = Path(f"{DATA}/exp3c"); Y = Path(f"{DATA}/exp3d")
F = Path(f"{DATA}/exp3f"); J = Path(f"{DATA}/exp3j"); L = Path(f"{DATA}/exp3l"); N = Path(f"{DATA}/exp3n")
BV = Path(f"{DATA}/tools/ldsc_bv")
RUN = N/"test" if STAGE == "check" else N/"runs"; RUN.mkdir(parents=True, exist_ok=True)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
TRAITS = {t: N/"sumstats"/f"{t}.sumstats.gz" for t in ["HEIGHT", "BMI", "SCZ", "BIP", "LDL"]}
OLD_HEIGHT = Path(f"{DATA}/exp3b/sumstats/PASS_Height1.sumstats")
OWN = {"HEIGHT": "stromal", "BMI": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "LDL": "liver"}
LAY = ["CAGE", "DNASE", "H3K27ac"]
LIN12 = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
def refs(t, model):
    g = OWN[t]; eff = [str(Y/"annot"/f"L5_{g}_{m}"/f"L5_{g}_{m}.") for m in LAY]
    if model == "T3all": return eff
    if model.startswith("T3_"): s = model[3:]; return [str(F/"annot"/f"T3_{s}_{g}"/f"T3_{s}_{g}.")]
    if model.startswith("M1_"): h = model[3:]; return [str(C/"annot"/f"L5_{h}"/f"L5_{h}.")]
    if model == "EA": return eff + [str(L/"annot"/f"A5_{g}"/f"A5_{g}.")]
    raise ValueError(model)
MODELS = [f"T3_{s:02d}" for s in range(1, 21)] + ["T3all"] + [f"M1_{h}" for h in LIN12] + ["EA"]
MAIN = [(t, m) for t in TRAITS for m in MODELS]; assert len(MAIN) == 170 and len(set(MAIN)) == 170
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
    cm = max(0.001 * ss.N.max(), 80) if len(ref_ld_cnames) > 1 else None   # estimate_h2와 같은 χ² 상한 (분할 모형)
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
def run(t, model, ss, g=None):
    out = RUN/f"{t}__{model}"
    if Path(f"{out}.npz").exists(): return dict(np.load(f"{out}.npz", allow_pickle=False))
    ST.update(sep=None, names=None, hsq=None)
    if g: OWN[t] = g
    sys.argv = ["ldsc", "--h2", str(ss), "--ref-ld-chr", ",".join([str(B/"baselineLD_v2.2"/"baselineLD.")] + refs(t, model)), "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."), "--out", str(out)]
    with contextlib.redirect_stdout(io.StringIO()): LD.main()
    for l in LOGS: l.log_fh.close()
    LOGS.clear(); h = ST["hsq"]
    h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log", encoding="utf-8").read()).group(1))
    r = dict(names=np.array(ST["names"]), coef=np.asarray(h.coef).ravel(), part_delete=np.asarray(h.part_delete_values), tot=float(np.asarray(h.tot).ravel()[0]),
             h2_log=h2, n_snp=int(ST["sep"][-1]), sep=np.asarray(ST["sep"]))
    assert r["part_delete"].shape == (200, len(r["coef"])), r["part_delete"].shape
    np.savez(f"{out}.npz", **r); return r
def rel(a, b): a, b = np.asarray(a, float), np.asarray(b, float); return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300))
t0 = time.time()
if STAGE == "check":
    OWN["OLDH"] = "stromal"
    for model, ref in [("T3all", J/"runs"/"HEIGHT__T3all.npz"), ("EA", L/"runs"/"HEIGHT__EA.npz")]:
        a = run("OLDH", model, OLD_HEIGHT); z = np.load(ref)
        r1, r2 = rel(a["coef"], z["coef"]), rel(a["part_delete"], z["part_delete"])
        P(f"  원래 키 {model}: 기존 실행과 계수 차이 {r1:.1e}, delete 값 차이 {r2:.1e} | 회귀 SNP {a['n_snp']:,} | h2 {a['h2_log']} | {(time.time() - t0) / 60:.1f}분")
        assert r1 < 1e-9 and r2 < 1e-9, "재현 실패 — 멈춘다"
    a = run("OLDH", "M1_stromal", OLD_HEIGHT); R = pd.read_csv(C/"sldsc"/"HEIGHT__M1__stromal.results", sep="\t")
    r3 = rel(a["coef"], R.Coefficient); P(f"  원래 키 M1 stromal: 3c .results와 계수 차이 {r3:.1e}"); assert r3 < 1e-9, "재현 실패 — 멈춘다"
    P("확인 통과")
else:
    for i, (t, model) in enumerate(MAIN):
        run(t, model, TRAITS[t]); el = (time.time() - t0) / 60
        msg = f"[{i + 1}/{len(MAIN)}] {t} {model} 완료 ({el:.0f}분, 남은 약 {el / (i + 1) * (len(MAIN) - i - 1):.0f}분)"
        P("  " + msg); (N/"PIPELINE_STATUS.txt").write_text(msg, encoding="utf-8")
(EXP/"results"/f"n02_{STAGE}.txt").write_text("\n".join(log), encoding="utf-8")
