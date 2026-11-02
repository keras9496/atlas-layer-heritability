"""실험 3r r02: 새 AD GWAS(N_AD, GCST90704647 proxy 제외판)의 S-LDSC (PROTOCOL_exp3r_2026-10-25, 2026-10-25 고정). 순차 실행, 한 프로세스.
 실행기는 3j j02 / 3n n02 / 3p p02 / 3q q02와 같다: 원래 ldsc(ldsc_bv)의 main()을 같은 회귀 인자로 부르고 jackknife 블록 경계만 3j 정렬 경계로 바꾼다.
 모형(주석은 모두 기존 것, 계통 = 골수계):
   T3_01–T3_20 (3f 재추출) + T3all (3d 전체 트랙)           — 층위 지도
   M1_<계통> 12개 (3c 계통 주석 L5_<계통>)                  — 계통 지도
   EA, EAM (3l: 효과 + 이진 활성도 A5 [+ 실측 BM])          — 활성도 조건화
   EAc, S, SA (3p: 연속 활성도 Ac3 / 층화 효과 S5 / S5 + A5)
   HL10, HL20 (3q: 활성도 상위 10 %·20 % 안 활성도 맞춤 효과 상·하 반분)
 단계 (python r02_sldsc.py check | main):
   check = 원래 AD GWAS(Bellenguez, APOE 제외)로 T3all·EA·EAc·HL10을 이 실행기로 다시 돌려 기존 npz와 계수·delete 값 상대 차이 < 1e-9이어야 한다. 어긋나면 멈춘다.
   main  = N_AD × 40회. 재실행 안전: npz가 있으면 건너뛴다.
출력: exp3r/runs[/test]/<형질>__<모형>.npz (names, coef, part_delete, tot, h2_log, n_snp, sep), .log; results/r02_<단계>.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import sys, re, io, time, contextlib
import numpy as np, pandas as pd
from pathlib import Path
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); B = X/"GRCh38"; C = Path(f"{DATA}/exp3c"); Y = Path(f"{DATA}/exp3d")
F = Path(f"{DATA}/exp3f"); J = Path(f"{DATA}/exp3j"); L = Path(f"{DATA}/exp3l"); W = Path(f"{DATA}/exp3e")
O = Path(f"{DATA}/exp3p"); Q = Path(f"{DATA}/exp3q"); R = Path(f"{DATA}/exp3r")
BV = Path(f"{DATA}/tools/ldsc_bv")
RUN = R/"test" if STAGE == "check" else R/"runs"; RUN.mkdir(parents=True, exist_ok=True)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
TRAITS = {"N_AD": R/"sumstats"/"N_AD.sumstats.gz", "AD": X/"bellenguez_noAPOE.sumstats.gz"}
G = "myeloid"; LAY = ["CAGE", "DNASE", "H3K27ac"]
LIN12 = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
def refs(model):
    g = G; eff = [str(Y/"annot"/f"L5_{g}_{m}"/f"L5_{g}_{m}.") for m in LAY]
    a5 = str(L/"annot"/f"A5_{g}"/f"A5_{g}."); s5 = str(O/"annot"/f"S5_{g}"/f"S5_{g}.")
    if model == "T3all": return eff
    if model.startswith("T3_"): return [str(F/"annot"/f"T3_{model[3:]}_{g}"/f"T3_{model[3:]}_{g}.")]
    if model.startswith("M1_"): h = model[3:]; return [str(C/"annot"/f"L5_{h}"/f"L5_{h}.")]
    if model.startswith("HL"): k = model[2:]; return [str(Q/"annot"/f"HL{k}_{g}"/f"HL{k}_{g}.")]
    return {"EA": eff + [a5], "EAM": eff + [a5, str(W/"annot"/f"BM_{g}"/f"BM_{g}.")], "EAc": eff + [str(O/"annot"/f"Ac3_{g}"/f"Ac3_{g}.")],
            "S": [s5], "SA": [s5, a5]}[model]
def nadd(model):
    if model.startswith("T3"): return 3
    if model.startswith("M1_"): return 1
    return {"EA": 6, "EAM": 9, "EAc": 6, "S": 3, "SA": 6, "HL10": 6, "HL20": 6}[model]
MODELS = [f"T3_{s:02d}" for s in range(1, 21)] + ["T3all"] + [f"M1_{h}" for h in LIN12] + ["EA", "EAM", "EAc", "S", "SA", "HL10", "HL20"]
MAIN = [("N_AD", m) for m in MODELS]; assert len(MAIN) == 40 and len(set(MAIN)) == 40
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
def run(t, model):
    out = RUN/f"{t}__{model}"
    if Path(f"{out}.npz").exists(): return dict(np.load(f"{out}.npz", allow_pickle=False))
    ST.update(sep=None, names=None, hsq=None)
    sys.argv = ["ldsc", "--h2", str(TRAITS[t]), "--ref-ld-chr", ",".join([str(B/"baselineLD_v2.2"/"baselineLD.")] + refs(model)), "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."), "--out", str(out)]
    with contextlib.redirect_stdout(io.StringIO()): LD.main()
    for l in LOGS: l.log_fh.close()
    LOGS.clear(); h = ST["hsq"]
    h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log", encoding="utf-8").read()).group(1))
    r = dict(names=np.array(ST["names"]), coef=np.asarray(h.coef).ravel(), part_delete=np.asarray(h.part_delete_values), tot=float(np.asarray(h.tot).ravel()[0]),
             h2_log=h2, n_snp=int(ST["sep"][-1]), sep=np.asarray(ST["sep"]))
    assert r["part_delete"].shape == (200, len(r["coef"])), r["part_delete"].shape
    assert len(r["coef"]) == 97 + nadd(model), (model, len(r["coef"]))
    np.savez(f"{out}.npz", **r); return r
def rel(a, b): a, b = np.asarray(a, float), np.asarray(b, float); return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300))
t0 = time.time()
if STAGE == "check":
    for model, ref in [("T3all", J/"runs"/"AD__T3all.npz"), ("EA", L/"runs"/"AD__EA.npz"), ("EAc", O/"runs"/"AD__EAc.npz"), ("HL10", Q/"runs"/"AD__HL10.npz")]:
        a = run("AD", model); z = np.load(ref)
        r1, r2 = rel(a["coef"], z["coef"]), rel(a["part_delete"], z["part_delete"])
        P(f"  원래 AD {model}: 기존 실행과 계수 차이 {r1:.1e}, delete 값 차이 {r2:.1e} | 회귀 SNP {a['n_snp']:,} | {(time.time() - t0) / 60:.1f}분")
        assert r1 < 1e-9 and r2 < 1e-9, "재현 실패 — 멈춘다"
    P("확인 통과")
else:
    for i, (t, model) in enumerate(MAIN):
        run(t, model); el = (time.time() - t0) / 60
        msg = f"[{i + 1}/{len(MAIN)}] {t} {model} 완료 ({el:.0f}분, 남은 약 {el / (i + 1) * (len(MAIN) - i - 1):.0f}분)"
        P("  " + msg); (R/"PIPELINE_STATUS.txt").write_text(msg, encoding="utf-8")
(EXP/"results"/f"r02_{STAGE}.txt").write_text("\n".join(log), encoding="utf-8")
