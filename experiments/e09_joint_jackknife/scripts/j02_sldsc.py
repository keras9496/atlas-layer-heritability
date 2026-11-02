"""실험 3j j02: 블록 정렬 S-LDSC (PROTOCOL_exp3j_2026-10-14). 순차 실행, 한 프로세스.
 원래 ldsc(ldsc_bv, 3f·3g·3h·3i와 같은 판)의 main()을 같은 회귀 인자로 부르고, jackknife 블록 경계만 j01 경계로 바꾼다
 (분할 유전율의 jackknife는 LstsqJackknifeFast 하나에서 블록을 정하므로 여기에 separators를 넘긴다).
 --overlap-annot, --frqfile-chr, --print-coefficients는 .results 출력(비율·농축)용이라 뺀다. 계수와 delete 값은 회귀 객체에서 직접 저장한다.
 기준선 LD 점수와 회귀 가중치는 한 번만 읽는다.
 단계 (python j02_sldsc.py check | dry | main):
   check = 기존 실행 재현. 3f·3g·3h·3i에서 하나씩, 블록을 바꾸지 않은 실행이 기존 출력(.results 계수 전체, .part_delete 전체)과 같아야 한다(상대 차이 < 1e-9).
           같은 작업을 정렬 블록으로 돌려 점 추정치가 같아야 한다(상대 차이 < 1e-9). 어긋나면 멈춘다.
   dry   = 가짜 형질 3개(AD Z 섞음, 결과 아님) × 골수계 T3all, 정렬 블록 — j03 형식 확인용
   main  = 383회: 3f T3_01–20 × 14형질 (280) + T3all 14, 3h 33, 3i 14, 3g 42 (이 순서)
 재실행 안전: npz가 있으면 건너뛴다.
출력: exp3j/runs[/test]/<형질>__<모형>.npz (names, coef, part_delete, tot, h2_log, n_snp, sep), .log; results/j02_<단계>.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import sys, os, re, io, time, contextlib, glob
import numpy as np, pandas as pd
from pathlib import Path
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); B = X/"GRCh38"; Y = Path(f"{DATA}/exp3d"); F = Path(f"{DATA}/exp3f")
G = Path(f"{DATA}/exp3g"); H = Path(f"{DATA}/exp3h"); I = Path(f"{DATA}/exp3i"); J = Path(f"{DATA}/exp3j")
BV = Path(f"{DATA}/tools/ldsc_bv"); SS = Path(f"{DATA}/exp3b/sumstats")
RUN = J/"test" if STAGE in ("check", "dry") else J/"runs"; RUN.mkdir(parents=True, exist_ok=True)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
TRAITS = {"AD": X/"bellenguez_noAPOE.sumstats.gz", "MONO": Path(f"{DATA}/exp3b/monocyte_count.sumstats.gz"), "EA": SS/"PASS_Years_of_Education2.sumstats",
          "SCZ": SS/"PASS_Schizophrenia.sumstats", "BIP": SS/"PASS_Bipolar_Disorder.sumstats", "NEURO": SS/"PASS_Neuroticism.sumstats", "BMI": SS/"PASS_BMI1.sumstats",
          "IBD": SS/"PASS_IBD.sumstats", "RA": SS/"PASS_Rheumatoid_Arthritis.sumstats", "LDL": SS/"PASS_LDL.sumstats", "HEIGHT": SS/"PASS_Height1.sumstats",
          "BMD": F/"sumstats"/"BMD.sumstats.gz", "NEUT": F/"sumstats"/"NEUT.sumstats.gz", "TC": F/"sumstats"/"TC.sumstats.gz"}   # 3f v02와 같다
OWN = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
       "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal", "BMD": "stromal", "NEUT": "myeloid", "TC": "liver"}
LAY = ["CAGE", "DNASE", "H3K27ac"]
def refs(t, model):
    g = OWN.get(t, "myeloid")
    if model == "T3all": return [str(Y/"annot"/f"L5_{g}_{m}"/f"L5_{g}_{m}.") for m in LAY]
    if model.startswith("T3all_"): h = model[6:]; return [str(Y/"annot"/f"L5_{h}_{m}"/f"L5_{h}_{m}.") for m in LAY]
    if model.startswith("T3_"): s = model[3:]; return [str(F/"annot"/f"T3_{s}_{g}"/f"T3_{s}_{g}.")]
    if model.startswith("G3_"): q = model[3:]; return [str(G/"annot"/f"G3_{q}_{g}"/f"G3_{q}_{g}.")]
    if model.startswith("M3_"): return [str(I/"annot"/model/f"{model}.")]
    raise ValueError(model)
H3 = sorted((Path(f).name.split("__")[0], Path(f).name.split("__")[1][:-8]) for f in glob.glob(str(H/"sldsc"/"*__T3all_*.results")))
assert len(H3) == 33, len(H3)   # 3h 프로토콜의 작업 수
MAIN = ([(t, f"T3_{s:02d}") for t in TRAITS for s in range(1, 21)] + [(t, "T3all") for t in TRAITS] + H3
        + [(t, f"M3_{OWN[t]}") for t in TRAITS] + [(t, f"G3_{q}") for t in TRAITS for q in ["L0p1", "L0p5", "L1p0"]])
assert len(MAIN) == 383 and len(set(MAIN)) == 383
# --- ldsc 불러오기와 패치 ---
sys.path.insert(0, str(BV))
import ldscore.ldsc as LD
from ldscore import sumstats as S, jackknife as jk, parse as ps
BL = np.load(J/"blocks.npz"); BOUNDS = BL["bounds"]; NB = len(BOUNDS) + 1
KEY = pd.read_parquet(J/"snp_key.parquet").drop_duplicates("SNP").set_index("SNP").key
ST = dict(mode="equal", sep=None, names=None, hsq=None)
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
    ST["sep"] = separators(snps) if ST["mode"] == "aligned" else np.floor(np.linspace(0, len(snps), 201)).astype(int)
    return out
S._read_ld_sumstats = rls
_LJF = jk.LstsqJackknifeFast
class LJF(_LJF):
    def __init__(self, x, y, n_blocks=None, separators=None):
        if ST["mode"] == "aligned":
            separators = list(ST["sep"]); assert separators[-1] == x.shape[0]; n_blocks = None
        _LJF.__init__(self, x, y, n_blocks, separators)
jk.LstsqJackknifeFast = LJF
_eh = S.estimate_h2
def eh(args, lg): r = _eh(args, lg); ST["hsq"] = r; return r
S.estimate_h2 = eh
LOGS = []; _li = LD.Logger.__init__
def li(self, fh): _li(self, fh); LOGS.append(self)
LD.Logger.__init__ = li
def run(t, model, mode, ss=None):
    out = RUN/f"{t}__{model}" if mode == "aligned" else RUN/f"{t}__{model}__equal"
    if Path(f"{out}.npz").exists(): return dict(np.load(f"{out}.npz", allow_pickle=False))
    ST.update(mode=mode, sep=None, names=None, hsq=None)
    sys.argv = ["ldsc", "--h2", str(ss or TRAITS[t]), "--ref-ld-chr", ",".join([str(B/"baselineLD_v2.2"/"baselineLD.")] + refs(t, model)),
                "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."), "--out", str(out)]
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
    old = {"3f": (F/"sldsc", "AD", "T3_01"), "3g": (G/"sldsc", "BIP", "G3_L1p0"), "3h": (H/"sldsc", H3[0][0], H3[0][1]), "3i": (I/"sldsc", "HEIGHT", "M3_stromal")}
    for k, (d, t, model) in old.items():
        R = pd.read_csv(d/f"{t}__{model}.results", sep="\t"); D0 = np.loadtxt(d/f"{t}__{model}.part_delete")
        e = run(t, model, "equal"); a = run(t, model, "aligned")
        r1, r2, r3 = rel(e["coef"], R.Coefficient), rel(e["part_delete"], D0), rel(a["coef"], e["coef"])
        nb = np.diff(a["sep"])
        P(f"  {k} {t} {model}: 등분 블록 재현 — 계수 {r1:.1e}, delete 값 {r2:.1e} | 정렬 블록 점 추정치 차이 {r3:.1e} | 회귀 SNP {a['n_snp']:,}, 정렬 블록 SNP {nb.min():,}–{nb.max():,}"
          f" | h2 {e['h2_log']} | {(time.time() - t0) / 60:.1f}분")
        assert r1 < 1e-9 and r2 < 1e-9 and r3 < 1e-9, "재현 실패 — 멈춘다"
    P("확인 통과")
elif STAGE == "dry":
    d0 = pd.read_csv(TRAITS["AD"], sep="\t")
    for i in range(1, 4):
        fk = J/"test"/f"FAKE{i}.sumstats.gz"
        if not fk.exists(): d = d0.copy(); d["Z"] = np.random.default_rng(20261014 + i).permutation(d.Z.to_numpy()); d.to_csv(fk, sep="\t", index=False, compression="gzip")
        r = run(f"FAKE{i}", "T3all", "aligned", ss=fk); P(f"  FAKE{i} T3all 완료 (h2 {r['h2_log']}, SNP {r['n_snp']:,}; {(time.time() - t0) / 60:.1f}분)")
else:
    for i, (t, model) in enumerate(MAIN):
        run(t, model, "aligned"); el = (time.time() - t0) / 60
        msg = f"[{i + 1}/{len(MAIN)}] {t} {model} 완료 ({el:.0f}분, 남은 약 {el / (i + 1) * (len(MAIN) - i - 1):.0f}분)"
        P("  " + msg); (J/"PIPELINE_STATUS.txt").write_text(msg, encoding="utf-8")
(EXP/"results"/f"j02_{STAGE}.txt").write_text("\n".join(log), encoding="utf-8")
