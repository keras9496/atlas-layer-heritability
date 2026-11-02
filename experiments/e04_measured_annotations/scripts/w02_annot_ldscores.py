"""실험 3e w02: 실측 계통 주석 (PROTOCOL_exp3e_2026-10-08) → 염색체별 LD 점수 (실험 3d y03과 같은 l2 호출).
 M_{g}_{m}: 용어마다 과반 규칙으로 활성 요소 → 계통 안 합집합 → 창(CAGE·DNASE ±250 bp, H3K27ac ±1,000 bp) 안의 참조 SNP = 1.   P_{g}_{m}: 창 없이 같은 규칙.
 FANTOM5: robust CAGE 피크(과반 라이브러리 TPM ≥ 1) ∪ 전사 인핸서(과반 라이브러리 used = 1).   ENCODE: 용어 실험의 과반에서 피크가 덮는 위치(실험 안 파일은 합집합).
 SNP: 1000G EUR hg38 참조 SNP 전체를 위치로 배정한다(희귀 SNP 포함).  좌표: BED 0 기반 [start, end), bim BP는 1 기반 → BP-1 ∈ [start-w, end+w).
 시험: W02_TEST=1 → chr22 만, exp3e/test/annot 에 쓴다.
출력: exp3e/annot/<이름>/<이름>.<chr>.{annot.gz,l2.ldscore.gz,l2.M,l2.M_5_50} (단일 30 + 묶음 BM_g·BMp_g 10), exp3e/measured_sd_common.csv, results/w02_annot.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import os, gzip, math, subprocess
import numpy as np, pandas as pd
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); Y = Path(f"{DATA}/exp3d"); W = Path(f"{DATA}/exp3e")
F5 = W/"raw"/"fantom5"; EN = W/"raw"/"encode"; PF = X/"GRCh38"/"plink_files"; LDSC = Path(f"{DATA}/tools/ldsc/ldsc.py")
TEST = os.environ.get("W02_TEST") == "1"; CHRS = [22] if TEST else list(range(1, 23)); AN = W/"test"/"annot" if TEST else W/"annot"; AN.mkdir(parents=True, exist_ok=True)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LING = ["myeloid", "CNS_tissue", "lymphoid", "liver", "stromal"]; MODS = ["CAGE", "DNASE", "H3K27ac"]; WIN = {"CAGE": 250, "DNASE": 250, "H3K27ac": 1000}
TS = pd.read_csv(W/"term_sources.csv"); LT = pd.read_csv(W/"fantom5_library_terms.csv").fillna(""); EF = pd.read_csv(W/"encode_files.csv")
# ---------- 참조 SNP (전역 좌표) ----------
BIM = {c: pd.read_csv(PF/f"1000G.EUR.hg38.{c}.bim", sep="\t", header=None, names=["CHR", "SNP", "CM", "BP", "a1", "a2"]) for c in CHRS}
OFF = {c: c * 10**9 for c in CHRS}
G = np.concatenate([OFF[c] + BIM[c].BP.to_numpy(np.int64) - 1 for c in CHRS]); ORD = np.argsort(G, kind="stable"); GS = G[ORD]; N = len(G)
P(f"참조 SNP {N:,}개 (염색체 {len(CHRS)}개)")
def cover(df, w):
    """BED 구간(chrom, start, end) ± w 가 덮는 참조 SNP (bool, bim 순서)."""
    ch = df.iloc[:, 0].astype(str).str.replace("chr", "", regex=False)
    ok = ch.isin([str(c) for c in CHRS]); d = df[ok.to_numpy()]; ch = ch[ok].astype(int).to_numpy()
    s = np.array([OFF[c] for c in ch], dtype=np.int64) + d.iloc[:, 1].to_numpy(np.int64) - w; e = np.array([OFF[c] for c in ch], dtype=np.int64) + d.iloc[:, 2].to_numpy(np.int64) + w
    o = np.argsort(s); s, e = s[o], e[o]
    if len(s) == 0: return np.zeros(N, bool)
    e = np.maximum.accumulate(e)                           # 정렬 후 누적 최대 끝 → 겹친 구간 병합과 같은 판정
    i = np.searchsorted(s, GS, side="right") - 1
    hit = np.zeros(N, bool); v = i >= 0; hit_s = np.zeros(len(GS), bool); hit_s[v] = GS[v] < e[i[v]]
    hit[ORD] = hit_s; return hit
def majority(n): return max(1, math.ceil(n / 2))
ANN = {}
# ---------- FANTOM5 CAGE ----------
SDRF = pd.read_excel(F5/"HumanSamples2.0.sdrf.xlsx"); SDRF["ff"] = "FF:" + SDRF["Source Name"].astype(str)
LT = LT.merge(SDRF[["ff", "Characteristics [Category]"]].rename(columns={"Characteristics [Category]": "cat"}), on="ff", how="left")
CAT = {"tissue": "tissues", "primary_cell": "primary cells"}   # 보완 1: 라이브러리 분류 = Atlas biosample_type (ENCODE 분류 규칙과 같게; 시간 경과·분획·세포주 제외)
CT = TS[TS["mod"] == "CAGE"].drop_duplicates(["term", "biosample_type"])
lib_of = {(t, b): LT[LT.atlas_terms.str.split(";").apply(lambda s: t in s) & (LT.cat == CAT[b])] for t, b in zip(CT.term, CT.biosample_type)}
need_cols = sorted(set(int(c) for v in lib_of.values() for c in v.col))
bed = pd.read_csv(F5/"hg38_fair+new_CAGE_peaks_phase1and2.bed.gz", sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "name"])
with gzip.open(F5/"hg38_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz", "rt") as fh:
    for k, line in enumerate(fh):
        if line.startswith("00Annotation"): hdr = line.rstrip("\n").split("\t"); skip = k; break
TPM = pd.read_csv(F5/"hg38_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz", sep="\t", skiprows=skip, usecols=[0] + need_cols, index_col=0, dtype={hdr[c]: np.float32 for c in need_cols})
TPM = TPM[~TPM.index.str.startswith("01STAT") & ~TPM.index.str.startswith("02STAT")]
P(f"FANTOM5 CAGE 피크 {len(TPM):,}개 × 라이브러리 {len(need_cols)}개, BED 이름 일치 {TPM.index.isin(bed.name).mean():.4f}")
bed = bed.set_index("name").loc[TPM.index].reset_index()[["chrom", "start", "end"]]
with gzip.open(F5/"F5.hg38.enhancers.expression.usage.matrix.gz", "rt") as fh: uh = fh.readline().rstrip("\n").split("\t")
uh = uh[1:] if uh[0] in ("", "id") else uh                   # 머리줄은 CNhs ID만(첫 열 이름 없음)
need_cnhs = sorted(set(c for v in lib_of.values() for c in v.cnhs) & set(uh))
USE = pd.read_csv(F5/"F5.hg38.enhancers.expression.usage.matrix.gz", sep="\t", header=None, skiprows=1, index_col=0, names=["id"] + uh,
                  usecols=["id"] + need_cnhs, dtype={c: np.int8 for c in need_cnhs})
eb = pd.read_csv(F5/"F5.hg38.enhancers.bed.gz", sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "name"]).set_index("name").loc[USE.index].reset_index()[["chrom", "start", "end"]]
P(f"FANTOM5 인핸서 {len(USE):,}개 × 라이브러리 {USE.shape[1]}개")
cage_rows = []
for g in LING:
    terms = TS[(TS.lineage == g) & (TS["mod"] == "CAGE")][["term", "biosample_type"]].itertuples(index=False); pk = np.zeros(len(TPM), bool); en = np.zeros(len(USE), bool)
    for t, b in terms:
        L = lib_of[(t, b)]
        if len(L): pk |= ((TPM[[hdr[c] for c in L.col]] >= 1).sum(1) >= majority(len(L))).to_numpy()
        u = [c for c in L.cnhs if c in USE.columns]
        if u: en |= (USE[u].sum(1) >= majority(len(u))).to_numpy()
        cage_rows.append(dict(lineage=g, term=t, n_lib=len(L), n_lib_enh=len(u)))
    el = pd.concat([bed[pk], eb[en]], ignore_index=True)
    ANN[f"M_{g}_CAGE"] = cover(el, WIN["CAGE"]); ANN[f"P_{g}_CAGE"] = cover(el, 0)
    P(f"  {g:11s} CAGE: 활성 피크 {pk.sum():,}, 활성 인핸서 {en.sum():,} → M {ANN[f'M_{g}_CAGE'].mean():.4f}, P {ANN[f'P_{g}_CAGE'].mean():.5f} (참조 SNP 비율)")
del TPM, USE
# ---------- ENCODE DNASE · H3K27ac ----------
def read_peaks(f):
    return pd.read_csv(EN/f"{f}.bed.gz", sep="\t", header=None, usecols=[0, 1, 2])
for m in ["DNASE", "H3K27ac"]:
    for g in LING:
        Mg = np.zeros(N, bool); Pg = np.zeros(N, bool); E = EF[(EF.lineage == g) & (EF["mod"] == m)]
        for t, Et in E.groupby("term"):
            exps = Et.groupby("experiment").file.apply(list); k = majority(len(exps)); cm = np.zeros(N, np.int16); cp = np.zeros(N, np.int16)
            for fs in exps:
                d = pd.concat([read_peaks(f) for f in fs], ignore_index=True); cm += cover(d, WIN[m]); cp += cover(d, 0)
            Mg |= cm >= k; Pg |= cp >= k
        ANN[f"M_{g}_{m}"] = Mg; ANN[f"P_{g}_{m}"] = Pg
        P(f"  {g:11s} {m:8s}: 용어 {E.term.nunique()}, 실험 {E.experiment.nunique()} → M {Mg.mean():.4f}, P {Pg.mean():.4f}")
NAMES = [f"{k}_{g}_{m}" for k in ["M", "P"] for g in LING for m in MODS]
A = pd.DataFrame({n: ANN[n].astype(np.int8) for n in NAMES})
# ---------- sd (MAF ≥ 5 %), 겹침 (기술) ----------
frq = pd.concat([pd.read_csv(PF/f"1000G.EUR.hg38.{c}.frq", sep=r"\s+", usecols=["SNP", "MAF"]) for c in CHRS], ignore_index=True)
snp = pd.concat([BIM[c].SNP for c in CHRS], ignore_index=True); assert (frq.SNP.to_numpy() == snp.to_numpy()).all(), "frq 와 bim 순서가 다르다"
com = (frq.MAF >= .05).to_numpy()
A[com].std(ddof=0).rename("sd_common").to_csv(W/("test_measured_sd_common.csv" if TEST else "measured_sd_common.csv"))
P("\n주석 크기 (MAF ≥ 5 % 참조 SNP 중 비율)\n" + A[com].mean().to_frame("frac").T.to_string(float_format=lambda v: f"{v:.4f}"))
L5 = {g: np.concatenate([pd.read_csv(Y/"annot"/f"L5_{g}_CAGE"/f"L5_{g}_CAGE.{c}.annot.gz", sep="\t").iloc[:, 0].to_numpy() for c in CHRS]).astype(bool) for g in LING}
P("\nAtlas L5_{g}_CAGE (MAF ≥ 5 %) 중 실측 주석 안 비율")
for g in LING:
    a = L5[g] & com; anyM = (A[[f"M_{g}_{m}" for m in MODS]].to_numpy().any(1))
    P(f"  {g:11s} L5 {a.sum():,}개: M_CAGE 안 {(a & ANN[f'M_{g}_CAGE']).sum() / a.sum():.3f}, 세 실측 중 하나 {(a & anyM).sum() / a.sum():.3f}, "
      f"P_CAGE 안 {(a & ANN[f'P_{g}_CAGE']).sum() / a.sum():.3f} | 비교: 전체 공통 SNP의 M_CAGE 비율 {ANN[f'M_{g}_CAGE'][com].mean():.3f}")
# ---------- LD 점수 ----------
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1"); i0 = 0
groups = {n: [n] for n in NAMES}
for g in LING: groups[f"BM_{g}"] = [f"M_{g}_{m}" for m in MODS]; groups[f"BMp_{g}"] = [f"P_{g}_{m}" for m in MODS]
for c in CHRS:
    n_c = len(BIM[c]); Ac = A.iloc[i0:i0 + n_c].reset_index(drop=True); i0 += n_c
    out_all = AN/"_all"/f"all.{c}"; (AN/"_all").mkdir(exist_ok=True)
    if not Path(f"{out_all}.l2.ldscore.gz").exists():
        Ac.to_csv(f"{out_all}.annot.gz", sep="\t", index=False, compression="gzip")
        r = subprocess.run(["python", "-W", "ignore", str(LDSC), "--l2", "--bfile", str(PF/f"1000G.EUR.hg38.{c}"), "--ld-wind-cm", "1", "--annot", f"{out_all}.annot.gz", "--thin-annot",
                            "--out", str(out_all), "--print-snps", str(X/"hm3_no_MHC.list.txt")], cwd=LDSC.parent, env=env, capture_output=True, text=True)
        if not Path(f"{out_all}.l2.ldscore.gz").exists(): raise SystemExit(f"chr{c} l2 실패\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}")
    L = pd.read_csv(f"{out_all}.l2.ldscore.gz", sep="\t"); M = open(f"{out_all}.l2.M").read().split(); M5 = open(f"{out_all}.l2.M_5_50").read().split()
    for gname, cols in groups.items():
        d = AN/gname; d.mkdir(exist_ok=True); idx = [NAMES.index(x) for x in cols]
        L[["CHR", "SNP", "BP"] + [f"{x}L2" for x in cols]].to_csv(d/f"{gname}.{c}.l2.ldscore.gz", sep="\t", index=False, compression="gzip")
        Ac[cols].to_csv(d/f"{gname}.{c}.annot.gz", sep="\t", index=False, compression="gzip")
        open(d/f"{gname}.{c}.l2.M", "w").write("\t".join(M[i] for i in idx)); open(d/f"{gname}.{c}.l2.M_5_50", "w").write("\t".join(M5[i] for i in idx))
    P(f"chr{c} 완료")
pd.DataFrame(cage_rows).to_csv(W/"fantom5_term_libraries.csv", index=False)
(EXP/"results"/("w02_annot_test.txt" if TEST else "w02_annot.txt")).write_text("\n".join(log), encoding="utf-8")
