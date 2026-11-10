"""실험 3s s01: FunGen-xQTL(NG00184.v1) 단일 맥락 SuSiE 결과 → 세포형 MaxCPP 주석 5개(Q5) → 염색체별 LD 점수 (PROTOCOL_exp3s_2026-10-25, 2026-10-25 고정).
 MaxCPP(Hormozdiari 2018 정의): SNP마다 묶음 안 모든 맥락·유전자의 PIP 최댓값, 보고되지 않은 SNP = 0. 연속 주석 [0, 1].
 묶음: Q_MIC = MiGA 미세아교 4부위 + ROSMAP 단일핵 미세아교(세 판) / Q_NEU = 단일핵 흥분·억제 뉴런 / Q_GLIA = 단일핵 성상세포·희소돌기세포·OPC /
       Q_BULK = bulk 뇌 8개(Knight-ADRC, MSBB 4, ROSMAP DLPFC·PCC·앞꼬리핵) / Q_MONO = ROSMAP 단핵구
 SNV 파일(*_snp_scfmAll*)만 쓴다. 1000G EUR hg38 bim(지금 주석과 같은 참조)에 염색체·위치(chromEnd, 1-based)·두 대립유전자(방향 무관)로 맞춘다. bim의 모든 SNP에 값을 준다.
 확인(어긋나면 멈춘다): (1) 맥락별 FunGen SNV 행(고유 변이) 중 bim에 맞은 비율 ≥ 50 %, (2) 값 범위 [0, 1].
 기술: 묶음별 MaxCPP > 0, > 0.1, > 0.5 SNP 수(대상 공통 SNP), 기준선 GTEx_eQTL_MaxCPP와의 상관(대상 공통 SNP).
 단계: python s01_annot_ldscores.py check | main
   check = 주석을 만들고 확인·기술, 염색체별 bim 순서 표를 exp3s/annot_bim/<c>.parquet, sd를 exp3s/annot_sd_common.csv에 저장
   main  = 저장한 표로 LD 점수(3l·3p와 같은 l2 호출, --ld-wind-cm 1, thin annot, HapMap3 no-MHC 출력)
출력: exp3s/annot/Q5/Q5.<c>.{annot.gz,l2.ldscore.gz,l2.M,l2.M_5_50}, results/s01_<단계>.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import os, sys, gzip, subprocess
import numpy as np, pandas as pd
from pathlib import Path
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); S3 = Path(f"{DATA}/exp3s"); RAW = S3/"raw"
PF = X/"GRCh38"/"plink_files"; LDSC = Path(f"{DATA}/tools/ldsc/ldsc.py"); AB = S3/"annot_bim"; AB.mkdir(parents=True, exist_ok=True); AN = S3/"annot"/"Q5"
(EXP/"results").mkdir(exist_ok=True)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
NAMES = ["Q_MIC", "Q_NEU", "Q_GLIA", "Q_BULK", "Q_MONO"]
def group_of(src, ctype):
    if src == "MiGA" or ctype == "Microglia": return "Q_MIC"
    if ctype in ("Excitatory neuron", "Inhibitory neuron"): return "Q_NEU"
    if ctype in ("Astrocyte", "Oligodendrocyte", "Oligodendrocyte progenitor cell"): return "Q_GLIA"
    if ctype == "Monocyte": return "Q_MONO"
    return "Q_BULK"
def bim(c): return pd.read_csv(PF/f"1000G.EUR.hg38.{c}.bim", sep="\t", header=None, usecols=[1, 3, 4, 5], names=["SNP", "pos", "a1", "a2"])

if STAGE == "check":
    M = pd.read_csv(RAW/"ADSP_FunGen_xQTL.v1.metadata.txt", sep="\t")
    M = M[M["File name"].str.contains("_snp_scfmAll") & M["File name"].str.contains(r"\.(?:eQTL|snuc-eQTL)\.", regex=True)].copy()
    M["ctx"] = M["File name"].str.split(".").str[1] + "." + M["File name"].str.split(".").str[3]          # 코호트.맥락 코드
    M["group"] = [group_of(s, c) for s, c in zip(M["Data Source"], M["cell type"])]
    P("맥락 → 묶음:")
    for (g, ctx, ct), n in M.groupby(["group", "ctx", "cell type"]).size().items(): P(f"  {g:7s} {ctx:40s} {ct} ({n}개 파일)")
    assert set(M.group) == set(NAMES), set(M.group)
    files = {p.name: p for p in RAW.rglob("*_snp_scfmAll*.bed.gz")}; miss = [f for f in M["File name"] if f not in files]; assert not miss, f"파일 없음 {miss[:3]}"
    parts = []
    for k, (fn, ctx, g) in enumerate(zip(M["File name"], M.ctx, M.group)):
        D = pd.read_csv(files[fn], sep="\t", usecols=["#chrom", "chromEnd", "ref", "alt", "PIP"], dtype={"#chrom": str})
        D = D[(D.ref.str.len() == 1) & (D.alt.str.len() == 1)]
        D = D.assign(chr=D["#chrom"].str.replace("chr", "", regex=False), pos=D.chromEnd.astype(np.int64),
                     key=np.where(D.ref < D.alt, D.ref + D.alt, D.alt + D.ref), ctx=ctx, group=g)[["chr", "pos", "key", "PIP", "ctx", "group"]]
        parts.append(D)
        if (k + 1) % 100 == 0: P(f"  파일 {k + 1}/{len(M)} 읽음")
    V = pd.concat(parts, ignore_index=True); del parts
    G = V.groupby(["group", "chr", "pos", "key"], sort=False).PIP.max().reset_index()
    P(f"파일 {len(M)}개 읽음, 변이–유전자–맥락 행 {len(V):,}. 묶음별 고유 SNV: " + ", ".join(f"{g} {(G.group == g).sum():,}" for g in NAMES))
    # ---- bim에 맞추기 (염색체·위치·대립유전자 쌍, 방향 무관) ----
    T = pd.read_parquet(X/"targets.parquet", columns=["rsid"]).rsid
    tabs, bkeys = {}, []
    for c in range(1, 23):
        Bm = bim(c); Bm["key"] = np.where(Bm.a1 < Bm.a2, Bm.a1 + Bm.a2, Bm.a2 + Bm.a1); Bm["row"] = np.arange(len(Bm))
        bkeys.append(Bm.assign(chr=str(c))[["chr", "pos", "key"]].drop_duplicates())
        W = G[G.chr == str(c)].merge(Bm[["pos", "key", "row"]], on=["pos", "key"])
        A = W.pivot_table(index="row", columns="group", values="PIP", aggfunc="max").reindex(index=np.arange(len(Bm)), columns=NAMES).fillna(0).astype(np.float32)
        tab = A.reset_index(drop=True); tab.insert(0, "SNP", Bm.SNP.to_numpy()); tab.to_parquet(AB/f"{c}.parquet", compression="zstd"); tabs[c] = tab
        P(f"  chr{c}: bim {len(Bm):,} SNP, 값 > 0 " + ", ".join(f"{g} {(tab[g] > 0).sum():,}" for g in NAMES))
    BK = pd.concat(bkeys, ignore_index=True).assign(inbim=True)
    # ---- 확인 (1): 맥락별 맞춘 비율 ----
    P("\n확인 (1): 맥락별 FunGen SNV(고유 변이) 중 bim에 맞은 비율")
    U = V[["ctx", "chr", "pos", "key"]].drop_duplicates().merge(BK, on=["chr", "pos", "key"], how="left"); U["inbim"] = U.inbim.fillna(False).astype(bool)
    bad = []
    for ctx, d in U.groupby("ctx"):
        fr = d.inbim.mean(); P(f"  {ctx:40s} {len(d):>9,} 변이, 맞음 {fr * 100:.1f} %")
        if fr < .5: bad.append(ctx)
    assert not bad, f"맞은 비율 50 % 미만 맥락: {bad} — 멈춘다"
    ALL = pd.concat(tabs.values(), ignore_index=True); v = ALL[NAMES].to_numpy(); assert v.min() >= 0 and v.max() <= 1, "값 범위"
    P(f"확인 (2) 값 범위 {v.min():.3g}–{v.max():.3g}"); P("확인 (1)–(2) 통과")
    # ---- 기술 ----
    TG = ALL.drop_duplicates("SNP").set_index("SNP").reindex(T).fillna(0)
    P(f"\n기술: 대상 공통 SNP {len(TG):,}개 중 MaxCPP > 0 / > 0.1 / > 0.5 (평균)")
    for g in NAMES: x = TG[g]; P(f"  {g:7s} {(x > 0).sum():>9,} / {(x > .1).sum():>7,} / {(x > .5).sum():>6,} (평균 {x.mean():.5f})")
    BLc = []
    for c in range(1, 23):
        b = pd.read_csv(X/"GRCh38"/"baselineLD_v2.2"/f"baselineLD.{c}.annot.gz", sep="\t", usecols=["SNP", "GTEx_eQTL_MaxCPP", "BLUEPRINT_H3K27acQTL_MaxCPP"]); BLc.append(b)
    BLc = pd.concat(BLc).drop_duplicates("SNP").set_index("SNP"); J_ = TG.join(BLc, how="inner")
    P(f"\n기술: 기준선 MaxCPP와의 상관(대상 공통 SNP 중 기준선에 있는 {len(J_):,}개)")
    for g in NAMES: P(f"  {g:7s} GTEx eQTL {np.corrcoef(J_[g], J_.GTEx_eQTL_MaxCPP)[0, 1]:+.3f} | BLUEPRINT H3K27ac {np.corrcoef(J_[g], J_.BLUEPRINT_H3K27acQTL_MaxCPP)[0, 1]:+.3f}")
    P("  묶음 사이 상관: " + "; ".join(f"{a}–{b} {np.corrcoef(TG[a], TG[b])[0, 1]:+.2f}" for i, a in enumerate(NAMES) for b in NAMES[i + 1:]))
    TG.std(ddof=0).rename("sd_common").to_csv(S3/"annot_sd_common.csv")
    (EXP/"results"/"s01_check.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit()

# ---------------- main: LD 점수 ----------------
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1"); AN.mkdir(parents=True, exist_ok=True); (S3/"_all").mkdir(exist_ok=True)
for c in range(1, 23):
    tab = pd.read_parquet(AB/f"{c}.parquet"); assert (tab.SNP.to_numpy() == bim(c).SNP.to_numpy()).all(), f"chr{c} bim 순서"
    A = tab[NAMES].astype(np.float32); out = AN/f"Q5.{c}"
    if not Path(f"{out}.l2.ldscore.gz").exists():
        A.to_csv(f"{out}.annot.gz", sep="\t", index=False, compression="gzip", float_format="%.6g")
        r = subprocess.run(["python", "-W", "ignore", str(LDSC), "--l2", "--bfile", str(PF/f"1000G.EUR.hg38.{c}"), "--ld-wind-cm", "1", "--annot", f"{out}.annot.gz", "--thin-annot",
                            "--out", str(out), "--print-snps", str(X/"hm3_no_MHC.list.txt")], cwd=LDSC.parent, env=env, capture_output=True, text=True)
        if not Path(f"{out}.l2.ldscore.gz").exists(): raise SystemExit(f"chr{c} l2 실패\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}")
    P(f"chr{c} 완료"); (S3/"PIPELINE_STATUS.txt").write_text(f"s01 chr{c} 완료", encoding="utf-8")
(EXP/"results"/"s01_main.txt").write_text("\n".join(log), encoding="utf-8")
