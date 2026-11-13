"""실험 3w w03: 활성도 구간 주석(PROTOCOL_exp3w 4.3)과 두 가지 집계(4.1 효과·활성도 주석의 겹침, 4.2 효과 크기 구간 SNP의 성질).
 점수는 3q q01과 같은 방법으로 원자료에서 다시 계산한다: 효과 e = 칸 안 트랙 |효과| 최댓값(3d L5와 같은 트랙), R_e = e의 백분위,
   활성도 R_a = 트랙별 대상 SNP 백분위의 최댓값을 다시 백분위로(3l A5, 3p Ac와 같다).
 새 주석(다섯 묶음 × 세 층위 × 네 구간 = 60개): AB1 = R_a ≥ 0.99, AB2 = 0.95 ≤ R_a < 0.99, AB3 = 0.90 ≤ R_a < 0.95, AB4 = 0.80 ≤ R_a < 0.90.
 확인(어긋나면 멈춘다): (1) AB1 ∪ AB2 = 3l A5, e 상위 5% = 3d L5 (전체 bim SNP), (2) 구간 크기 1, 4, 5, 10% ± 0.05%p.
 집계: (4.1) L5(효과 상위 5%) SNP 가운데 R_a ≥ 0.95, 0.90, 0.80인 비율. (4.2) 효과 크기 구간(R_e ≥ 0.999, 0.995–0.999, 0.99–0.995, 0.95–0.99)과 전체 SNP에 대해
   TSS 1 kb 안 비율·거리 중앙값(3k의 거리 표), 같은 층위 R_a의 중앙값과 R_a ≥ 0.90 비율, baselineLD 주석 일곱과의 겹침, MAF 평균.
 단계: python -I -X utf8 w03_actbin.py check | main <k> <K>   (check = 점수·주석·확인·집계, 주석 표 저장 / main = 염색체 조각의 LD 점수)
 출력: exp3w/annot/AB_<g>/..., exp3w/annot_table.parquet, exp3w/annot_sd_common.csv, results/w03_check.txt, w03_overlap.csv, w03_bins.csv
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import os, sys, subprocess
import numpy as np, pandas as pd
import pyarrow.parquet as pq
from pathlib import Path
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); C3 = Path(f"{DATA}/exp3c"); Y = Path(f"{DATA}/exp3d")
K3 = Path(f"{DATA}/exp3k"); L = Path(f"{DATA}/exp3l"); O = Path(f"{DATA}/exp3w"); AN = O/"annot"; AN.mkdir(parents=True, exist_ok=True); (EXP/"results").mkdir(exist_ok=True)
PF = X/"GRCh38"/"plink_files"; LDSC = Path(f"{DATA}/tools/ldsc/ldsc.py")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LING = ["myeloid", "CNS_tissue", "lymphoid", "liver", "stromal"]; MODS = ["CAGE", "DNASE", "H3K27ac"]
ACT = {"CAGE": "CAGE_ACTIVE", "DNASE": "DNASE_ACTIVE", "CHIP_HISTONE": "CHIP_HISTONE_ACTIVE"}
ABIN = {"AB1": (.99, 1.01), "AB2": (.95, .99), "AB3": (.90, .95), "AB4": (.80, .90)}; WANT = {"AB1": .01, "AB2": .04, "AB3": .05, "AB4": .10}
EBIN = {"0-0.1%": (.999, 1.01), "0.1-0.5%": (.995, .999), "0.5-1%": (.99, .995), "1-5%": (.95, .99)}
BASECOLS = ["Promoter_UCSC", "TSS_Hoffman", "Coding_UCSC", "Enhancer_Hoffman", "H3K27ac_Hnisz", "DHS_Trynka", "Conserved_LindbladToh"]
TL = pd.read_csv(C3/"track_lineage.csv")
EFF = {(g, m): (TL.scorer + "::" + TL.key)[(TL.lineage == g) & (TL["mod"] == m)].tolist() for g in LING for m in MODS}
ACC = {(g, m): (TL.scorer.map(ACT) + "::" + TL.key)[(TL.lineage == g) & (TL["mod"] == m)].tolist() for g in LING for m in MODS}
T = pd.read_parquet(X/"targets.parquet")[["chr", "rsid", "atlas_variant", "maf"]]
NAMES = [f"{b}_{g}_{m}" for g in LING for b in ABIN for m in MODS]                       # 묶음 안에서 구간이 먼저(3o B12와 같은 순서)
GROUPS = {f"AB_{g}": [f"{b}_{g}_{m}" for b in ABIN for m in MODS] for g in LING}
TAB = O/"annot_table.parquet"
def read_cols(src, cols):
    V = []
    for f in sorted(src.glob("part_*.parquet")):
        have = set(pq.read_schema(f).names); miss = [c for c in cols if c not in have]; assert not miss, f"{f.name}: 열 없음 {miss[:3]}"
        V.append(pd.read_parquet(f, columns=["atlas_variant"] + cols).set_index("atlas_variant").astype(np.float32))
    V = pd.concat(V); V = V[~V.index.duplicated()]
    return T[["atlas_variant"]].join(V, on="atlas_variant").drop(columns="atlas_variant")
def bim_snps(c): return pd.read_csv(PF/f"1000G.EUR.hg38.{c}.bim", sep="\t", header=None, usecols=[1], names=["SNP"]).SNP

if STAGE == "check":
    # ---------- 효과 e, R_e (3q q01과 같다) ----------
    need = sorted(set(c for v in EFF.values() for c in v)); parts = sorted((C3/"atlas").glob("part_*.parquet")); P(f"효과: 트랙 열 {len(need)}개, 조각 {len(parts)}개")
    V = []
    for f in parts:
        have = set(pq.read_schema(f).names); cols = [c for c in need if c in have]
        D = pd.read_parquet(f, columns=["atlas_variant"] + cols).set_index("atlas_variant").astype(np.float32)
        V.append(pd.DataFrame({f"{g}_{m}": (D[[c for c in v if c in have]].max(1) if any(c in have for c in v) else np.nan) for (g, m), v in EFF.items()}, index=D.index))
    V = pd.concat(V); V = V[~V.index.duplicated()]
    Ev = T[["atlas_variant"]].join(V, on="atlas_variant").drop(columns="atlas_variant").reset_index(drop=True); del V
    assert int(Ev.isna().sum().sum()) == 0, "효과 결측"
    Re = Ev.rank(pct=True)
    # ---------- 활성도 R_a (3q q01과 같다) ----------
    Ra = {}
    for g in LING:
        cols = sorted(set(c for m in MODS for c in ACC[(g, m)])); D = read_cols(L/"atlas_active", cols); assert int(D.isna().sum().sum()) == 0, "활성도 결측"
        for m in MODS:
            sc = None
            for c in ACC[(g, m)]:
                r = D[c].rank(pct=True).to_numpy(); sc = r if sc is None else np.maximum(sc, r)
            Ra[f"{g}_{m}"] = pd.Series(sc).rank(pct=True).to_numpy()
        del D; P(f"활성도 {g} 완료")
    Ra = pd.DataFrame(Ra)
    # ---------- 주석 ----------
    ann = pd.DataFrame({"rsid": T.rsid.to_numpy()}); chk = {}
    P("주석 크기(대상 SNP %)")
    for g in LING:
        for m in MODS:
            k = f"{g}_{m}"; ra = Ra[k].to_numpy()
            chk[f"A5_{k}"] = (ra >= .95).astype(np.int8); chk[f"L5_{k}"] = (Re[k].to_numpy() >= .95).astype(np.int8)
            for b, (lo, hi) in ABIN.items():
                a = ((ra >= lo) & (ra < hi)).astype(np.int8); ann[f"{b}_{k}"] = a; assert abs(a.mean() - WANT[b]) <= .0005, (b, k, a.mean())
            assert np.array_equal((ann[f"AB1_{k}"] | ann[f"AB2_{k}"]).to_numpy(), chk[f"A5_{k}"])
            P(f"  {k:20s} " + " ".join(f"{b} {ann[f'{b}_{k}'].mean() * 100:.3f}" for b in ABIN))
    # ---------- 확인 (1): 기존 주석 재현 (전체 bim SNP) ----------
    Cq = pd.DataFrame(chk); Cq["rsid"] = T.rsid.to_numpy(); Cq = Cq.drop_duplicates("rsid").set_index("rsid"); bad = {"A5": 0, "L5": 0}
    for c in range(1, 23):
        b = bim_snps(c); R = Cq.reindex(b).fillna(0).astype(np.int8).reset_index(drop=True)
        for g in LING:
            ref = pd.read_csv(L/"annot"/f"A5_{g}"/f"A5_{g}.{c}.annot.gz", sep="\t")
            for m in MODS: bad["A5"] += int((ref[f"A5_{g}_{m}"].to_numpy() != R[f"A5_{g}_{m}"].to_numpy()).sum())
            for m in MODS:
                n = f"L5_{g}_{m}"; ref = pd.read_csv(Y/"annot"/n/f"{n}.{c}.annot.gz", sep="\t").iloc[:, 0].to_numpy(); bad["L5"] += int((ref != R[n].to_numpy()).sum())
    P(f"확인 (1) 3l A5 재현: 다른 칸 {bad['A5']}개 / 3d L5 재현: 다른 칸 {bad['L5']}개"); assert bad["A5"] == 0 and bad["L5"] == 0, "재현 실패 — 멈춘다"
    P("확인 (1)·(2) 통과")
    # ---------- 집계 4.1: 효과 주석 SNP가 활성도 상위에 드는 비율 ----------
    rows = []
    for g in LING:
        for m in MODS:
            k = f"{g}_{m}"; e5 = chk[f"L5_{k}"] == 1; ra = Ra[k].to_numpy(); a5 = ra >= .95
            rows.append(dict(group=g, layer=m, n_effect=int(e5.sum()), in_top5=float((ra[e5] >= .95).mean()), in_top10=float((ra[e5] >= .90).mean()), in_top20=float((ra[e5] >= .80).mean()),
                             jaccard_top5=float((e5 & a5).sum() / (e5 | a5).sum()), median_Ra=float(np.median(ra[e5]))))
    OV = pd.DataFrame(rows); OV.to_csv(EXP/"results"/"w03_overlap.csv", index=False)
    P("\n집계 4.1: 효과 상위 5% SNP 가운데 같은 칸의 활성도 상위 5% / 10% / 20%에 드는 비율")
    for r in OV.itertuples(): P(f"  {r.group:11s} {r.layer:8s} {r.in_top5 * 100:5.1f}% / {r.in_top10 * 100:5.1f}% / {r.in_top20 * 100:5.1f}%  (Jaccard {r.jaccard_top5:.2f})")
    P("  범위: " + ", ".join(f"{m} 상위 10% 안 {OV[OV.layer == m].in_top10.min() * 100:.0f}–{OV[OV.layer == m].in_top10.max() * 100:.0f}%" for m in MODS))
    # ---------- 집계 4.2: 효과 크기 구간 SNP의 성질 ----------
    dist = pd.read_parquet(K3/"target_tss_dist.parquet").drop_duplicates("rsid").set_index("rsid").dist.reindex(T.rsid).to_numpy(); assert not np.isnan(dist).any()
    BA = pd.concat([pd.read_csv(X/"GRCh38"/"baselineLD_v2.2"/f"baselineLD.{c}.annot.gz", sep="\t", usecols=["SNP"] + BASECOLS) for c in range(1, 23)]).drop_duplicates("SNP").set_index("SNP")
    BAt = BA.reindex(T.rsid); inb = BAt[BASECOLS[0]].notna().to_numpy(); P(f"\nbaselineLD 주석이 있는 대상 SNP: {inb.mean() * 100:.2f}%"); BAv = BAt.to_numpy(float); maf = T.maf.to_numpy()
    rows = []
    def desc(sel, **kw):
        s2 = sel & inb
        return dict(kw, n=int(sel.sum()), tss_1kb=float((dist[sel] <= 1000).mean()), tss_dist_median=float(np.median(dist[sel])), maf_mean=float(maf[sel].mean()),
                    **{c: float(np.nanmean(BAv[s2, i])) for i, c in enumerate(BASECOLS)})
    rows.append(desc(np.ones(len(T), bool), group="all", layer="all", bin="all SNPs", Ra_median=.5, Ra_top10=.10))
    for g in LING:
        for m in MODS:
            k = f"{g}_{m}"; re_ = Re[k].to_numpy(); ra = Ra[k].to_numpy()
            for b, (lo, hi) in EBIN.items():
                sel = (re_ >= lo) & (re_ < hi); rows.append(desc(sel, group=g, layer=m, bin=b, Ra_median=float(np.median(ra[sel])), Ra_top10=float((ra[sel] >= .90).mean())))
    BN = pd.DataFrame(rows); BN.to_csv(EXP/"results"/"w03_bins.csv", index=False)
    P("집계 4.2: 효과 크기 구간 SNP — TSS 1 kb 안 비율 / 같은 층위 활성도 상위 10% 안 비율 / Promoter_UCSC / Conserved (전체 SNP: " +
      f"{BN.tss_1kb[0] * 100:.1f}% / 10.0% / {BN.Promoter_UCSC[0] * 100:.1f}% / {BN.Conserved_LindbladToh[0] * 100:.1f}%)")
    for g in LING:
        for m in ["CAGE", "DNASE"]:
            d = BN[(BN.group == g) & (BN.layer == m)].set_index("bin")
            P(f"  {g:11s} {m:6s} " + " | ".join(f"{b}: {d.tss_1kb[b] * 100:4.1f} / {d.Ra_top10[b] * 100:4.1f} / {d.Promoter_UCSC[b] * 100:4.1f} / {d.Conserved_LindbladToh[b] * 100:4.1f}" for b in ("0-0.1%", "1-5%")))
    ann = ann.drop_duplicates("rsid").set_index("rsid")[NAMES]; ann.to_parquet(TAB, compression="zstd"); ann.std(ddof=0).rename("sd_common").to_csv(O/"annot_sd_common.csv")
    P(f"\n주석 표 저장: {TAB} ({len(ann):,} SNP × {ann.shape[1]})")
    (EXP/"results"/"w03_check.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit()

# ---------------- main: LD 점수 (염색체 조각) ----------------
k, K = int(sys.argv[2]), int(sys.argv[3])
ann = pd.read_parquet(TAB); assert sorted(ann.columns) == sorted(NAMES), "주석 표 열이 다르다"; ann = ann[NAMES]
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
for c in list(range(1, 23))[k::K]:
    out_all = AN/"_all"/f"all.{c}"; (AN/"_all").mkdir(exist_ok=True)
    b = bim_snps(c); Ac = ann.reindex(b).fillna(0).astype(np.int8).reset_index(drop=True)
    if not Path(f"{out_all}.l2.ldscore.gz").exists():
        Ac.to_csv(f"{out_all}.annot.gz", sep="\t", index=False, compression="gzip")
        r = subprocess.run(["python", "-W", "ignore", str(LDSC), "--l2", "--bfile", str(PF/f"1000G.EUR.hg38.{c}"), "--ld-wind-cm", "1", "--annot", f"{out_all}.annot.gz", "--thin-annot",
                            "--out", str(out_all), "--print-snps", str(X/"hm3_no_MHC.list.txt")], cwd=LDSC.parent, env=env, capture_output=True, text=True)
        if not Path(f"{out_all}.l2.ldscore.gz").exists(): raise SystemExit(f"chr{c} l2 실패\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}")
    Ld = pd.read_csv(f"{out_all}.l2.ldscore.gz", sep="\t"); M = open(f"{out_all}.l2.M").read().split(); M5 = open(f"{out_all}.l2.M_5_50").read().split()
    for gname, cols in GROUPS.items():
        d = AN/gname; d.mkdir(exist_ok=True); idx = [NAMES.index(x) for x in cols]
        Ld[["CHR", "SNP", "BP"] + [f"{x}L2" for x in cols]].to_csv(d/f"{gname}.{c}.l2.ldscore.gz", sep="\t", index=False, compression="gzip")
        Ac[cols].to_csv(d/f"{gname}.{c}.annot.gz", sep="\t", index=False, compression="gzip")
        open(d/f"{gname}.{c}.l2.M", "w").write("\t".join(M[i] for i in idx)); open(d/f"{gname}.{c}.l2.M_5_50", "w").write("\t".join(M5[i] for i in idx))
    P(f"chr{c} 완료")
(EXP/"results"/f"w03_main_{k}of{K}.txt").write_text("\n".join(log), encoding="utf-8")
