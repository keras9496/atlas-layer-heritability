"""실험 3l l02: 예측 활성도 주석 A5_{g}_{m} → 염색체별 LD 점수 (PROTOCOL_exp3l_2026-10-17, 2026-10-18 고정).
 활성도 점수: l01의 ACTIVE 원점수(max(ΣALT, ΣREF)), 주 지도 전체 트랙 주석(3d L5_{g}_{m})과 같은 트랙.
 SNP 점수 = 트랙마다 대상 SNP에서 매긴 백분위(pandas rank pct, 3d와 같은 기본값)의 최댓값 → 다시 백분위 → 상위 5 % = A5_{g}_{m}.
 계통 5개(골수계, 뇌 조직, 림프계, 간, 결합조직) × 층위 3(CAGE, DNase, H3K27ac) = 15개. 묶음 A5_{g} = 3층위.
 확인(어긋나면 멈춘다): (1) 대상 SNP 전부에 활성도 값이 있다(결측 0). (2) 같은 코드로 3c 효과 점수에서 다시 만든 L5가 3d L5 주석과 모든 SNP에서 같다.
 단계: python l02_annot_ldscores.py check | main
   check = 확인 (2)만 (3c 자료로, 활성도 자료 불필요) / main = 활성도 주석 + 확인 (1)·(2) + LD 점수(3d·3k와 같은 l2 호출)
출력: exp3l/annot/A5_{g}/A5_{g}.<chr>.{annot.gz,l2.ldscore.gz,l2.M,l2.M_5_50}, exp3l/annot_sd_common.csv, results/l02_<단계>.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import os, sys, subprocess
import numpy as np, pandas as pd
import pyarrow.parquet as pq
from pathlib import Path
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); C3 = Path(f"{DATA}/exp3c"); Y = Path(f"{DATA}/exp3d"); L = Path(f"{DATA}/exp3l")
PF = X/"GRCh38"/"plink_files"; LDSC = Path(f"{DATA}/tools/ldsc/ldsc.py"); AN = L/"annot"; AN.mkdir(parents=True, exist_ok=True); (EXP/"results").mkdir(exist_ok=True)
LING = ["myeloid", "CNS_tissue", "lymphoid", "liver", "stromal"]; MODS = ["CAGE", "DNASE", "H3K27ac"]
ACT = {"CAGE": "CAGE_ACTIVE", "DNASE": "DNASE_ACTIVE", "CHIP_HISTONE": "CHIP_HISTONE_ACTIVE"}
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
TL = pd.read_csv(C3/"track_lineage.csv")
EFF = {(g, m): (TL.scorer + "::" + TL.key)[(TL.lineage == g) & (TL["mod"] == m)].tolist() for g in LING for m in MODS}            # 3d y03와 같은 열
ACC = {(g, m): (TL.scorer.map(ACT) + "::" + TL.key)[(TL.lineage == g) & (TL["mod"] == m)].tolist() for g in LING for m in MODS}
T = pd.read_parquet(X/"targets.parquet")[["chr", "rsid", "atlas_variant"]]
def read_cols(src, cols):
    """조각들에서 필요한 열만 읽어 atlas_variant 색인(중복 제거) 표로, 대상 SNP 순서로 붙인다 (3d y03와 같은 방식)"""
    V = []
    for f in sorted(src.glob("part_*.parquet")):
        have = set(pq.read_schema(f).names); miss = [c for c in cols if c not in have]; assert not miss, f"{f.name}: 열 없음 {miss[:3]}"
        V.append(pd.read_parquet(f, columns=["atlas_variant"] + cols).set_index("atlas_variant").astype(np.float32))
    V = pd.concat(V); V = V[~V.index.duplicated()]
    return T[["atlas_variant"]].join(V, on="atlas_variant").drop(columns="atlas_variant")
def top5(s): return (s.rank(pct=True) >= .95).astype(np.int8)
# ---------- 확인 (2): 효과 점수 → L5 재현 ----------
ann_eff = {}
for g in LING:
    cols = sorted(set(c for m in MODS for c in EFF[(g, m)])); D = read_cols(C3/"atlas", cols)
    for m in MODS: ann_eff[f"L5_{g}_{m}"] = top5(D[EFF[(g, m)]].max(1)).to_numpy()
    P(f"효과 점수 L5 재계산: {g} (트랙 {len(cols)})")
E = pd.DataFrame(ann_eff); E["rsid"] = T.rsid.to_numpy(); E = E.drop_duplicates("rsid").set_index("rsid")
bad = 0
for c in range(1, 23):
    b = pd.read_csv(PF/f"1000G.EUR.hg38.{c}.bim", sep="\t", header=None, usecols=[1], names=["SNP"]).SNP
    R = E.reindex(b).fillna(0).astype(np.int8).reset_index(drop=True)
    for n in E.columns:
        ref = pd.read_csv(Y/"annot"/n/f"{n}.{c}.annot.gz", sep="\t").iloc[:, 0].to_numpy(); bad += int((ref != R[n].to_numpy()).sum())
P(f"확인 (2) L5 재현: 다른 칸 {bad}개 (15개 주석 × 전체 bim SNP)"); assert bad == 0, "L5 재현 실패 — 멈춘다"
if STAGE == "check":
    (EXP/"results"/"l02_check.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit()
# ---------- 활성도 주석 ----------
nparts = len(list((L/"atlas_active").glob("part_*.parquet"))); assert nparts == 583, f"l01 조각 {nparts}/583"
ann, raw_desc = {}, []
for g in LING:
    cols = sorted(set(c for m in MODS for c in ACC[(g, m)])); D = read_cols(L/"atlas_active", cols)
    na = int(D.isna().sum().sum()); P(f"활성도 {g}: 트랙 {len(cols)}, 결측 칸 {na}"); assert na == 0, "확인 (1) 결측 — 멈춘다"
    for m in MODS:
        sc = None
        for c in ACC[(g, m)]:                                # 트랙별 백분위의 최댓값 (메모리 때문에 트랙 하나씩)
            r = D[c].rank(pct=True).to_numpy(); sc = r if sc is None else np.maximum(sc, r)
        ann[f"A5_{g}_{m}"] = top5(pd.Series(sc)).to_numpy()
        raw_desc.append(dict(lineage=g, layer=m, n_tracks=len(ACC[(g, m)]), raw_median=float(D[ACC[(g, m)]].median().median()), raw_max=float(D[ACC[(g, m)]].max().max())))
    del D
A = pd.DataFrame(ann); A["rsid"] = T.rsid.to_numpy(); A = A.drop_duplicates("rsid").set_index("rsid")
NAMES = list(ann); sd = A[NAMES].std(ddof=0).rename("sd_common"); sd.to_csv(L/"annot_sd_common.csv")
P("활성도 주석 크기(대상 SNP 중 %): " + ", ".join(f"{n} {A[n].mean() * 100:.2f}" for n in NAMES))
P("효과 L5와의 겹침 (같은 계통·층위, 대상 SNP): " + ", ".join(f"{g}_{m} Jaccard {((A[f'A5_{g}_{m}'] & E[f'L5_{g}_{m}']).sum() / (A[f'A5_{g}_{m}'] | E[f'L5_{g}_{m}']).sum()):.2f}" for g in LING for m in MODS))
pd.DataFrame(raw_desc).to_csv(EXP/"results"/"l02_raw_scale.csv", index=False)
# ---------- LD 점수 (3d·3k와 같은 호출) ----------
groups = {f"A5_{g}": [f"A5_{g}_{m}" for m in MODS] for g in LING}
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
for c in range(1, 23):
    b = pd.read_csv(PF/f"1000G.EUR.hg38.{c}.bim", sep="\t", header=None, usecols=[1], names=["SNP"]).SNP
    Ac = A.reindex(b).fillna(0).astype(np.int8).reset_index(drop=True)[NAMES]; out_all = AN/"_all"/f"all.{c}"; (AN/"_all").mkdir(exist_ok=True)
    if not Path(f"{out_all}.l2.ldscore.gz").exists():
        Ac.to_csv(f"{out_all}.annot.gz", sep="\t", index=False, compression="gzip")
        r = subprocess.run(["python", "-W", "ignore", str(LDSC), "--l2", "--bfile", str(PF/f"1000G.EUR.hg38.{c}"), "--ld-wind-cm", "1", "--annot", f"{out_all}.annot.gz", "--thin-annot",
                            "--out", str(out_all), "--print-snps", str(X/"hm3_no_MHC.list.txt")], cwd=LDSC.parent, env=env, capture_output=True, text=True)
        if not Path(f"{out_all}.l2.ldscore.gz").exists(): raise SystemExit(f"chr{c} l2 실패\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}")
    Ld = pd.read_csv(f"{out_all}.l2.ldscore.gz", sep="\t"); M = open(f"{out_all}.l2.M").read().split(); M5 = open(f"{out_all}.l2.M_5_50").read().split()
    for gname, cols in groups.items():
        d = AN/gname; d.mkdir(exist_ok=True); idx = [NAMES.index(x) for x in cols]
        Ld[["CHR", "SNP", "BP"] + [f"{x}L2" for x in cols]].to_csv(d/f"{gname}.{c}.l2.ldscore.gz", sep="\t", index=False, compression="gzip")
        Ac[cols].to_csv(d/f"{gname}.{c}.annot.gz", sep="\t", index=False, compression="gzip")
        open(d/f"{gname}.{c}.l2.M", "w").write("\t".join(M[i] for i in idx)); open(d/f"{gname}.{c}.l2.M_5_50", "w").write("\t".join(M5[i] for i in idx))
    P(f"chr{c} LD 점수 완료")
(EXP/"results"/"l02_main.txt").write_text("\n".join(log), encoding="utf-8")
