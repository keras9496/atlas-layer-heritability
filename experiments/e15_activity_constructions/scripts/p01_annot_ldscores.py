"""실험 3p p01: 연속 활성도 Ac_{g}_{m}과 활성도 층화 효과 S5_{g}_{m} (PROTOCOL_exp3p_2026-10-23, 2026-10-23 고정) → 염색체별 LD 점수.
 효과 점수 e = 칸(계통 g × 층위 m) 안 트랙 |효과| 최댓값 (3o o01과 같은 읽기, 3d L5와 같은 트랙) → 대상 SNP 백분위 R_e.
 활성도 점수 a = l01 ACTIVE 원점수의 트랙별 대상 SNP 백분위의 최댓값 (3l l02와 같은 계산) → 대상 SNP 백분위 R_a.
 Ac = R_a (값 없음 = 0). S5 = R_a로 나눈 100구간(⌊100·R_a⌋, 최댓값은 99구간) 안에서 e의 백분위 ≥ 0.95.
 확인(어긋나면 멈춘다): (1) R_a ≥ 0.95 = 3l A5 (전체 bim SNP), (2) R_e ≥ 0.95 = 3d L5, (3) S5 크기 5.00 % ± 0.05 %p, 100구간 각각 5 % ± 0.5 %p, (4) Ac ∈ [0, 1].
 기술: S5–L5 Jaccard, 활성도 10분위별 L5·S5 비율.
 단계: python p01_annot_ldscores.py check | main
   check = 점수·주석·확인·기술, 주석 표를 exp3p/annot_table.parquet로 저장 / main = 저장한 주석 표로 LD 점수(3l·3o와 같은 l2 호출)
출력: exp3p/annot/{Ac3,S5}_{g}/..., exp3p/annot_sd_common.csv, results/p01_<단계>.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import os, sys, subprocess
import numpy as np, pandas as pd
import pyarrow.parquet as pq
from pathlib import Path
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); C3 = Path(f"{DATA}/exp3c"); Y = Path(f"{DATA}/exp3d")
L = Path(f"{DATA}/exp3l"); O = Path(f"{DATA}/exp3p"); AN = O/"annot"; AN.mkdir(parents=True, exist_ok=True); (EXP/"results").mkdir(exist_ok=True)
PF = X/"GRCh38"/"plink_files"; LDSC = Path(f"{DATA}/tools/ldsc/ldsc.py")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LING = ["myeloid", "CNS_tissue", "lymphoid", "liver", "stromal"]; MODS = ["CAGE", "DNASE", "H3K27ac"]
ACT = {"CAGE": "CAGE_ACTIVE", "DNASE": "DNASE_ACTIVE", "CHIP_HISTONE": "CHIP_HISTONE_ACTIVE"}
NSTR = 100
TL = pd.read_csv(C3/"track_lineage.csv")
EFF = {(g, m): (TL.scorer + "::" + TL.key)[(TL.lineage == g) & (TL["mod"] == m)].tolist() for g in LING for m in MODS}             # 3d·3o와 같은 열
ACC = {(g, m): (TL.scorer.map(ACT) + "::" + TL.key)[(TL.lineage == g) & (TL["mod"] == m)].tolist() for g in LING for m in MODS}    # 3l l02와 같은 열
T = pd.read_parquet(X/"targets.parquet")[["chr", "rsid", "atlas_variant"]]
NAMES = [f"Ac_{g}_{m}" for g in LING for m in MODS] + [f"S5_{g}_{m}" for g in LING for m in MODS]
TAB = O/"annot_table.parquet"

def read_cols(src, cols):
    """조각들에서 필요한 열만 읽어 atlas_variant 색인(중복 제거) 표로, 대상 SNP 순서로 붙인다 (3l l02와 같다)"""
    V = []
    for f in sorted(src.glob("part_*.parquet")):
        have = set(pq.read_schema(f).names); miss = [c for c in cols if c not in have]; assert not miss, f"{f.name}: 열 없음 {miss[:3]}"
        V.append(pd.read_parquet(f, columns=["atlas_variant"] + cols).set_index("atlas_variant").astype(np.float32))
    V = pd.concat(V); V = V[~V.index.duplicated()]
    return T[["atlas_variant"]].join(V, on="atlas_variant").drop(columns="atlas_variant")

def bim_snps(c): return pd.read_csv(PF/f"1000G.EUR.hg38.{c}.bim", sep="\t", header=None, usecols=[1], names=["SNP"]).SNP

if STAGE == "check":
    # ---------- 효과 점수 e와 R_e (3o o01과 같다) ----------
    need = sorted(set(c for v in EFF.values() for c in v)); parts = sorted((C3/"atlas").glob("part_*.parquet")); P(f"효과: 트랙 열 {len(need)}개, 조각 {len(parts)}개")
    V = []
    for f in parts:
        have = set(pq.read_schema(f).names); cols = [c for c in need if c in have]
        D = pd.read_parquet(f, columns=["atlas_variant"] + cols).set_index("atlas_variant").astype(np.float32)
        V.append(pd.DataFrame({f"{g}_{m}": (D[[c for c in v if c in have]].max(1) if any(c in have for c in v) else np.nan) for (g, m), v in EFF.items()}, index=D.index))
    V = pd.concat(V); V = V[~V.index.duplicated()]
    Ev = T[["atlas_variant"]].join(V, on="atlas_variant").drop(columns="atlas_variant").reset_index(drop=True); del V
    P(f"대상 {len(T):,}, 효과 값 있음 {int(Ev[f'{LING[0]}_{MODS[0]}'].notna().sum()):,}")
    Re = Ev.rank(pct=True)
    # ---------- 활성도 점수 a와 R_a (3l l02와 같다) ----------
    nparts = len(list((L/"atlas_active").glob("part_*.parquet"))); assert nparts == 583, f"l01 조각 {nparts}/583"
    Ra = {}
    for g in LING:
        cols = sorted(set(c for m in MODS for c in ACC[(g, m)])); D = read_cols(L/"atlas_active", cols)
        na = int(D.isna().sum().sum()); P(f"활성도 {g}: 트랙 {len(cols)}, 결측 칸 {na}"); assert na == 0, "활성도 결측 — 멈춘다"
        for m in MODS:
            sc = None
            for c in ACC[(g, m)]:
                r = D[c].rank(pct=True).to_numpy(); sc = r if sc is None else np.maximum(sc, r)
            Ra[f"{g}_{m}"] = pd.Series(sc).rank(pct=True).to_numpy()
        del D
    Ra = pd.DataFrame(Ra)
    # ---------- 주석 ----------
    ann = pd.DataFrame({"rsid": T.rsid.to_numpy()}); chk = {}
    strata_dev = []
    for g in LING:
        for m in MODS:
            k = f"{g}_{m}"; ra = Ra[k]; e = Ev[k]
            ann[f"Ac_{k}"] = ra.fillna(0).astype(np.float32).to_numpy()
            ok = ra.notna() & e.notna(); st = np.minimum(np.floor(ra[ok] * NSTR), NSTR - 1).astype(int)
            within = e[ok].groupby(st.to_numpy()).rank(pct=True)
            s5 = np.zeros(len(T), np.int8); s5[np.flatnonzero(ok.to_numpy())] = (within.to_numpy() >= .95).astype(np.int8); ann[f"S5_{k}"] = s5
            frac = pd.Series(s5[ok.to_numpy()]).groupby(st.to_numpy()).mean(); strata_dev.append((k, float(s5.mean()), float(frac.min()), float(frac.max()), int(st.nunique())))
            chk[f"A5_{k}"] = (ra >= .95).astype(np.int8).to_numpy(); chk[f"L5_{k}"] = (Re[k] >= .95).astype(np.int8).to_numpy()
    # ---------- 확인 (1)·(2): 3l A5, 3d L5 재현 (전체 bim SNP) ----------
    C = pd.DataFrame(chk); C["rsid"] = T.rsid.to_numpy(); C = C.drop_duplicates("rsid").set_index("rsid")
    bad = {"A5": 0, "L5": 0}
    for c in range(1, 23):
        b = bim_snps(c); R = C.reindex(b).fillna(0).astype(np.int8).reset_index(drop=True)
        for g in LING:
            ref = pd.read_csv(L/"annot"/f"A5_{g}"/f"A5_{g}.{c}.annot.gz", sep="\t")
            for m in MODS: bad["A5"] += int((ref[f"A5_{g}_{m}"].to_numpy() != R[f"A5_{g}_{m}"].to_numpy()).sum())
            for m in MODS:
                n = f"L5_{g}_{m}"; ref = pd.read_csv(Y/"annot"/n/f"{n}.{c}.annot.gz", sep="\t").iloc[:, 0].to_numpy(); bad["L5"] += int((ref != R[n].to_numpy()).sum())
    P(f"확인 (1) 3l A5 재현: 다른 칸 {bad['A5']}개 / 확인 (2) 3d L5 재현: 다른 칸 {bad['L5']}개 (15개 주석 × 전체 bim SNP)")
    assert bad["A5"] == 0 and bad["L5"] == 0, "재현 실패 — 멈춘다"
    # ---------- 확인 (3)·(4) ----------
    P("확인 (3) S5 크기와 100구간별 비율 (대상 SNP %; 구간 최소–최대, 구간 수)")
    worst = 0.0
    for k, tot, lo, hi, ns in strata_dev:
        P(f"  S5_{k:20s} {tot * 100:.3f} % | 구간 {lo * 100:.2f}–{hi * 100:.2f} % | {ns}구간")
        assert abs(tot - .05) <= .0005 and lo >= .045 and hi <= .055 and ns == NSTR, f"S5 크기 확인 실패: {k}"
    acv = ann[[n for n in NAMES if n.startswith("Ac_")]].to_numpy(); assert acv.min() >= 0 and acv.max() <= 1, "Ac 범위"
    P(f"확인 (4) Ac 범위 {acv.min():.3g}–{acv.max():.3g}"); P("확인 (1)–(4) 통과")
    # ---------- 기술 ----------
    def jac(a, b): return (a & b).sum() / max((a | b).sum(), 1)
    P("\n기술: S5와 L5의 겹침 (같은 계통·층위, 대상 SNP Jaccard)")
    for g in LING: P(f"  {g:11s} " + ", ".join(f"{m} {jac(ann[f'S5_{g}_{m}'].to_numpy(), chk[f'L5_{g}_{m}']):.2f}" for m in MODS))
    P("\n기술: 활성도 10분위별 주석 비율 (%; L5 → S5), 1분위(낮음) … 10분위(높음)")
    for g in LING:
        for m in MODS:
            k = f"{g}_{m}"; dec = np.minimum(np.floor(Ra[k].fillna(0).to_numpy() * 10), 9).astype(int)
            l5 = pd.Series(chk[f"L5_{k}"]).groupby(dec).mean() * 100; s5 = pd.Series(ann[f"S5_{k}"].to_numpy()).groupby(dec).mean() * 100
            P(f"  {k:20s} L5 " + " ".join(f"{v:4.1f}" for v in l5) + " | S5 " + " ".join(f"{v:4.1f}" for v in s5))
    ann = ann.drop_duplicates("rsid").set_index("rsid"); ann.to_parquet(TAB, compression="zstd")
    ann.std(ddof=0).rename("sd_common").to_csv(O/"annot_sd_common.csv")
    P(f"\n주석 표 저장: {TAB} ({len(ann):,} SNP × {ann.shape[1]})")
    (EXP/"results"/"p01_check.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit()

# ---------------- main: LD 점수 ----------------
ann = pd.read_parquet(TAB); assert sorted(ann.columns) == sorted(NAMES), "주석 표 열이 다르다"; ann = ann[NAMES]   # check 단계는 칸마다 Ac·S5 순서로 저장한다
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
groups = {**{f"Ac3_{g}": [f"Ac_{g}_{m}" for m in MODS] for g in LING}, **{f"S5_{g}": [f"S5_{g}_{m}" for m in MODS] for g in LING}}
CONT = [n for n in NAMES if n.startswith("Ac_")]
for c in range(1, 23):
    out_all = AN/"_all"/f"all.{c}"; (AN/"_all").mkdir(exist_ok=True)
    b = bim_snps(c); Ac = ann.reindex(b).fillna(0).reset_index(drop=True)
    for n in NAMES: Ac[n] = Ac[n].astype(np.float32) if n in CONT else Ac[n].astype(np.int8)
    if not Path(f"{out_all}.l2.ldscore.gz").exists():
        Ac.to_csv(f"{out_all}.annot.gz", sep="\t", index=False, compression="gzip", float_format="%.6g")
        r = subprocess.run(["python", "-W", "ignore", str(LDSC), "--l2", "--bfile", str(PF/f"1000G.EUR.hg38.{c}"), "--ld-wind-cm", "1", "--annot", f"{out_all}.annot.gz", "--thin-annot",
                            "--out", str(out_all), "--print-snps", str(X/"hm3_no_MHC.list.txt")], cwd=LDSC.parent, env=env, capture_output=True, text=True)
        if not Path(f"{out_all}.l2.ldscore.gz").exists(): raise SystemExit(f"chr{c} l2 실패\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}")
    Ld = pd.read_csv(f"{out_all}.l2.ldscore.gz", sep="\t"); M = open(f"{out_all}.l2.M").read().split(); M5 = open(f"{out_all}.l2.M_5_50").read().split()
    for gname, cols in groups.items():
        d = AN/gname; d.mkdir(exist_ok=True); idx = [NAMES.index(x) for x in cols]
        Ld[["CHR", "SNP", "BP"] + [f"{x}L2" for x in cols]].to_csv(d/f"{gname}.{c}.l2.ldscore.gz", sep="\t", index=False, compression="gzip")
        Ac[cols].to_csv(d/f"{gname}.{c}.annot.gz", sep="\t", index=False, compression="gzip", float_format="%.6g")
        open(d/f"{gname}.{c}.l2.M", "w").write("\t".join(M[i] for i in idx)); open(d/f"{gname}.{c}.l2.M_5_50", "w").write("\t".join(M5[i] for i in idx))
    P(f"chr{c} 완료"); (O/"PIPELINE_STATUS.txt").write_text(f"p01 chr{c} 완료", encoding="utf-8")
(EXP/"results"/"p01_main.txt").write_text("\n".join(log), encoding="utf-8")
