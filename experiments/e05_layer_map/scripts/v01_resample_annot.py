"""실험 3f v01: 트랙 균형 재추출 주석 (PROTOCOL_exp3f_2026-10-09) → 염색체별 LD 점수 (실험 3d y03과 같은 l2 호출).
 단위: DNASE·H3K27ac = 트랙, CAGE = biosample(가닥 +/− 트랙 2개, 키에서 '|가닥'을 뗀 것). 칸(계통 g × 층위 m)마다 단위 k = 3개 비복원 추출(단위 ≤ 3이면 전부).
 시드 20개: numpy default_rng(20261009) 하나로 시드 → 계통 → 층위 순서. 주석 L5r{s}_{g}_{m} = 뽑힌 트랙 |ALT − REF| 최댓값 → 대상 SNP 백분위 → 상위 5 %.
 묶음 T3_{s}_{g} = [CAGE, DNASE, H3K27ac] (100개). MAF < 5 % 또는 Atlas 값 없는 SNP = 0 (3d와 같다).
 시험: V01_TEST=1 → exp3c/test 조각, chr22 만, exp3f/test/annot 에 쓴다.
출력: exp3f/annot/T3_{s}_{g}/..., exp3f/annot_sd_common.csv, exp3f/resample_units.csv, results/v01_annot.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess, os
import numpy as np, pandas as pd
import pyarrow.parquet as pq
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); C3 = Path(f"{DATA}/exp3c"); F = Path(f"{DATA}/exp3f")
PF = X/"GRCh38"/"plink_files"; TEST = os.environ.get("V01_TEST") == "1"; LDSC = Path(f"{DATA}/tools/ldsc/ldsc.py")
SRC = C3/"test" if TEST else C3/"atlas"; AN = F/"test"/"annot" if TEST else F/"annot"; AN.mkdir(parents=True, exist_ok=True); (EXP/"results").mkdir(exist_ok=True)
CHRS = [22] if TEST else range(1, 23)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LING = ["myeloid", "CNS_tissue", "lymphoid", "liver", "stromal"]; MODS = ["CAGE", "DNASE", "H3K27ac"]; K, NSEED = 3, 20
TL = pd.read_csv(C3/"track_lineage.csv"); TL["col"] = TL.scorer + "::" + TL.key
TL["unit"] = [k.split("|")[0] if m == "CAGE" else k for k, m in zip(TL.key, TL["mod"])]
UNITS = {(g, m): TL[(TL.lineage == g) & (TL["mod"] == m)].groupby("unit").col.apply(list).sort_index() for g in LING for m in MODS}
P("단위 수 (CAGE = biosample, 나머지 = 트랙)\n" + pd.DataFrame({m: {g: len(UNITS[(g, m)]) for g in LING} for m in MODS}).to_string())
rng = np.random.default_rng(20261009); GROUPS, rec = {}, []
for s in range(1, NSEED + 1):
    for g in LING:
        for m in MODS:
            u = UNITS[(g, m)]; names = list(u.index)
            pick = sorted(rng.choice(names, K, replace=False).tolist()) if len(names) > K else names
            GROUPS[f"L5r{s:02d}_{g}_{m}"] = [c for n in pick for c in u[n]]
            rec.append(dict(seed=s, lineage=g, mod=m, n_units=len(names), n_picked=len(pick), units=";".join(pick)))
for g in LING:
    for m in MODS: GROUPS[f"FULL_{g}_{m}"] = [c for v in UNITS[(g, m)] for c in v]   # Jaccard 기술용 (LD 점수 안 만듦)
pd.DataFrame(rec).to_csv(F/("test_resample_units.csv" if TEST else "resample_units.csv"), index=False)
NAMES = [n for n in GROUPS if n.startswith("L5r")]; need = sorted(set(c for v in GROUPS.values() for c in v))
parts = sorted(SRC.glob("part_*.parquet")); P(f"주석 {len(NAMES)}개 (+ 기술용 전체 트랙 15개), 트랙 열 {len(need)}개, 조각 {len(parts)}개 읽기")
V = []
for f in parts:
    have = set(pq.read_schema(f).names); cols = [c for c in need if c in have]
    D = pd.read_parquet(f, columns=["atlas_variant"] + cols).set_index("atlas_variant").astype(np.float32)
    V.append(pd.DataFrame({n: (D[[c for c in v if c in have]].max(1) if any(c in have for c in v) else np.nan) for n, v in GROUPS.items()}, index=D.index))
V = pd.concat(V); V = V[~V.index.duplicated()]
T = pd.read_parquet(X/"targets.parquet")[["rsid", "atlas_variant"]]; T = T.join(V, on="atlas_variant")
P(f"대상 {len(T):,}, Atlas 값 있음 {int(T[NAMES[0]].notna().sum()):,}")
ann = pd.DataFrame({"rsid": T.rsid})
for n in GROUPS: ann[n] = (T[n].rank(pct=True) >= .95).astype(np.int8)
P("\n재추출 주석과 전체 트랙 주석의 겹침 (Jaccard, 시드 20개 중앙값 [최소–최대])")
for g in LING:
    for m in MODS:
        J = [((ann[f"L5r{s:02d}_{g}_{m}"] & ann[f"FULL_{g}_{m}"]).sum() / max((ann[f"L5r{s:02d}_{g}_{m}"] | ann[f"FULL_{g}_{m}"]).sum(), 1)) for s in range(1, NSEED + 1)]
        P(f"  {g:11s} {m:8s} {np.median(J):.2f} [{min(J):.2f}–{max(J):.2f}]")
ann = ann.drop_duplicates("rsid").set_index("rsid")[NAMES]
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
groups = {f"T3_{s:02d}_{g}": [f"L5r{s:02d}_{g}_{m}" for m in MODS] for s in range(1, NSEED + 1) for g in LING}
for c in CHRS:
    out_all = AN/"_all"/f"all.{c}"; (AN/"_all").mkdir(exist_ok=True)
    b = pd.read_csv(PF/f"1000G.EUR.hg38.{c}.bim", sep="\t", header=None, names=["CHR", "SNP", "CM", "BP", "a1", "a2"])
    Ac = ann.reindex(b.SNP).fillna(0).astype(np.int8).reset_index(drop=True)
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
ann.std(ddof=0).rename("sd_common").to_csv(F/("test_annot_sd_common.csv" if TEST else "annot_sd_common.csv"))
(EXP/"results"/("v01_annot_test.txt" if TEST else "v01_annot.txt")).write_text("\n".join(log), encoding="utf-8")
