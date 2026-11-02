"""실험 3i i01: 세포형 일치 주석 (PROTOCOL_exp3i_2026-10-13) → 염색체별 LD 점수 (3f v01·3g u01과 같은 l2 호출).
 묶음 × 층위마다 단위 하나(CAGE = 라이브러리 하나의 두 가닥 트랙, DNase·H3K27ac = 트랙 하나). 후보가 여럿이면 표에 먼저 적은 용어 → hCAGE를 LQhCAGE보다 먼저 → 키 정렬 순서.
 칸(계통 g × 층위 m)마다 선택 단위들의 |log2 ALT − log2 REF| 최댓값 → 대상 SNP 백분위 → 상위 5 % = M5_{g}_{m}. 묶음 M3_{g} = [CAGE, DNASE, H3K27ac].
 확인: 같은 방식의 전체 트랙 상위 5 % sd가 3d annot_sd_common.csv의 L5_{g}_{m}과 같아야 한다(차이 < 1e-9, 아니면 멈춘다).
 기술: 선택 단위 표, 주석 크기, 일치 주석과 전체 트랙 주석(L5)의 Jaccard, 같은 계통 안 층위 간 Jaccard.
 시험: I01_TEST=1 → exp3c/test 조각, chr22 만, exp3i/test/annot 에 쓴다 (5 % 확인 생략).
출력: exp3i/annot/M3_{g}/..., exp3i/annot_sd_common.csv, exp3i/selected_units.csv, results/i01_annot.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess, os, re
import numpy as np, pandas as pd
import pyarrow.parquet as pq
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); C3 = Path(f"{DATA}/exp3c"); Y = Path(f"{DATA}/exp3d"); I = Path(f"{DATA}/exp3i")
PF = X/"GRCh38"/"plink_files"; TEST = os.environ.get("I01_TEST") == "1"; LDSC = Path(f"{DATA}/tools/ldsc/ldsc.py")
SRC = C3/"test" if TEST else C3/"atlas"; AN = I/"test"/"annot" if TEST else I/"annot"; AN.mkdir(parents=True, exist_ok=True); (EXP/"results").mkdir(exist_ok=True)
CHRS = [22] if TEST else range(1, 23)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LING = ["myeloid", "CNS_tissue", "lymphoid", "liver", "stromal"]; MODS = ["CAGE", "DNASE", "H3K27ac"]
def same(t): return {m: [t] for m in MODS}
GROUPS = {   # 프로토콜 짝짓기 표 (용어 순서 = 우선순위)
    "myeloid": {"monocyte": {"CAGE": ["CL:0002057"], "DNASE": ["CL:0001054"], "H3K27ac": ["CL:0001054"]}},
    "CNS_tissue": {"spinal cord": same("UBERON:0002240"),
                   "hippocampus": {"CAGE": ["UBERON:0001954"], "DNASE": ["UBERON:0001954"], "H3K27ac": ["UBERON:0002305"]},
                   "caudate": {"CAGE": ["UBERON:0001873"], "DNASE": ["UBERON:0002626"], "H3K27ac": ["UBERON:0001873"]},
                   "dlPFC": {"CAGE": ["UBERON:0002702"], "DNASE": ["UBERON:0009834"], "H3K27ac": ["UBERON:0009834"]}},
    "stromal": {"dermal fibroblast": same("CL:0002551"), "lung fibroblast": same("CL:0002553"), "MSC": same("CL:0000134")},
    "liver": {"hepatocyte": same("CL:0000182"), "liver tissue": same("UBERON:0002107")},
    "lymphoid": {"NK": same("CL:0000623"), "naive CD4 T": same("CL:0000895"),
                 "Treg": {"CAGE": ["CL:0002677"], "DNASE": ["CL:0000792"], "H3K27ac": ["CL:0000792"]},
                 "spleen": same("UBERON:0002106"), "thymus": same("UBERON:0002370")}}
TL = pd.read_csv(C3/"track_lineage.csv"); TL["col"] = TL.scorer + "::" + TL.key
TL["term"] = TL.key.map(lambda k: (re.search(r"(CL|UBERON|EFO|NTR|CLO):\d+", k) or [None])[0])
TL["lib"] = TL.key.str.replace(r"\|[+\-.]$", "", regex=True)   # CAGE 라이브러리 = 가닥 표시를 뗀 키
sel, CELLS = [], {}
for g, gs in GROUPS.items():
    for m in MODS:
        cols = []
        for gname, terms in gs.items():
            cand = TL[(TL.lineage == g) & (TL["mod"] == m) & TL.term.isin(terms[m])].copy()
            if cand.empty: raise SystemExit(f"{g} {gname} {m}: 후보 없음 — 짝짓기 표와 트랙 목록이 맞지 않는다")
            cand["o_term"] = cand.term.map(terms[m].index); cand["o_lq"] = cand.key.str.startswith("LQ").astype(int)
            first = cand.sort_values(["o_term", "o_lq", "key"]).iloc[0]
            unit = cand[cand.lib == first.lib] if m == "CAGE" else cand[cand.key == first.key]
            cols += unit.col.tolist()
            sel.append(dict(lineage=g, group=gname, layer=m, unit=first.lib if m == "CAGE" else first.key, biosample=first.biosample_name, n_tracks=len(unit), n_candidates=len(cand)))
        CELLS[(g, m)] = cols
S = pd.DataFrame(sel)
P("선택 단위 (묶음 × 층위마다 하나)\n" + S.to_string(index=False))
P("\n층위별 단위 수 (계통 × 층위)\n" + S.groupby(["lineage", "layer"]).size().unstack().to_string())
for g in LING: assert S[S.lineage == g].groupby("layer").size().nunique() == 1   # 세 층위 단위 수가 같아야 한다
FULL = {(g, m): TL[(TL.lineage == g) & (TL["mod"] == m)].col.tolist() for g in LING for m in MODS}
need = sorted(set(c for v in list(CELLS.values()) + list(FULL.values()) for c in v)); parts = sorted(SRC.glob("part_*.parquet")); P(f"트랙 열 {len(need)}개, 조각 {len(parts)}개 읽기")
V = []
for f in parts:
    have = set(pq.read_schema(f).names); cols = [c for c in need if c in have]
    D = pd.read_parquet(f, columns=["atlas_variant"] + cols).set_index("atlas_variant").astype(np.float32)
    def mx(v): return D[[c for c in v if c in have]].max(1) if any(c in have for c in v) else np.nan
    V.append(pd.DataFrame({**{f"M_{g}_{m}": mx(v) for (g, m), v in CELLS.items()}, **{f"F_{g}_{m}": mx(v) for (g, m), v in FULL.items()}}, index=D.index))
V = pd.concat(V); V = V[~V.index.duplicated()]
T = pd.read_parquet(X/"targets.parquet")[["rsid", "atlas_variant"]]; T = T.join(V, on="atlas_variant")
P(f"대상 {len(T):,}, Atlas 값 있음 {int(T[f'M_{LING[0]}_CAGE'].notna().sum()):,}")
R = T[[c for c in T.columns if c.startswith(("M_", "F_"))]].rank(pct=True)
ann = pd.DataFrame({"rsid": T.rsid}); chk = pd.DataFrame({"rsid": T.rsid})
for g in LING:
    for m in MODS: ann[f"M5_{g}_{m}"] = (R[f"M_{g}_{m}"] >= .95).astype(np.int8); chk[f"L5_{g}_{m}"] = (R[f"F_{g}_{m}"] >= .95).astype(np.int8)
ann = ann.drop_duplicates("rsid").set_index("rsid"); chk = chk.drop_duplicates("rsid").set_index("rsid")
NAMES = list(ann.columns)
if not TEST:
    ref = pd.read_csv(Y/"annot_sd_common.csv", index_col=0).sd_common; mine = chk.std(ddof=0); dif = (mine - ref[mine.index]).abs()
    P(f"\n확인: 다시 만든 전체 트랙 상위 5 % sd 대 3d L5 sd — 최대 차이 {dif.max():.2e}")
    if dif.max() >= 1e-9: raise SystemExit("상위 5 % 재현 실패 — 트랙 묶음 또는 백분위 방식이 3d와 다르다. 멈춘다.")
def jac(a, b): return (a & b).sum() / max((a | b).sum(), 1)
P("\n주석 크기와 전체 트랙 주석(L5)과의 Jaccard")
for g in LING: P(f"  {g:11s} " + "  ".join(f"{m} {int(ann[f'M5_{g}_{m}'].sum()):,} (J {jac(ann[f'M5_{g}_{m}'], chk[f'L5_{g}_{m}']):.2f})" for m in MODS))
P("\n같은 계통 안 층위 간 Jaccard (CAGE–DNase / CAGE–H3K27ac / DNase–H3K27ac): 일치 주석 | 전체 트랙")
for g in LING:
    P(f"  {g:11s} " + " / ".join(f"{jac(ann[f'M5_{g}_{a}'], ann[f'M5_{g}_{b}']):.2f}" for a, b in [("CAGE", "DNASE"), ("CAGE", "H3K27ac"), ("DNASE", "H3K27ac")])
      + "  |  " + " / ".join(f"{jac(chk[f'L5_{g}_{a}'], chk[f'L5_{g}_{b}']):.2f}" for a, b in [("CAGE", "DNASE"), ("CAGE", "H3K27ac"), ("DNASE", "H3K27ac")]))
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
groups = {f"M3_{g}": [f"M5_{g}_{m}" for m in MODS] for g in LING}
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
ann.std(ddof=0).rename("sd_common").to_csv(I/("test_annot_sd_common.csv" if TEST else "annot_sd_common.csv"))
S.to_csv(I/("test_selected_units.csv" if TEST else "selected_units.csv"), index=False)
(EXP/"results"/("i01_annot_test.txt" if TEST else "i01_annot.txt")).write_text("\n".join(log), encoding="utf-8")
