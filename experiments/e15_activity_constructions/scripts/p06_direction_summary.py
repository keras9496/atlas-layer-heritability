"""실험 3p p06: CAGE − DNase 방향의 종합 정리 (기술, 새 계산 없음).
 기준: 3o 주석 기준 5개(연속, 구간 4개)와 3p 활성도 조건화 4개(EA, EAc, S, EAM; 새 GWAS는 EAM 없음). 비교 기준 = 전체 트랙 5 % 모형(E)의 c1 부호.
출력: results/p06_direction_summary.csv
"""
import pandas as pd, numpy as np
from pathlib import Path
R = Path(__file__).resolve().parents[2]
O = pd.read_csv(R/"e14_continuous_bins/results/o03_contrasts_long.csv"); C = pd.read_csv(R/"e15_activity_constructions/results/p04_contrasts.csv")
rows = []
for (t, g), d in C.groupby(["trait", "gwas"]):
    ref = float(d.set_index("model").loc["E", "c1"]); sgn = np.sign(ref)
    o = O[(O.trait == t) & (O.gwas == g) & O.criterion.isin(["CONT", "B1", "B2", "B3", "B4"])].set_index("criterion").c1
    a = d.set_index("model").c1.drop("E").reindex(["EA", "EAc", "S", "SA", "EAM"]).dropna()
    crit = int((np.sign(o) == sgn).sum()); act = {m: bool(np.sign(v) == sgn) for m, v in a.drop("SA", errors="ignore").items()}
    rows.append(dict(trait=t, gwas=g, ref_c1=round(ref, 2), criteria_same=f"{crit}/{len(o)}", **{f"{m}_same": v for m, v in act.items()},
                     **{f"{m}_c1": round(float(v), 2) for m, v in a.items()}, total_same=f"{crit + sum(act.values())}/{len(o) + len(act)}"))
D = pd.DataFrame(rows); D.to_csv(R/"e15_activity_constructions/results/p06_direction_summary.csv", index=False)
print(D[["trait", "gwas", "ref_c1", "criteria_same", "EA_c1", "EAc_c1", "S_c1", "EAM_c1", "total_same"]].to_string(index=False))
