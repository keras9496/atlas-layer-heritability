# -*- coding: utf-8 -*-
"""그림 공용 모듈: 결과 파일의 자리, 색, 이름. 그림 스크립트는 experiments/*/results 의 결과 파일만 읽는다(새 적합 없음).
 층위 색(파랑·주황·초록)과 형질 범주 색 다섯은 색각 이상 대비와 명도 검사를 거친 값이다.
"""
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = Path(__file__).resolve().parents[1] / "experiments"      # 분석 단계의 폴더(결과 파일은 <단계>/results)
plt.rcParams.update({"font.family": "Arial", "font.size": 8.5, "axes.titlesize": 9.5, "axes.labelsize": 8.5,
                     "axes.unicode_minus": False, "svg.fonttype": "none", "pdf.fonttype": 42})

# ------------------------------------------------------------------ 색
COL = {"CAGE": "#2f66b0", "DNASE": "#c25e12", "H3K27ac": "#3b9466"}          # 층위
COL_OTHER = {"H3K4me1": "#8fb9a0", "H3K4me3": "#5f7f6c", "RNA": "#9a9a9a"}   # 모식도에서만 쓰는 나머지 층위(글자 표지와 함께)
SHORT = {"CAGE": "CAGE", "DNASE": "DNase", "H3K27ac": "H3K27ac"}
LAY = ["CAGE", "DNASE", "H3K27ac"]
CATS = ["Brain", "AD", "Blood & immune", "Lipid", "Skeletal"]
CATCOL = dict(zip(CATS, ["#7b4fb0", "#b8323f", "#00918f", "#b08a12", "#c2408f"]))   # 형질 범주
INK, MUTED, GRID = "#222222", "#666666", "#d0d0d0"

# ------------------------------------------------------------------ 이름
TN = {"AD": "Alzheimer's disease", "MONO": "Monocyte count", "EA": "Educational attainment", "SCZ": "Schizophrenia", "BIP": "Bipolar disorder",
      "NEURO": "Neuroticism", "BMI": "BMI", "IBD": "IBD", "RA": "Rheumatoid arthritis", "LDL": "LDL cholesterol", "HEIGHT": "Height",
      "BMD": "Bone mineral density", "NEUT": "Neutrophil count", "TC": "Total cholesterol"}
LN = {"CNS_tissue": "Brain tissue", "neural_cells": "Neural cells", "myeloid": "Myeloid", "lymphoid": "Lymphoid", "liver": "Liver", "GI": "Gastrointestinal",
      "cardiovascular": "Cardiovascular", "skeletal_muscle": "Skeletal muscle", "adipose": "Adipose", "stromal": "Stromal/mesenchymal",
      "reproductive_endocrine": "Reproductive/endocrine", "epithelial_barrier": "Epithelial/barrier"}
LIN12 = list(LN)
