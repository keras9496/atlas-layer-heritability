"""Repository-wide paths.

DATA is the root of the data folder described in data_layout.md. Every script builds its paths from it, so a
clone on another machine only needs the environment variable AGA_DATA (for example AGA_DATA=/mnt/data/atlas).
Shell runners use AGA_DATA_POSIX for the same folder in POSIX form (Git Bash on Windows: /d/... for D:/...).
"""
import os
from pathlib import Path

DATA = os.environ.get("AGA_DATA", "D:/AD_GWAS_Atlas_data").replace("\\", "/").rstrip("/")
REPO = Path(__file__).resolve().parent
EXPERIMENTS = REPO/"experiments"

# Analysis stages in the order they were run (folder name, manuscript name of the stage, original working-folder name)
STAGES = [
    ("e01_targets_sldsc", "Target SNPs, initial AD analysis and AD robustness (3, 3b)", "실험3_SLDSC"),
    ("e02_lineage_map", "Lineage map (3c)", "실험3c_형질세포지도"),
    ("e03_single_layer_map", "Single-layer map and eQTL comparison (3d, 3d-A)", "실험3d_모달리티지도"),
    ("e04_measured_annotations", "Measured-annotation conditioning (3e)", "실험3e_실측주석조건부"),
    ("e05_layer_map", "Layer map with resampling and validation (3f)", "실험3f_층위이질성"),
    ("e06_thresholds", "Threshold sensitivity (3g)", "실험3g_역치민감도"),
    ("e07_shared_lineages", "Comparison within shared lineages (3h)", "실험3h_계통고정비교"),
    ("e08_cell_type_matched", "Cell-type-matched layers (3i)", "실험3i_세포형일치"),
    ("e09_joint_jackknife", "Joint jackknife across traits (3j)", "실험3j_형질간공분산"),
    ("e10_proximal_distal", "Proximal and distal SNPs (3k)", "실험3k_근위원위"),
    ("e11_activity", "Activity conditioning (3l)", "실험3l_활성도조건부"),
    ("e12_diagnostics", "Diagnostics and influential regions (3m)", "실험3m_기존산출물점검"),
    ("e13_larger_gwas", "Larger GWAS for five traits (3n)", "실험3n_최신GWAS"),
    ("e14_continuous_bins", "Continuous and non-overlapping annotations (3o)", "실험3o_연속구간"),
    ("e15_activity_constructions", "Activity constructions (3p)", "실험3p_활성도구성"),
    ("e16_within_active", "Effects within predicted active elements (3q)", "실험3q_활성요소내효과"),
    ("e17_larger_ad_gwas", "Larger AD GWAS (3r)", "실험3r_더큰AD_GWAS"),
    ("e18_cell_type_eqtl", "Cell-type eQTL conditioning (3s)", "실험3s_세포형eQTL조건화"),
]
