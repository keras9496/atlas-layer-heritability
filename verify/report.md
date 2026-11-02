# Reproduction check

Reference: `atlas-layer-heritability/experiments`  
Reproduction: `atlas-layer-heritability/repro_run/experiments`  
Tolerance: rtol 1e-06, atol 1e-09

| Status | Files |
|---|---|
| same bytes | 67 |
| match | 66 |
| match (timing/path lines only) | 33 |
| note written by hand (not a pipeline output) | 17 |
| match (manuscript copy has added notes) | 4 |
| differ (expected, no value differs) | 2 |
| output of a script not included (AlphaGenome API) | 1 |

| File | Status | Detail |
|---|---|---|
| e01_targets_sldsc/results/RESULTS_EXP3_2026-09-27.md | match (manuscript copy has added notes) | paired with RESULTS_EXP3_2026-10-30.md; 179 reproduced lines all present in the manuscript's 227 |
| e01_targets_sldsc/results/RESULTS_EXP3b_2026-09-28.md | note written by hand (not a pipeline output) |  |
| e01_targets_sldsc/results/x01_targets.txt | match | 2 lines |
| e01_targets_sldsc/results/x03_annot.txt | match | 30 lines |
| e01_targets_sldsc/results/x04_sldsc.csv | same bytes |  |
| e01_targets_sldsc/results/x04_sldsc.txt | match | 94 lines |
| e01_targets_sldsc/results/x06_prepare_3b.txt | differ (expected, no value differs) | the manuscript's copy is the 3-line log of a partial re-run of x06 (only its last block); a full run writes all blocks, with the same values; 7 changed lines of 3 (after masking timestamps/durations/paths) |
| e01_targets_sldsc/results/x07_3b.csv | same bytes |  |
| e01_targets_sldsc/results/x07_3b.txt | match | 66 lines |
| e01_targets_sldsc/results/x07_3b_block_loo.csv | same bytes |  |
| e01_targets_sldsc/results/x07_3b_diffs.csv | same bytes |  |
| e02_lineage_map/results/RESULTS_EXP3c_2026-10-02.md | match (manuscript copy has added notes) | paired with RESULTS_EXP3c_2026-10-30.md; 106 reproduced lines all present in the manuscript's 165 |
| e02_lineage_map/results/z01_lineages.txt | output of a script not included (AlphaGenome API) |  |
| e02_lineage_map/results/z03_annot.txt | match | 41 lines |
| e02_lineage_map/results/z04_sldsc.csv | same bytes |  |
| e02_lineage_map/results/z04_sldsc.txt | match (timing/path lines only) | 92 lines |
| e02_lineage_map/results/z05_map_tau_star.csv | same bytes |  |
| e03_single_layer_map/results/RESULTS_EXP3d_2026-10-03.md | match (manuscript copy has added notes) | paired with RESULTS_EXP3d_2026-10-31.md; 99 reproduced lines all present in the manuscript's 147 |
| e03_single_layer_map/results/RESULTS_EXP3d_2026-10-04.md | match (manuscript copy has added notes) | paired with RESULTS_EXP3d_2026-11-02.md; 130 reproduced lines all present in the manuscript's 227 |
| e03_single_layer_map/results/RESULTS_EXP3d_A_eqtl_2026-10-04.md | note written by hand (not a pipeline output) |  |
| e03_single_layer_map/results/y03_annot.txt | match | 53 lines |
| e03_single_layer_map/results/y04_A.csv | same bytes |  |
| e03_single_layer_map/results/y04_A.txt | match (timing/path lines only) | 82 lines |
| e03_single_layer_map/results/y04_A_diffs.csv | same bytes |  |
| e03_single_layer_map/results/y04_B.csv | same bytes |  |
| e03_single_layer_map/results/y04_B.txt | match (timing/path lines only) | 26 lines |
| e03_single_layer_map/results/y05_D0_tau_star.csv | same bytes |  |
| e03_single_layer_map/results/y05_D1_map_tau_star.csv | same bytes |  |
| e03_single_layer_map/results/y05_D2_profiles.csv | same bytes |  |
| e03_single_layer_map/results/y06_eqtl_compare.csv | same bytes |  |
| e03_single_layer_map/results/y06_eqtl_compare.txt | match | 63 lines |
| e04_measured_annotations/results/RESULTS_EXP3e_2026-10-08.md | note written by hand (not a pipeline output) |  |
| e04_measured_annotations/results/w01_sources.txt | match | 28 lines |
| e04_measured_annotations/results/w02_annot.txt | match | 51 lines |
| e04_measured_annotations/results/w02_annot_test.txt | match | 30 lines |
| e04_measured_annotations/results/w03_diffs.csv | same bytes |  |
| e04_measured_annotations/results/w03_diffs_dryrun.csv | same bytes |  |
| e04_measured_annotations/results/w03_sldsc.csv | same bytes |  |
| e04_measured_annotations/results/w03_sldsc.txt | match (timing/path lines only) | 377 lines |
| e04_measured_annotations/results/w03_sldsc_dryrun.csv | same bytes |  |
| e04_measured_annotations/results/w03_sldsc_dryrun.txt | match (timing/path lines only) | 8 lines |
| e04_measured_annotations/results/w04_EH1.csv | same bytes |  |
| e04_measured_annotations/results/w04_Sm.csv | same bytes |  |
| e04_measured_annotations/results/w04_summary.txt | match | 58 lines |
| e04_measured_annotations/results/w05_own_layer.csv | same bytes |  |
| e04_measured_annotations/results/w05_own_layer.txt | match | 19 lines |
| e05_layer_map/results/RESULTS_EXP3f_2026-10-10.md | note written by hand (not a pipeline output) |  |
| e05_layer_map/results/v01_annot.txt | match | 48 lines |
| e05_layer_map/results/v01_annot_test.txt | match | 27 lines |
| e05_layer_map/results/v02_dry.txt | match (timing/path lines only) | 4 lines |
| e05_layer_map/results/v02_main.txt | match (timing/path lines only) | 232 lines |
| e05_layer_map/results/v02_val.txt | match (timing/path lines only) | 100 lines |
| e05_layer_map/results/v03_dry_layers.csv | same bytes |  |
| e05_layer_map/results/v03_dry_tests.txt | match | 6 lines |
| e05_layer_map/results/v03_main_layers.csv | same bytes |  |
| e05_layer_map/results/v03_main_tests.txt | match | 41 lines |
| e05_layer_map/results/v03_val_layers.csv | same bytes |  |
| e05_layer_map/results/v03_val_tests.txt | match | 27 lines |
| e05_layer_map/results/v04_validation_sumstats.txt | match (timing/path lines only) | 9 lines |
| e05_layer_map/results/v05_posthoc_loo.txt | match | 15 lines |
| e06_thresholds/results/RESULTS_EXP3g_2026-10-11.md | note written by hand (not a pipeline output) |  |
| e06_thresholds/results/u01_annot.txt | match | 50 lines |
| e06_thresholds/results/u01_annot_test.txt | match | 27 lines |
| e06_thresholds/results/u02_dry.txt | match (timing/path lines only) | 4 lines |
| e06_thresholds/results/u02_main.txt | match (timing/path lines only) | 43 lines |
| e06_thresholds/results/u03_dry_layers.csv | same bytes |  |
| e06_thresholds/results/u03_dry_tests.txt | match | 6 lines |
| e06_thresholds/results/u03_main_layers.csv | same bytes |  |
| e06_thresholds/results/u03_main_tests.txt | match | 132 lines |
| e07_shared_lineages/results/h01_dry.txt | match | 2 lines |
| e07_shared_lineages/results/h01_main.txt | match (timing/path lines only) | 48 lines |
| e07_shared_lineages/results/h02_dry_profiles.csv | same bytes |  |
| e07_shared_lineages/results/h02_dry_tests.txt | match | 2 lines |
| e07_shared_lineages/results/h02_main_profiles.csv | same bytes |  |
| e07_shared_lineages/results/h02_main_tests.txt | match | 86 lines |
| e07_shared_lineages/results/RESULTS_EXP3h_2026-10-13.md | note written by hand (not a pipeline output) |  |
| e08_cell_type_matched/results/i01_annot.txt | match | 96 lines |
| e08_cell_type_matched/results/i01_annot_test.txt | match | 73 lines |
| e08_cell_type_matched/results/i02_dry.txt | match | 2 lines |
| e08_cell_type_matched/results/i02_main.txt | match (timing/path lines only) | 15 lines |
| e08_cell_type_matched/results/i03_dry_layers.csv | same bytes |  |
| e08_cell_type_matched/results/i03_dry_tests.txt | match | 2 lines |
| e08_cell_type_matched/results/i03_main_layers.csv | same bytes |  |
| e08_cell_type_matched/results/i03_main_tests.txt | match | 32 lines |
| e08_cell_type_matched/results/RESULTS_EXP3i_2026-10-13.md | note written by hand (not a pipeline output) |  |
| e09_joint_jackknife/results/j01_blocks.txt | match | 3 lines |
| e09_joint_jackknife/results/j02_check.txt | match (timing/path lines only) | 5 lines |
| e09_joint_jackknife/results/j02_dry.txt | match (timing/path lines only) | 3 lines |
| e09_joint_jackknife/results/j02_main.txt | match (timing/path lines only) | 383 lines |
| e09_joint_jackknife/results/j03_dry_tests.txt | match | 3 lines |
| e09_joint_jackknife/results/j03_main_summary.csv | same bytes |  |
| e09_joint_jackknife/results/j03_main_tests.txt | match | 65 lines |
| e09_joint_jackknife/results/RESULTS_EXP3j_2026-10-14.md | note written by hand (not a pipeline output) |  |
| e10_proximal_distal/results/k01_annot.txt | match | 32 lines |
| e10_proximal_distal/results/k01_annot_test.txt | match | 10 lines |
| e10_proximal_distal/results/k02_dry.txt | match | 1 lines |
| e10_proximal_distal/results/k02_main.txt | match (timing/path lines only) | 14 lines |
| e10_proximal_distal/results/k03_dry_layers.csv | same bytes |  |
| e10_proximal_distal/results/k03_dry_tests.txt | match | 4 lines |
| e10_proximal_distal/results/k03_main_layers.csv | same bytes |  |
| e10_proximal_distal/results/k03_main_tests.txt | match | 50 lines |
| e10_proximal_distal/results/RESULTS_EXP3k_2026-10-15.md | note written by hand (not a pipeline output) |  |
| e11_activity/results/l02_check.txt | match | 6 lines |
| e11_activity/results/l02_main.txt | match | 35 lines |
| e11_activity/results/l02_raw_scale.csv | same bytes |  |
| e11_activity/results/l03_dry.txt | match (timing/path lines only) | 1 lines |
| e11_activity/results/l03_main.txt | match (timing/path lines only) | 42 lines |
| e11_activity/results/l03_pretest.txt | match (timing/path lines only) | 3 lines |
| e11_activity/results/l04_dry_tests.txt | match | 2 lines |
| e11_activity/results/l04_main_summary.csv | same bytes |  |
| e11_activity/results/l04_main_tests.txt | match | 62 lines |
| e11_activity/results/l04_main_vif.csv | same bytes |  |
| e11_activity/results/l04_pretest_tests.txt | match | 4 lines |
| e11_activity/results/l05_values.csv | same bytes |  |
| e11_activity/results/RESULTS_EXP3l_2026-10-19.md | note written by hand (not a pipeline output) |  |
| e12_diagnostics/results/m01_trait_info.csv | same bytes |  |
| e12_diagnostics/results/m02_tests.txt | match | 52 lines |
| e12_diagnostics/results/m02_threshold_contrasts.csv | same bytes |  |
| e12_diagnostics/results/m03_coef_corr.csv | same bytes |  |
| e12_diagnostics/results/m03_collinearity.csv | same bytes |  |
| e12_diagnostics/results/m03_single_vs_joint.csv | same bytes |  |
| e12_diagnostics/results/m03_tests.txt | match | 77 lines |
| e12_diagnostics/results/m04_exclusion_sets.json | match | 65 lines |
| e12_diagnostics/results/m04a_lobo.csv | same bytes |  |
| e12_diagnostics/results/m04a_tests.txt | match | 33 lines |
| e12_diagnostics/results/m04a_top_blocks.csv | same bytes |  |
| e12_diagnostics/results/m04b_check.txt | match (timing/path lines only) | 2 lines |
| e12_diagnostics/results/m04b_summary.csv | same bytes |  |
| e12_diagnostics/results/m04b_tests.txt | match (timing/path lines only) | 44 lines |
| e12_diagnostics/results/m05_convergence.csv | same bytes |  |
| e12_diagnostics/results/m05_tests.txt | match | 16 lines |
| e12_diagnostics/results/m06_decomp.csv | same bytes |  |
| e12_diagnostics/results/RESULTS_EXP3m_2026-10-17.md | note written by hand (not a pipeline output) |  |
| e13_larger_gwas/results/n01_munge.txt | differ (expected, no value differs) | the manuscript's copy repeats the LDSC munge summary on the same line (the summary was appended twice by the original run); the values are identical; 4 changed lines of 15 (after masking timestamps/durations/paths) |
| e13_larger_gwas/results/n02_check.txt | match (timing/path lines only) | 4 lines |
| e13_larger_gwas/results/n02_main.txt | match (timing/path lines only) | 170 lines |
| e13_larger_gwas/results/n03_activity.csv | same bytes |  |
| e13_larger_gwas/results/n03_contrasts.csv | same bytes |  |
| e13_larger_gwas/results/n03_layers.csv | same bytes |  |
| e13_larger_gwas/results/n03_lineage.csv | same bytes |  |
| e13_larger_gwas/results/n03_tests.txt | match | 42 lines |
| e13_larger_gwas/results/n03_trait_info.csv | same bytes |  |
| e13_larger_gwas/results/RESULTS_EXP3n_2026-10-21.md | note written by hand (not a pipeline output) |  |
| e14_continuous_bins/results/o01_annot.txt | match | 49 lines |
| e14_continuous_bins/results/o02_dry.txt | match | 3 lines |
| e14_continuous_bins/results/o02_main.txt | match (timing/path lines only) | 38 lines |
| e14_continuous_bins/results/o03_contrasts_long.csv | same bytes |  |
| e14_continuous_bins/results/o03_pairwise.csv | same bytes |  |
| e14_continuous_bins/results/o03_stability.csv | same bytes |  |
| e14_continuous_bins/results/o03_tests.txt | match | 81 lines |
| e14_continuous_bins/results/o03_vif.csv | same bytes |  |
| e14_continuous_bins/results/RESULTS_EXP3o_2026-10-21.md | note written by hand (not a pipeline output) |  |
| e15_activity_constructions/results/p01_check.txt | match (timing/path lines only) | 52 lines |
| e15_activity_constructions/results/p01_main.txt | match | 22 lines |
| e15_activity_constructions/results/p02_dry.txt | match | 4 lines |
| e15_activity_constructions/results/p02_main.txt | match (timing/path lines only) | 57 lines |
| e15_activity_constructions/results/p03_dry_tests.txt | match | 4 lines |
| e15_activity_constructions/results/p03_main_summary.csv | same bytes |  |
| e15_activity_constructions/results/p03_main_tests.txt | match | 137 lines |
| e15_activity_constructions/results/p03_main_vif.csv | same bytes |  |
| e15_activity_constructions/results/p04_contrasts.csv | same bytes |  |
| e15_activity_constructions/results/p04_values.csv | same bytes |  |
| e15_activity_constructions/results/p06_direction_summary.csv | same bytes |  |
| e15_activity_constructions/results/RESULTS_EXP3p_2026-10-23.md | note written by hand (not a pipeline output) |  |
| e16_within_active/results/q01_check.txt | match (timing/path lines only) | 48 lines |
| e16_within_active/results/q01_main.txt | match | 22 lines |
| e16_within_active/results/q02_dry.txt | match | 3 lines |
| e16_within_active/results/q02_main.txt | match (timing/path lines only) | 38 lines |
| e16_within_active/results/q03_dry_tests.txt | match | 3 lines |
| e16_within_active/results/q03_main_summary.csv | same bytes |  |
| e16_within_active/results/q03_main_tests.txt | match | 64 lines |
| e16_within_active/results/q03_main_vif.csv | same bytes |  |
| e16_within_active/results/q03_values.csv | same bytes |  |
| e16_within_active/results/RESULTS_EXP3q_2026-10-23.md | note written by hand (not a pipeline output) |  |
| e17_larger_ad_gwas/results/r01_munge.txt | match (timing/path lines only) | 3 lines |
| e17_larger_ad_gwas/results/r02_check.txt | match (timing/path lines only) | 5 lines |
| e17_larger_ad_gwas/results/r02_main.txt | match (timing/path lines only) | 40 lines |
| e17_larger_ad_gwas/results/r03_summary.csv | same bytes |  |
| e17_larger_ad_gwas/results/r03_tests.txt | match | 32 lines |
| e17_larger_ad_gwas/results/RESULTS_EXP3r_2026-10-25.md | note written by hand (not a pipeline output) |  |
| e18_cell_type_eqtl/results/RESULTS_EXP3s_2026-10-25.md | note written by hand (not a pipeline output) |  |
| e18_cell_type_eqtl/results/s01_check.txt | match | 111 lines |
| e18_cell_type_eqtl/results/s01_main.txt | match | 22 lines |
| e18_cell_type_eqtl/results/s02_dry.txt | match | 3 lines |
| e18_cell_type_eqtl/results/s02_main.txt | match (timing/path lines only) | 40 lines |
| e18_cell_type_eqtl/results/s03_dry_tests.txt | match | 3 lines |
| e18_cell_type_eqtl/results/s03_main_tests.txt | match | 75 lines |
| e18_cell_type_eqtl/results/s03_summary.csv | same bytes |  |
| e18_cell_type_eqtl/results/s03_values.csv | same bytes |  |
| e18_cell_type_eqtl/results/s03_vif.csv | same bytes |  |
