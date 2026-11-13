#!/usr/bin/env bash
# Run the whole analysis in order (Git Bash on Windows, or bash on Linux/macOS).
#
#   AGA_DATA=D:/AD_GWAS_Atlas_repro AGA_DATA_POSIX=/d/AD_GWAS_Atlas_repro bash run_all.sh            # everything, in order
#   bash run_all.sh list                                                                               # print the steps
#   bash run_all.sh e05_                                                                               # only steps whose name starts with e05_
#
# Each step is "name | requires | working dir | command". A step is skipped when its done-marker exists
# ($AGA_DATA/_status/<name>.done), waits while a required step has no marker, and is claimed with a lock file so that a
# second copy of this script, started in another terminal, takes the steps the first one has not claimed (the
# single-layer-map fits of e03 take two days and can run beside the rest; the fits of e19-e22 are split into pieces for the
# same purpose). Logs: $AGA_DATA/_logs/<name>.log. A runner that is killed leaves the lock file of its step behind: delete
# $AGA_DATA/_status/*.lock before starting again.
# The Atlas query scripts (x02, z01, z02, l01) are not in this repository; their outputs are inputs (data_layout.md).
set -u
cd "$(dirname "$0")"
export PYTHONIOENCODING=utf-8
export AGA_DATA="${AGA_DATA:-D:/AD_GWAS_Atlas_data}"
export AGA_DATA_POSIX="${AGA_DATA_POSIX:-$(printf '%s' "$AGA_DATA" | sed -E 's#^([A-Za-z]):#/\L\1#')}"
ST="$AGA_DATA_POSIX/_status"; LG="$AGA_DATA_POSIX/_logs"; mkdir -p "$ST" "$LG"
E=experiments; P="python -W ignore -u"

STEPS=(
# ---- e01: target SNPs, AD summary statistics, initial AD analysis (3), AD robustness (3b)
"e01_x00_sumstats      |                         | $E/e01_targets_sldsc        | $P scripts/x00_prepare_sumstats.py"
"e01_x01_targets       | e01_x00_sumstats        | $E/e01_targets_sldsc        | $P scripts/x01_targets.py"
"e01_x03_annot         | e01_x01_targets         | $E/e01_targets_sldsc        | $P scripts/x03_annot_ldscores.py"
"e01_x04_sldsc         | e01_x03_annot           | $E/e01_targets_sldsc        | $P scripts/x04_sldsc.py"
"e01_x05_summary       | e01_x04_sldsc           | $E/e01_targets_sldsc        | $P scripts/x05_summary.py"
"e01_x06_prepare_3b    | e01_x04_sldsc           | $E/e01_targets_sldsc        | $P scripts/x06_prepare_3b.py"
"e01_x07_sldsc_3b      | e01_x06_prepare_3b      | $E/e01_targets_sldsc        | $P scripts/x07_sldsc_3b.py"
# ---- e02: lineage map (3c)
"e02_z03_annot         | e01_x03_annot           | $E/e02_lineage_map          | $P scripts/z03_annot_ldscores.py"
"e02_z04_sldsc         | e02_z03_annot e01_x06_prepare_3b | $E/e02_lineage_map | $P scripts/z04_sldsc.py"
"e02_z05_summary       | e02_z04_sldsc           | $E/e02_lineage_map          | $P scripts/z05_summary.py"
# ---- e03: lineage x layer annotations (needed by most later stages), single-layer map fits (long), eQTL comparison (3d, 3d-A)
"e03_y03_annot         | e02_z03_annot           | $E/e03_single_layer_map     | $P scripts/y03_annot_ldscores.py"
"e03_y04_A             | e03_y03_annot e02_z04_sldsc | $E/e03_single_layer_map | $P scripts/y04_sldsc.py A"
"e03_y05_A             | e03_y04_A               | $E/e03_single_layer_map     | $P scripts/y05_summary.py"
"e03_y04_B             | e03_y05_A               | $E/e03_single_layer_map     | $P scripts/y04_sldsc.py B"
"e03_y05_B             | e03_y04_B               | $E/e03_single_layer_map     | $P scripts/y05_summary.py"
"e03_y06_eqtl          | e03_y05_B               | $E/e03_single_layer_map     | $P scripts/y06_eqtl_compare.py"
# ---- e04: measured annotations (3e); w04 and w05 need the single-layer fits of e03
"e04_w01_sources       |                         | $E/e04_measured_annotations | $P scripts/w01_measured_sources.py"
"e04_w02_annot         | e04_w01_sources e03_y03_annot | $E/e04_measured_annotations | $P scripts/w02_annot_ldscores.py"
"e04_w02_test          | e04_w01_sources e03_y03_annot | $E/e04_measured_annotations | env W02_TEST=1 $P scripts/w02_annot_ldscores.py"
"e04_w03_dry           | e04_w02_annot           | $E/e04_measured_annotations | env W03_DRYRUN=1 $P scripts/w03_sldsc.py"
"e04_w03_sldsc         | e04_w03_dry             | $E/e04_measured_annotations | env W03_DRYRUN=0 $P scripts/w03_sldsc.py"
# ---- e05: layer map with resampling, validation traits, leave-one-trait-out (3f)
"e05_v01_resample      | e03_y03_annot           | $E/e05_layer_map            | $P scripts/v01_resample_annot.py"
"e05_v01_test          | e03_y03_annot           | $E/e05_layer_map            | env V01_TEST=1 $P scripts/v01_resample_annot.py"
"e05_v02_dry           | e05_v01_resample        | $E/e05_layer_map            | $P scripts/v02_sldsc.py dry"
"e05_v03_dry           | e05_v02_dry             | $E/e05_layer_map            | $P scripts/v03_tests.py dry"
"e05_v02_main          | e05_v03_dry             | $E/e05_layer_map            | $P scripts/v02_sldsc.py main"
"e05_v03_main          | e05_v02_main            | $E/e05_layer_map            | $P scripts/v03_tests.py main"
"e05_v04_val_sumstats  | e01_x00_sumstats        | $E/e05_layer_map            | $P scripts/v04_validation_sumstats.py"
"e05_v02_val           | e05_v03_main e05_v04_val_sumstats | $E/e05_layer_map  | $P scripts/v02_sldsc.py val"
"e05_v03_val           | e05_v02_val             | $E/e05_layer_map            | $P scripts/v03_tests.py val"
"e05_v05_loo           | e05_v03_val             | $E/e05_layer_map            | $P scripts/v05_posthoc_loo.py"
# ---- e04 continued: own-layer comparison (Fig. 5a) and summary
"e04_w05_own_layer     | e04_w03_sldsc e05_v03_val e03_y04_B | $E/e04_measured_annotations | $P scripts/w05_own_layer.py"
"e04_w04_summary       | e04_w05_own_layer e03_y06_eqtl | $E/e04_measured_annotations | $P scripts/w04_summary.py"
# ---- e06: thresholds (3g)
"e06_u01_annot         | e05_v03_main            | $E/e06_thresholds           | $P scripts/u01_threshold_annot.py"
"e06_u01_test          | e05_v03_main            | $E/e06_thresholds           | env U01_TEST=1 $P scripts/u01_threshold_annot.py"
"e06_u02_dry           | e06_u01_annot           | $E/e06_thresholds           | $P scripts/u02_sldsc.py dry"
"e06_u03_dry           | e06_u02_dry             | $E/e06_thresholds           | $P scripts/u03_tests.py dry"
"e06_u02_main          | e06_u03_dry             | $E/e06_thresholds           | $P scripts/u02_sldsc.py main"
"e06_u03_main          | e06_u02_main            | $E/e06_thresholds           | $P scripts/u03_tests.py main"
# ---- e07: shared lineages (3h); e08: cell-type-matched layers (3i)
"e07_h01_dry           | e05_v03_val             | $E/e07_shared_lineages      | $P scripts/h01_sldsc.py dry"
"e07_h02_dry           | e07_h01_dry             | $E/e07_shared_lineages      | $P scripts/h02_tests.py dry"
"e07_h01_sldsc         | e05_v03_val             | $E/e07_shared_lineages      | $P scripts/h01_sldsc.py main"
"e07_h02_tests         | e07_h01_sldsc           | $E/e07_shared_lineages      | $P scripts/h02_tests.py main"
"e08_i01_annot         | e05_v03_val             | $E/e08_cell_type_matched    | $P scripts/i01_matched_annot.py"
"e08_i01_test          | e05_v03_val             | $E/e08_cell_type_matched    | env I01_TEST=1 $P scripts/i01_matched_annot.py"
"e08_i02_dry           | e08_i01_annot           | $E/e08_cell_type_matched    | $P scripts/i02_sldsc.py dry"
"e08_i03_dry           | e08_i02_dry             | $E/e08_cell_type_matched    | $P scripts/i03_tests.py dry"
"e08_i02_main          | e08_i03_dry             | $E/e08_cell_type_matched    | $P scripts/i02_sldsc.py main"
"e08_i03_main          | e08_i02_main            | $E/e08_cell_type_matched    | $P scripts/i03_tests.py main"
# ---- e09: aligned blocks and joint jackknife (3j)
"e09_j01_blocks        | e05_v03_val             | $E/e09_joint_jackknife      | $P scripts/j01_blocks.py"
"e09_j02_check        | e09_j01_blocks e05_v03_main e06_u03_main e07_h02_tests e08_i03_main | $E/e09_joint_jackknife | $P scripts/j02_sldsc.py check"
"e09_j02_dry          | e09_j01_blocks e06_u03_main e08_i03_main | $E/e09_joint_jackknife | $P scripts/j02_sldsc.py dry"
"e09_j03_dry          | e09_j02_dry e07_h02_tests | $E/e09_joint_jackknife    | $P scripts/j03_tests.py dry"
"e09_j02_sldsc         | e09_j01_blocks e06_u03_main e08_i03_main | $E/e09_joint_jackknife | $P scripts/j02_sldsc.py main"
"e09_j03_tests         | e09_j02_sldsc e07_h02_tests | $E/e09_joint_jackknife  | $P scripts/j03_tests.py main"
# ---- e10: proximal and distal SNPs (3k)
"e10_k01_annot         | e09_j03_tests           | $E/e10_proximal_distal      | $P scripts/k01_split_annot.py"
"e10_k01_test          | e09_j03_tests           | $E/e10_proximal_distal      | env K01_TEST=1 $P scripts/k01_split_annot.py"
"e10_k02_dry           | e10_k01_annot           | $E/e10_proximal_distal      | $P scripts/k02_sldsc.py dry"
"e10_k03_dry           | e10_k02_dry             | $E/e10_proximal_distal      | $P scripts/k03_tests.py dry"
"e10_k02_sldsc         | e10_k01_annot           | $E/e10_proximal_distal      | $P scripts/k02_sldsc.py main"
"e10_k03_tests         | e10_k02_sldsc           | $E/e10_proximal_distal      | $P scripts/k03_tests.py main"
# ---- e11: predicted activity (3l)
"e11_l02_check        | e09_j03_tests           | $E/e11_activity             | $P scripts/l02_annot_ldscores.py check"
"e11_l02_annot         | e09_j03_tests           | $E/e11_activity             | $P scripts/l02_annot_ldscores.py main"
"e11_l03_pretest      | e11_l02_annot e08_i03_main | $E/e11_activity          | $P scripts/l03_sldsc.py pretest"
"e11_l04_pretest      | e11_l03_pretest         | $E/e11_activity             | $P scripts/l04_tests.py pretest"
"e11_l03_dry           | e11_l02_annot           | $E/e11_activity             | $P scripts/l03_sldsc.py dry"
"e11_l04_dry           | e11_l03_dry             | $E/e11_activity             | $P scripts/l04_tests.py dry"
"e11_l03_main          | e11_l04_dry             | $E/e11_activity             | $P scripts/l03_sldsc.py main"
"e11_l04_main          | e11_l03_main            | $E/e11_activity             | $P scripts/l04_tests.py main"
"e11_l05_export        | e11_l04_main            | $E/e11_activity             | $P scripts/l05_fig_export.py"
# ---- e12: diagnostics (3m)
"e12_m01_trait_info    | e11_l05_export          | $E/e12_diagnostics          | $P scripts/m01_trait_info.py"
"e12_m02_thresholds    | e11_l05_export          | $E/e12_diagnostics          | $P scripts/m02_threshold_contrasts.py"
"e12_m03_collinearity  | e11_l05_export e03_y04_B | $E/e12_diagnostics         | $P scripts/m03_collinearity.py"
"e12_m04a_influence    | e11_l05_export          | $E/e12_diagnostics          | $P scripts/m04a_block_influence.py"
"e12_m04b_check       | e12_m04a_influence      | $E/e12_diagnostics          | $P scripts/m04b_exclusion_refit.py check"
"e12_m04b_refit        | e12_m04a_influence      | $E/e12_diagnostics          | $P scripts/m04b_exclusion_refit.py main"
"e12_m05_convergence   | e11_l05_export          | $E/e12_diagnostics          | $P scripts/m05_resample_convergence.py"
"e12_m06_export        | e12_m04b_refit          | $E/e12_diagnostics          | $P scripts/m06_supp_export.py"
# ---- e13: larger GWAS (3n); e14: continuous and bins (3o)
"e13_n01_munge         | e01_x00_sumstats        | $E/e13_larger_gwas          | $P scripts/n01_munge.py"
"e13_n02_check        | e13_n01_munge e12_m06_export | $E/e13_larger_gwas     | $P scripts/n02_sldsc.py check"
"e13_n02_sldsc         | e13_n01_munge e12_m06_export | $E/e13_larger_gwas     | $P scripts/n02_sldsc.py main"
"e13_n03_tests         | e13_n02_sldsc           | $E/e13_larger_gwas          | $P scripts/n03_tests.py"
"e14_o01_annot         | e13_n03_tests           | $E/e14_continuous_bins      | $P scripts/o01_annot_ldscores.py"
"e14_o02_dry           | e14_o01_annot           | $E/e14_continuous_bins      | $P scripts/o02_sldsc.py dry"
"e14_o02_sldsc         | e14_o01_annot           | $E/e14_continuous_bins      | $P scripts/o02_sldsc.py main"
"e14_o03_tests         | e14_o02_sldsc           | $E/e14_continuous_bins      | $P scripts/o03_tests.py"
# ---- e15: activity constructions (3p); e16: within active elements (3q)
"e15_p01_check        | e14_o03_tests           | $E/e15_activity_constructions | $P scripts/p01_annot_ldscores.py check"
"e15_p01_annot         | e15_p01_check           | $E/e15_activity_constructions | $P scripts/p01_annot_ldscores.py main"
"e15_p02_dry           | e15_p01_annot           | $E/e15_activity_constructions | $P scripts/p02_sldsc.py dry"
"e15_p03_dry           | e15_p02_dry             | $E/e15_activity_constructions | $P scripts/p03_tests.py dry"
"e15_p02_main          | e15_p03_dry             | $E/e15_activity_constructions | $P scripts/p02_sldsc.py main"
"e15_p03_main          | e15_p02_main            | $E/e15_activity_constructions | $P scripts/p03_tests.py main"
"e15_p04_export        | e15_p03_main            | $E/e15_activity_constructions | $P scripts/p04_fig_export.py"
"e15_p06_direction     | e15_p04_export          | $E/e15_activity_constructions | $P scripts/p06_direction_summary.py"
"e16_q01_check        | e15_p06_direction       | $E/e16_within_active        | $P scripts/q01_annot_ldscores.py check"
"e16_q01_annot         | e16_q01_check           | $E/e16_within_active        | $P scripts/q01_annot_ldscores.py main"
"e16_q02_dry           | e16_q01_annot           | $E/e16_within_active        | $P scripts/q02_sldsc.py dry"
"e16_q03_dry           | e16_q02_dry             | $E/e16_within_active        | $P scripts/q03_tests.py dry"
"e16_q02_main          | e16_q03_dry             | $E/e16_within_active        | $P scripts/q02_sldsc.py main"
"e16_q03_main          | e16_q02_main            | $E/e16_within_active        | $P scripts/q03_tests.py main"
# ---- e17: larger AD GWAS (3r); e18: cell-type eQTLs (3s)
"e17_r01_munge         | e01_x00_sumstats        | $E/e17_larger_ad_gwas       | $P scripts/r01_munge.py"
"e17_r02_check        | e17_r01_munge e16_q03_main | $E/e17_larger_ad_gwas    | $P scripts/r02_sldsc.py check"
"e17_r02_sldsc         | e17_r01_munge e16_q03_main | $E/e17_larger_ad_gwas    | $P scripts/r02_sldsc.py main"
"e17_r03_tests         | e17_r02_sldsc           | $E/e17_larger_ad_gwas       | $P scripts/r03_tests.py"
"e18_s01_check        | e17_r03_tests           | $E/e18_cell_type_eqtl       | $P scripts/s01_annot_ldscores.py check"
"e18_s01_annot         | e18_s01_check           | $E/e18_cell_type_eqtl       | $P scripts/s01_annot_ldscores.py main"
"e18_s02_dry           | e18_s01_annot           | $E/e18_cell_type_eqtl       | $P scripts/s02_sldsc.py dry"
"e18_s03_dry           | e18_s02_dry             | $E/e18_cell_type_eqtl       | $P scripts/s03_tests.py dry"
"e18_s02_main          | e18_s03_dry             | $E/e18_cell_type_eqtl       | $P scripts/s02_sldsc.py main"
"e18_s03_main          | e18_s02_main            | $E/e18_cell_type_eqtl       | $P scripts/s03_tests.py main"
# ---- e19-e22 and the figures are appended below
)

# ---- e19-e22: benchmark of 14 traits and the panel of 106 GWAS. These stages read the annotations, LD scores and run files of
# the stages above and start when those are done (G). The S-LDSC fits are split into pieces (name_<k>of<K>: every K-th job,
# starting at job k), as they were run; start several copies of this script to run pieces side by side. The number of pieces is
# part of the name of the log that each piece writes to results/.
PU="python -X utf8 -W ignore -u"
# G: the steps above that no other step requires. When they are done, all of e01-e18 is done.
names=""; reqs=" "; for s in "${STEPS[@]}"; do IFS='|' read -r n r d c <<< "$s"; names+="$n "; reqs+="$r "; done
G=""; for n in $names; do [[ "$reqs" == *" $n "* ]] || G+="$n "; done; G=${G% }; unset names reqs
add() { STEPS+=("$1 | $2 | $3 | $4"); }                                                                       # name, requires, dir, command
pieces() { local k; for ((k = 0; k < $2; k++)); do add "$1_${k}of$2" "$3" "$4" "$5 $k $2"; done; }           # name, K, requires, dir, command
allof() { local k s=""; for ((k = 0; k < $2; k++)); do s+="$1_${k}of$2 "; done; printf '%s' "$s"; }           # the names of the K pieces
# ---- e19: heritability share and enrichment (t01); three-layer models of the 14 traits in every tissue group (t08) (3t)
D=$E/e19_enrichment_all_groups
add e19_t01_check    "$G"               $D "$PU scripts/t01_enrichment.py check"
add e19_t01_fit      "e19_t01_check"    $D "$PU scripts/t01_enrichment.py fit"
add e19_t01_extract  "e19_t01_fit"      $D "$PU scripts/t01_enrichment.py extract"
add e19_t08_fit      "$G"               $D "$PU scripts/t08_matrix.py fit"
# ---- e20: panel of 106 GWAS: trait table, tissue-group and all-track three-layer models, tissue-group map, resampled models (3u)
D=$E/e20_gwas_panel
add e20_u00_inputs   "$G"               $D "$PU scripts/u00_inputs.py"
add e20_u01_prep     "e20_u00_inputs"   $D "$PU scripts/u01_fit.py prep"
add e20_u01_check    "e20_u01_prep"     $D "$PU scripts/u01_fit.py check"
pieces e20_u01_fit 8 "e20_u01_check"    $D "$PU scripts/u01_fit.py fit"
add e20_u02_lineage  "$(allof e20_u01_fit 8)" $D "$PU scripts/u02_lineage.py"
pieces e20_u01_fitR 12 "e20_u02_lineage" $D "$PU scripts/u01_fit.py fitR"
# fitS lists 120 models that fitR also lists (six GWAS analysed in the stromal group). It starts when fitR is done and skips them,
# so that no model is fitted by two copies of this script at the same time.
pieces e20_u01_fitS 12 "$(allof e20_u01_fitR 12)" $D "$PU scripts/u01_fit.py fitS"
add e20_u03_pairs    "e20_u02_lineage"  $D "$PU scripts/u03_pairs.py"
# ---- e21: panel of 106 GWAS: 72 annotations fitted one at a time, types of signal (3v)
D=$E/e21_panel_single_layer_map
add e21_v01_prepld   "$G"               $D "$PU scripts/v01_fit72.py prepld"
add e21_v01_check    "e21_v01_prepld"   $D "$PU scripts/v01_fit72.py check"
pieces e21_v01_fit 10 "e21_v01_check e20_u01_prep" $D "$PU scripts/v01_fit72.py fit"
add e21_v02_types    "$(allof e21_v01_fit 10)" $D "$PU scripts/v02_types.py all"
# ---- e22: benchmark of 14 traits; enriched GWAS-tissue-group pairs of the panel (3w). The panel pieces fit the resampled models
# that e20 has not fitted, so they wait for all fitR and fitS pieces of e20.
D=$E/e22_benchmark14
add e22_w01_check    "$G"               $D "$PU scripts/w01_fit.py check"
pieces e22_w01_bench 10 "e22_w01_check" $D "$PU scripts/w01_fit.py fit bench"
pieces e22_w01_panel 10 "e22_w01_check e20_u02_lineage $(allof e20_u01_fitR 12) $(allof e20_u01_fitS 12)" $D "$PU scripts/w01_fit.py fit panel"
pieces e22_w01_panelpd 6 "e22_w01_check e20_u02_lineage" $D "$PU scripts/w01_fit.py fit panelpd"
add e22_w02_check    "e21_v01_prepld"   $D "$PU scripts/w02_fit72.py check"
pieces e22_w02_fit 5 "e22_w02_check"    $D "$PU scripts/w02_fit72.py fit"
add e22_w03_check    "$G"               $D "$PU scripts/w03_actbin.py check"
pieces e22_w03_main 6 "e22_w03_check"   $D "$PU scripts/w03_actbin.py main"
pieces e22_w01_actbin 7 "e22_w01_check $(allof e22_w03_main 6)" $D "$PU scripts/w01_fit.py fit actbin"
add e22_w04_bench    "$(allof e22_w01_bench 10) $(allof e22_w02_fit 5) e19_t08_fit" $D "$PU scripts/w04_bench.py"
add e22_w05_panel    "$(allof e22_w01_panel 10)" $D "$PU scripts/w05_panel.py"
add e22_w06_actbin   "$(allof e22_w01_actbin 7) e22_w04_bench" $D "$PU scripts/w06_actbin.py"
add e22_w07_central  "e22_w05_panel"    $D "$PU scripts/w07_central.py"
add e22_w08_pd       "$(allof e22_w01_panelpd 6) e22_w05_panel" $D "$PU scripts/w08_pd.py"
add e22_w09_weighted "e22_w08_pd"       $D "$PU scripts/w09_weighted.py"
add e22_w09b_total   "e22_w09_weighted" $D "$PU scripts/w09b_total.py"
add e22_w10_constructions "e22_w09_weighted" $D "$PU scripts/w10_constructions.py"
# ---- figures and the values drawn in them
add figures "e19_t01_extract e20_u03_pairs e21_v02_types e22_w06_actbin e22_w07_central e22_w09b_total e22_w10_constructions" figures "$PU make_figures.py"

trim() { echo "$1" | sed -E 's/^ +| +$//g'; }
if [ "${1:-}" = "list" ]; then for s in "${STEPS[@]}"; do IFS='|' read -r n r d c <<< "$s"; printf "%-26s needs: %s\n" "$(trim "$n")" "$(trim "$r")"; done; exit 0; fi
# The track table that the excluded script z01_lineages.py writes is an input (data_layout.md); put the copy of the repository in place when it is absent.
[ -f "$AGA_DATA_POSIX/exp3c/track_lineage.csv" ] || { mkdir -p "$AGA_DATA_POSIX/exp3c" && cp "$E/e02_lineage_map/inputs/track_lineage.csv" "$AGA_DATA_POSIX/exp3c/track_lineage.csv"; }
PREFIX="${1:-}"; ME="$$@$(hostname)"
# Repeated passes over the list: run the first step that is not done, not claimed and whose requirements are done.
# When nothing can run but another runner holds a lock, wait and pass again; when nothing can run and no lock exists, stop.
while :; do
  ran=0; pending=0; done0=$(ls "$ST"/*.done 2> /dev/null | wc -l)
  for s in "${STEPS[@]}"; do
    IFS='|' read -r n r d c <<< "$s"; n=$(trim "$n"); r=$(trim "$r"); d=$(trim "$d"); c=$(trim "$c")
    [ -n "$PREFIX" ] && [[ "$n" != "$PREFIX"* ]] && continue
    [ -f "$ST/$n.done" ] && continue
    pending=1
    [ -f "$ST/$n.lock" ] && continue
    ready=1; for q in $r; do [ -f "$ST/$q.done" ] || ready=0; done
    [ $ready -eq 1 ] || continue
    ( set -C; echo "$ME" > "$ST/$n.lock" ) 2> /dev/null || continue          # with noclobber the lock is created atomically: only one runner can take the step
    echo "$(date '+%F %T') $n: start ($c)"
    echo "=== $(date '+%F %T') $n: $c" >> "$LG/$n.log"            # append, so the log of a failed attempt survives a restart
    if (cd "$d" && eval "$c") >> "$LG/$n.log" 2>&1; then date '+%F %T' > "$ST/$n.done"; rm -f "$ST/$n.lock"; echo "$(date '+%F %T') $n: done"
    else rm -f "$ST/$n.lock"; echo "$(date '+%F %T') $n: FAILED (see $LG/$n.log)"; exit 1; fi
    ran=1; break
  done
  [ $pending -eq 0 ] && { echo "$(date '+%F %T') all steps done"; exit 0; }
  [ $ran -eq 1 ] && continue
  if ls "$ST"/*.lock > /dev/null 2>&1; then sleep 120
  elif [ "$(ls "$ST"/*.done 2> /dev/null | wc -l)" != "$done0" ]; then continue          # a step of another runner finished during this pass: scan again
  else echo "$(date '+%F %T') nothing runnable and no other runner active (a required step failed or is outside this prefix) -> stop"; exit 1; fi
done
