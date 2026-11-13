# Tissue-group and regulatory-layer partitioning of complex-trait heritability with AlphaGenome Atlas predictions

Analysis code and result files.

AlphaGenome Atlas allelic-effect predictions for 5.83 million common SNPs are grouped into 12 tissue groups and six
regulatory layers (CAGE, DNase, H3K27ac, H3K4me1, H3K4me3, RNA-seq), turned into annotations (the top 5% of SNPs by
predicted effect) and used to partition the heritability of complex traits with stratified LD score regression (S-LDSC)
conditional on baselineLD v2.2. The analyses cover 14 traits, six of them with an earlier and a later GWAS, and a public
panel of 106 GWAS. The repository holds the scripts in the order in which they were run, the result files they wrote, and
the scripts that draw the figures from those result files. It is meant for reproducing the analyses.

## Contents

| Path | What it is |
|---|---|
| `experiments/e01_…` to `e22_…` | One folder per analysis stage: `scripts/`, `results/` (the result files written by the scripts) and, in three stages, `inputs/` (small tables that are inputs rather than downloads: locus lists, track–lineage table) |
| `run_all.sh` | Runs every step in order with dependencies, done-markers and logs; resumable; several copies can share the work |
| `config.py` | The single place where the data root is set (`AGA_DATA`) |
| `data_layout.md` | What must be downloaded, where it goes, and what the pipeline writes |
| `CODE_INDEX.md` | Which scripts stand behind each section of the Results, what each script does, and which result files each figure and table is made from |
| `figures/` | The figure scripts (`fig1.py` … `fig5.py`, `figS1.py` … `figS5.py`, run together by `make_figures.py`), which read only the result files; `figures/out/` holds the figures (PNG) and the values drawn in them (`numbers_<figure>.json`) |
| `verify/compare_results.py` | Compares a reproduction run with the result files in this repository; `verify/report.md` and `verify/report_e19_e22.md` are the reports of the runs described under *Reproduction* |
| `tools/` | `ldsc-python3_local_patch.diff` (local changes to LDSC), `make_run_copy.py` (makes the copy to run in) |
| `requirements.txt` | Package versions |

Stages e01–e18 build the annotations and LD scores and analyse 11 traits and three further traits with the GWAS first
used. Stages e19–e22 add the later GWAS of six traits to a benchmark of 14 traits (the later GWAS where there are two)
and analyse the panel of 106 GWAS; they reuse the annotations, LD scores and run files of the earlier stages. The scripts
and file names call the tissue groups lineages.

## Running

1. Install Python 3.13 with the packages in `requirements.txt`, LDSC (two forks, see `data_layout.md`) and plink2.
   The pipeline was run on Windows only, from Git Bash; one script addresses the plink2 binary as `tools/plink2.exe`.
2. Lay out the data as described in `data_layout.md` under one folder and point `AGA_DATA` at it. The AlphaGenome Atlas
   scores are the only input that cannot be downloaded from a public archive; they are queried from the Atlas API
   (the query scripts are not part of this repository, see below).
3. Make a copy of the repository to run in, so that the result files in the repository stay untouched:

       python tools/make_run_copy.py            # -> repro_run/ with empty results/ folders

4. Run the pipeline from that copy. In Git Bash on Windows:

       cd repro_run
       AGA_DATA=D:/data AGA_DATA_POSIX=/d/data bash run_all.sh        # everything, in order
       bash run_all.sh list                                           # the steps and their dependencies
       bash run_all.sh e05_                                           # one stage

   The pipeline has 220 steps (analysis steps, the self-checks and dry runs that preceded each stage, the pieces into
   which the fits of e20–e22 are split, and the figure step). Each step writes a done-marker under `$AGA_DATA/_status/`
   and a log under `$AGA_DATA/_logs/`, so an interrupted run continues where it stopped, and the S-LDSC scripts skip fits
   whose output exists. A runner that is killed leaves the lock file of its step behind; delete
   `$AGA_DATA/_status/*.lock` before starting again. A second copy of `run_all.sh` started in another terminal takes the steps the first has not
   claimed: the single-layer-map fits of `e03` take about two days and can run beside the rest, and the fits of e20–e22
   are split into 5 to 12 pieces each so that as many copies can run them side by side. The steps of e01–e18 took about
   118 hours of computation on one workstation (24 cores, 96 GB; most steps use a single core) and those of e19–e22 about
   88 hours summed over the pieces (about ten hours with ten copies of the runner). On that Windows workstation LDSC
   very occasionally (3 of about 2,000 fits in e01–e18) failed while reading intact LD score files; `y04_sldsc.py`
   retries such a fit, and for any other script simply starting `run_all.sh` again repeats the step, which skips the
   fits already done.

5. Check the result: `python verify/compare_results.py repro_run` compares every file under `experiments/*/results/`
   and `figures/out/` of the run with the files of this repository and writes a report (`--out` sets its path;
   `--data-root` names the data root so that it is masked in the logs).

The figure scripts set Arial (regular, bold and italic) and were run with matplotlib 3.10.8. They check the positions of
the text they draw and stop when labels overlap or leave the canvas, so with another font or another matplotlib version
they may stop at such a check, or draw with a substitute font where a script does not check; matplotlib then prints a
warning that the font was not found. With the same font and version they reproduce the PNG files in `figures/out/` byte
for byte.

All random choices (track resampling, bins, block alignment) use fixed seeds, and LD scores and S-LDSC use no random
numbers, so a run with the same inputs and package versions reproduces the result files: many of them byte for byte, and
the others in content, with numbers that agree within a relative tolerance of 1e-6 (S-LDSC fits that are computed again
can differ from the stored ones in the last digits; see *Reproduction*).

## Reproduction

**Stages e01–e18.** The pipeline was run once more in full, from a copy made as above and against the downloaded inputs
only, on 2026-10-30 to 2026-11-02 (one power cut in between; the run resumed from its markers). `verify/report.md` is the
comparison with the result files in the repository:

| Files | Outcome |
|---|---|
| 133 result tables and text outputs | identical (67 byte-identical; 66 identical in content, numbers within rtol 1e-6) |
| 33 run logs | identical except for timestamps, elapsed times and the data-root path |
| 2 run logs | differ for a documented reason without any value differing (see the report) |
| 1 output of the Atlas query script | not a pipeline output |

**Stages e19–e22 and the figures.** All steps of e19–e22 and the figure step were run on 2026-11-13, from a copy made as
above and with `run_all.sh`, in the data root written by the full run of e01–e18 (inputs and the outputs of e01–e18
only), so that every annotation, LD score and S-LDSC fit of these stages was computed anew. `verify/report_e19_e22.md`
is the comparison with the files in the repository:

| Files | Outcome |
|---|---|
| 32 result tables | identical (13 byte-identical; 19 identical in content, numbers within rtol 1e-6) |
| 23 text outputs | identical line by line (16) or except for elapsed times and the data-root path (7) |
| 86 logs of the pieces of the split fitting and LD-score steps | identical line by line (6) or except for timestamps, elapsed times and the data-root path (80) |
| 10 figures (PNG) | byte-identical |
| 11 files of the values drawn in the figures | identical (6 byte-identical; 5 identical in content, numbers within rtol 1e-6) |

The S-LDSC fits computed anew differ from the stored ones in the last digits (relative differences below 1e-9 in the
tables computed from them and below 1e-6 in the tables of the mixed model), so the tables are identical in content and
not always byte for byte.

## Data

No data are stored in this repository. Everything the pipeline reads is listed in `data_layout.md` with its source;
in short:

| Input | Where to get it |
|---|---|
| AlphaGenome Atlas allelic-effect and activity scores | AlphaGenome Atlas API (Google DeepMind); the query scripts are not included, see below |
| baselineLD v2.2 annotations, LD scores, 1000G EUR plink files and HapMap3 weights (hg38) | [Zenodo 10515792](https://doi.org/10.5281/zenodo.10515792) |
| 1000 Genomes high-coverage GRCh38 genotypes | plink2 resources |
| GWAS summary statistics of the 14 traits (earlier and later GWAS) | GWAS Catalog, Price lab processed GWAS, consortium downloads (`n00_download.sh`, `r00_download.sh`); listed in `data_layout.md` |
| Panel of 106 GWAS in LDSC format | [Zenodo 10515792](https://doi.org/10.5281/zenodo.10515792), `sumstats_indep107.tgz` |
| FANTOM5 CAGE, ENCODE DNase and H3K27ac peaks | downloaded by `e04/scripts/w01_measured_sources.py` |
| FunGen-xQTL fine-mapped eQTLs | NIAGADS NG00184 (`s00_download.sh`) |
| GENCODE v46 | GENCODE |

The small tables that are inputs rather than downloads (locus lists, the track–lineage table) are in the `inputs/`
folders of e01, e02 and e12. The result files under `experiments/*/results/` are the tables and text summaries from which
the figures and tables are made (about 4 MB in all); the S-LDSC run files, LD scores and annotations they were computed from
are not included and are rebuilt by the pipeline.

## Not included

- **AlphaGenome Atlas API access code.** Four scripts that query the Atlas are left out: `x02_atlas_query.py`,
  `z01_lineages.py` (track metadata; its lineage rules are in `e02_lineage_map/scripts/lineage_rules.py` verbatim),
  `z02_atlas_query.py` and `l01_atlas_active.py`. The scorers, windows and client version are described in the Methods
  of the manuscript.
- **Data.** Atlas predictions and annotations derived from them, LD scores, S-LDSC run files and GWAS summary
  statistics are not included (`.gitignore`); see *Data* above. Redistribution of Atlas-derived annotations is
  subject to the AlphaGenome terms of use.
- **LDSC itself.** LD scores and munging used [abdenlab/ldsc-python3](https://github.com/abdenlab/ldsc-python3) at
  commit b6176a4 with the local changes in `tools/ldsc-python3_local_patch.diff`; regression used the
  [brainvisa/ldscore](https://github.com/brainvisa/ldscore) fork at commit a09a869, unmodified. Both are GPL-3.0.
- **The manuscript.** Its text and the scripts that lay out its tables are not part of this repository; `CODE_INDEX.md`
  lists the result files behind every figure and table.

The docstrings of some scripts name the dated analysis protocol of their stage (`PROTOCOL_*`); the protocol documents
are not part of this repository, and the design of each stage is described in the Methods of the manuscript.

Code comments and the text outputs of the scripts are in Korean; this documentation is in English.

## Versions

The first public version of this repository (commit `e82ca8b`) holds stages e01–e18 and the figure script that went with
them (`paper/`). The current version adds stages e19–e22 and the scripts in `figures/`, which replace that figure script.
The scripts and result files of e01–e18 are unchanged, apart from one line of the header comment or docstring in nine
scripts (six `run_pipeline.sh` files and the three summary scripts of e01–e03).

## License

MIT (see `LICENSE`). LDSC, which this code calls but does not include, is GPL-3.0.
