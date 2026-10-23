# Adaptations: where the kit does not run episode code exactly as a learner would

The rule is learner fidelity. Episode commands run verbatim, with the same module loads, paths, flags and SLURM headers, in the real Negishi environment. `build_kit.py` extracts every block the kit runs from the episodes; nothing is hand-copied. This file lists every place where running unattended forces a difference, and why. If an episode command fails as written, the run records the failure; the kit never works around it.

## Where the learner directory lives

| What | Learner | Kit | Why |
|---|---|---|---|
| `$SCRATCH`, `$RCAC_SCRATCH` (episode 2, setup.md) | `/scratch/negishi/<user>` | `$KIT_SCRATCH_BASE/runs/<RUN_ID>/scratch` | For `aseethar`, the real scratch could hold staged or earlier copies. Pointing both variables at a fresh test root makes every `${RCAC_SCRATCH}/scrna_workshop` path land in a fresh copy. `lib.sh` refuses to run if a learner directory resolves inside `$STAGED`. |
| `$USER` in the R episodes' `set-paths` | user name | `aseethar/scrna-kit/runs/<RUN_ID>/scratch` (continuous) or `.../cold/<ep>/scratch` (cold) | Episodes 3 to 8 build `work_dir` from `"/scratch/negishi/", Sys.getenv("USER")`, not from `RCAC_SCRATCH`. Setting `USER` (inside the R session only) to the test root's path below `/scratch/negishi/` makes the unchanged code resolve to the test copy. For a real learner both forms give the same path; preflight checks `$RCAC_SCRATCH = /scratch/negishi/$USER`. |

## Episode 2 (shell, SLURM)

| Episode block | Learner | Kit |
|---|---|---|
| `mkdir`, `rsync`, whitelist extraction, `ls` and `zcat \| head` checks, `cat << 'EOF' > ...` for the three job scripts | login node, typed in order | Step `02-prep`: the same blocks in order in a learner-like shell on the login node, run by `submit_all.sh`. The episode's `rsync` repeats the setup.md copy, so it is a near no-op, as for a learner. The three written scripts are compared with the episode text (`kit_check_written`). |
| `sbatch --test-only` for each of the three scripts (patch 20) | login node, in class | Step `02-prep`, as written, after the scripts are written. Submits nothing. |
| `sbatch run_cellranger.sh`, `sbatch .../build_star_index.sh`, `sbatch .../run_starsolo.sh` | after the workshop (patch 20; the instructor runs them beforehand); STARsolo after the learner sees the index finish | The written scripts are submitted from `scrna_workshop/` with only submission options added: `--parsable` (to read the job ID), `--dependency=afterok` (STARsolo after the index), `--kill-on-invalid-dep=yes`, and `--qos` only if `KIT_QOS` is set. Script content is never changed. The `squeue` line in the submit block is not run. |
| Cell Ranger array example (callout, placeholder samples) | template | `sbatch --test-only` only (header validation), in preflight and `submit_all.sh --dry-run`. Not run. |
| Challenge 1 solution fragment (v2 flags) | reading | Not run; it is a fragment. |
| Output listing and comparison blocks (`ls .../outs/`, `ls .../Solo.out/Gene/filtered/`, cell and gene counts) | after both jobs finish | Step `02-compare`, a 1-core batch job after both jobs (`afterany`, so a STARsolo failure still leaves the Cell Ranger half). Kit extras copy `metrics_summary.csv`, `web_summary.html`, STARsolo `Summary.csv`, `Log.final.out` and `Barcodes.stats`. They count shared barcodes, and compare the fresh matrix (md5 of file and content) with the staged `filtered_feature_bc_matrix/`. |

Learner-like shell: the generated sessions run in a child `bash` without `set -e`, `set -u`, `pipefail` or `set -E`, as in a terminal. An `ERR` trap logs every failing top-level command as `### KIT-ERR`, and the step fails if there is one. `set -E` is left out on purpose: the variant-calling kit saw it fire inside Lmod's shell functions (SIGPIPE, rc 141). Lmod's `module` function is used as exported by `sbatch`. If it is not exported, the session sources Lmod's init file, which gives the same function a login shell has.

`$TMPDIR` is namespaced per job (`$RUN/tmp/<jobid>`) in kit jobs only, because Negishi's `$TMPDIR` is a node-wide `/tmp`. The learner scripts run with whatever SLURM gives them.

## Episodes 3 to 8 (R)

Each episode is one fresh R session in the configured R:

- **image (default):** `apptainer exec --cleanenv` of the OOD image, with a new empty `HOME` and `R_LIBS_USER` (a new learner), plus the host library bind if `OOD_HOST_LIB` is set.
- **module:** `module load $R_MODULE`.
- **current:** the R on `PATH` in an OOD RStudio terminal or an srun shell, used by `run_R_chain.sh --r current`.

Resources match the OOD form in setup.md: `OOD_CORES` and `OOD_TIME` on `OOD_ACCOUNT`.

A learner may keep one RStudio session across episodes; a fresh session per episode is the stricter case. It is what a learner gets after a disconnect, or when they start an episode cold. The episode 6 `%>%` finding comes from this.

`R/run_blocks.R` runs the chunks in order with `source(echo = TRUE, print.eval = TRUE)`, which prints commands and values as the console does. Differences from pasting into RStudio:

- an error stops the rest of that chunk (pasting would continue line by line), and the run continues with the next chunk so later failures are visible too;
- warnings and messages are printed where they occur and logged with the chunk label;
- each chunk's plots go to PNG files (`<ep>_<label>_NN.png`, 1400 x 1000 px, 150 dpi) instead of the Plots pane;
- `VmHWM` and `VmRSS` are recorded after every chunk.

Challenge, solution and callout chunks run in order. A failure in one of them is reported, and makes the step WARN rather than FAIL.

| Chunk | Kit | Why |
|---|---|---|
| episode 3 `challenge1` | not run | Fill-in-the-blank template with literal `???`; `solution1` runs the filled-in code. |
| plain ```` ```r ```` fences (Seurat v5 vs v4 syntax in episode 3, SCTransform in episode 4, the sub-clustering list item in episode 6, the sessionInfo tip in episode 8) | not run | Display-only examples, not knitr chunks. |
| episode 4 `install-glmgampoi` | run as written | Learners run it; it needs Bioconductor egress and installs into the fresh `R_LIBS_USER`. |

### Continuous and cold passes

- **Continuous** (`r03` to `r08`): one learner copy, the RDS chain as written. Episode 3 reads the matrix its `set-paths` names; `build_kit.py` writes that path to `generated/ep3_matrix_path.txt`. With patch 01 applied (2026-10-16), that path is the staged copy in `scrna_workshop/filtered_feature_bc_matrix/`, so `r03` starts at once. If the episode pointed at `cellranger_output/...` again, the build would make `r03` wait for `02-cellranger`. `r07` starts from SeuratData and has no dependency.
- **Cold** (`c03` to `c06`, `c08`): each episode in its own new directory, so the parallel jobs never share files. Each is seeded with what a learner who skipped the previous episode needs, copied from `$STAGED`:
  - `c03`: the staged `filtered_feature_bc_matrix/`, copied to the path in `generated/ep3_matrix_path.txt`, where episode 3 looks;
  - `c04`, `c05`, `c06`, `c08`: every staged `rds/<name>.rds` the episode reads and does not save itself, copied into `scrna_workshop/` where `readRDS(paste0(work_dir, ...))` looks. `build_kit.py` derives the list from the episode's `readRDS()` and `saveRDS()` calls and writes it to `generated/cold_seed.tsv` (`c08` gets `ifnb_annotated.rds` and `pbmc_annotated.rds`); the build stops if a file is not listed under `rds/` in `staged_manifest.tsv`.

  The episodes' cold-start callouts (patch 05, applied 2026-10-22) tell learners to make the same copy. There is no `c07`: episode 7 starts from SeuratData in both passes.

### Kit-only additions (never change learner objects, never print to the learner console)

Rows marked `extra` in `generated/rplans/*.tsv` run after a named chunk, in a child environment of the learner's global environment:

- `kit_anchor()` writes the anchor values in `anchors.tsv` to `records/anchors/<step>.tsv`;
- `kit_handoff()` copies each saved handoff RDS to `$RUN/handoff/<step>/`, the source for `publish_rds.sh`;
- `kit_marker_check()` records which manual label has the highest mean expression of each canonical marker, and the cluster-by-label table;
- `kit_table()` records the SingleR vs manual cross-table;
- `kit_top_genes()` records the top 20 DE genes;

The step `env-egress` is not episode code. It runs `R/env_egress.R` in the same R environment, with a fresh `HOME` and user library. It records:

- R and Bioconductor versions, and `packageVersion()` of every package the episode code uses (derived by `build_kit.py` from `library()`, `::`, `install()` and `test.use` in the extracted chunks);
- proxy variables and cache locations;
- each network call the episodes make, timed and reported pass or fail: `celldex::MonacoImmuneData()`, `HumanPrimaryCellAtlasData()`, `SeuratData::InstallData("ifnb")`, `clusterProfiler::enrichGO()`, `enrichKEGG()`, and BiocManager.
