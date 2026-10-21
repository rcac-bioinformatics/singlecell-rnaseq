# Negishi test kit (scRNA-seq lesson)

Runs every episode of the lesson on Negishi as a new learner would, in the real environment: episode 2 with `ml biocontainers` and the `rcac-rnaseq` account, episodes 3 to 8 in the Open OnDemand RStudio (bioconductor) R. It brings back one tarball with every output, figure, timing, warning, and anchor number the lesson needs. The structure, naming, and run loop follow the rnaseq-analysis kit (delivered 2026-10-06). A few safeguards come from the variant-calling kit: frozen anchors, `build_kit.py --check` before submitting, the bind-variable guard, the sacct summary, timestamped tarballs, checksummed fetches, and laptop self-tests.

`build_kit.py` generates every learner command from the episodes. R chunks are taken by label, shell blocks by a content anchor, and the build fails if an anchor stops matching. `docs/ADAPTATIONS.md` lists every place where running unattended differs from a learner. The staged data on Depot is only read. Each run works in a fresh copy under `$KIT_SCRATCH_BASE` (default `/scratch/negishi/$USER/scrna-kit`). Only `publish_staged.sh` and `publish_rds.sh` write to `/depot/workshop/data/scrna_workshop/`, and only after a listing and a typed confirmation.

## Before the first run: fill in config.sh

Preflight FAILs while any of these is still `FILL-IN`:

| Variable | What | Where to find it |
|---|---|---|
| `OOD_SIF` | the RStudio (bioconductor) image the OOD app runs (image mode, the default) | the OOD app's form or `submit.yml.erb`, or ask the OOD admins |
| `OOD_HOST_LIB` | host package library the app binds into the image (empty if none) | same |
| `R_MODULE` | module string giving the same R (module mode, `KIT_R_MODE=module`) | `module spider r-bioconductor` |
| `OOD_ACCOUNT` | account learners pick on the OOD form | learners/setup.md says "provided by instructor" |

`KIT_R_MODE=image` (default) runs each R episode with `apptainer exec --cleanenv` of `OOD_SIF`, with a fresh `HOME` and `R_LIBS_USER`, as rnaseq-analysis did. `KIT_R_MODE=module` loads `R_MODULE` instead.

## Command order

Run everything as `aseethar` from a Negishi login node, in `tmux` or `screen`.

| Step | Command | Where | Wall time |
|---|---|---|---|
| 0 | send the kit (below) | laptop | 1 min |
| 1 | `bash ~/singlecell-rnaseq/negishi-run/00_preflight.sh` | login node | 3 to 10 min |
| 2 | `bash ~/singlecell-rnaseq/negishi-run/01_learner_setup.sh 2026-10-20` | login node | 20 to 40 min (copies the staged tree, about 45 GB) |
| 3 | `bash ~/singlecell-rnaseq/negishi-run/submit_all.sh --run 2026-10-20 --dry-run` | login node | 1 min |
| 4 | `bash ~/singlecell-rnaseq/negishi-run/submit_all.sh --run 2026-10-20` | login node | 5 to 15 min (02-prep runs the episode 2 login-node commands), then jobs run for hours |
| 5 | `docs/ood_manual_check.md` (OOD form, and one episode from an RStudio terminal) | browser | 30 min |
| 6 | `bash ~/singlecell-rnaseq/negishi-run/collect.sh --run 2026-10-20` | login node | 1 to 3 min; rerun any time |
| 7 | bring the tarball back (below) | laptop | 1 min |

Step 1 must show no FAIL before you continue. WARN rows are carried into the summary.

What runs after submission (generated/steps.tsv):

| Steps | What | Expected |
|---|---|---|
| `02-cellranger`, `02-starindex` then `02-starsolo`, `02-compare` | the learner's own `run_cellranger.sh`, `build_star_index.sh`, `run_starsolo.sh` (64 cores each), then the episode's comparison blocks and the metrics | Cell Ranger about 7 h (the episode's callout), STAR index and STARsolo about 1.5 h together, plus queue waits |
| `r03` to `r06`, `r07`, `r08`, `08-sceasy` | continuous pass: episodes 3 to 8 as written, each a fresh R session; 03 reads the staged matrix the learner copied in setup (patch 01) | 1 to 15 min per episode (local pre-screen); episode 6 longer without presto |
| `c03` to `c06`, `c08` | cold pass: each episode in its own directory, seeded from the staged `filtered_feature_bc_matrix/` and `rds/` | in parallel, at once |
| `env-egress` | R and Bioconductor versions, `packageVersion()` of every package the code uses, then each network call (Monaco, HPCA, ifnb, enrichGO, enrichKEGG) with proxy and cache paths | 10 to 20 min |

Watch the queue with `squeue -u $USER -o '%.10i %.24j %.8T %.10M %R'`. Kit jobs are named `kit-<step>`; the learner scripts keep their own names (`cellranger`, `star_index`, `starsolo`).

### Running the R episodes without sbatch

From an OOD RStudio terminal, or from an `srun` shell after loading the R module, the same episodes run in that session's R:

```bash
bash ~/singlecell-rnaseq/negishi-run/run_R_chain.sh --run 2026-10-20 --pass cold --episodes 03,04 --r current
bash ~/singlecell-rnaseq/negishi-run/run_R_chain.sh --run 2026-10-20 --episodes env --r current
```

With `--pass continuous`, a failed episode skips the rest of its RDS chain (03 to 06, or 07 to 08), because the next one would have no input. `--r image` or `--r module` uses the configured environment instead of the current R.

## Sending the kit and returning the results

From the laptop. The kit reads the episodes, `learners/setup.md`, and `config.yaml`:

```bash
cd ~/svn/singlecell-rnaseq
python3 negishi-run/build_kit.py --check          # must say "matches a fresh build"
bash negishi-run/tests/run_all.sh                 # laptop self-tests, about 1 min
rsync -av --relative ./episodes ./learners ./config.yaml ./negishi-run aseethar@negishi.rcac.purdue.edu:singlecell-rnaseq/
```

`collect.sh` prints the command to bring the tarball back. Each collect writes a new timestamped tarball, so nothing is overwritten. Extract it into `negishi-run/records/`, which is gitignored. Do not extract it into the repo root: sandpaper would render its `SUMMARY.md`.

```bash
rsync -avP aseethar@negishi.rcac.purdue.edu:/scratch/negishi/aseethar/scrna-kit/scrna-kit-2026-10-20-20261021-0930.tar.gz ~/svn/singlecell-rnaseq/negishi-run/records/
tar -xzf ~/svn/singlecell-rnaseq/negishi-run/records/scrna-kit-2026-10-20-20261021-0930.tar.gz -C ~/svn/singlecell-rnaseq/negishi-run/records/
```

Start with `SUMMARY.md`. It covers:

- step status, with sacct wait, elapsed time and MaxRSS;
- the R environment and the per-call egress results;
- anchors (continuous and cold) against the episodes' output blocks;
- every R error and warning, verbatim;
- output blocks whose numbers the console no longer prints;
- the figures captured;
- runtime against the schedule in `learners/setup.md`;
- a per-episode summary table.

## Options and recovery

`submit_all.sh --run ID [--only all|shell|r|cold] [--from STEP] [--force] [--dry-run]`

- **Reruns are safe.** A step that already finished in this run is skipped: its marker says PASS or WARN, or, for a learner script, its outputs exist. `--force` reruns it. `submit_all.sh` refuses while jobs from an earlier submission of the same run are still queued.
- **A failed step cancels its dependents** (`--kill-on-invalid-dep=yes`). The logs are:
  - `runs/ID/logs/<step>.<jobid>.out`;
  - `logs/session-<step>.log` for shell steps;
  - `logs/R-<step>.log` for R steps, with per-chunk consoles in `records/console/<step>/`;
  - the learner scripts' own logs in `runs/ID/scratch/scrna_workshop/*_<jobid>.out`.

  Fix the cause, then resubmit from that step with `--from STEP`.
- **Before submitting anything,** `submit_all.sh` runs `build_kit.py --check`. It refuses if `SINGULARITY_BIND` or `APPTAINER_BIND` is set, because learners do not have them.
- **If the episodes change after the kit was generated,** preflight warns and `submit_all.sh` refuses. Rebuild on Negishi with `python3 negishi-run/build_kit.py` (Python 3.6 is enough) and start a new run ID.

## Staged data on Depot

`STAGED=/depot/workshop/data/scrna_workshop` (config.sh, learners/setup.md, episode 2). `staged_manifest.tsv` lists exactly what it must hold:

| Path | Source |
|---|---|
| `README.txt` | `staged-top/README.txt`: what each folder is for, which episode uses it, where it came from |
| `fastq/` | 6 FASTQs (the I1 files are optional) from 10x (`sources.tsv`) |
| `reference/refdata-gex-GRCh38-2024-A/` | from 10x (`sources.tsv`) |
| `filtered_feature_bc_matrix/` | the kit's own Cell Ranger run, or an earlier run of the episode 2 script whose counts equal the episode 3 output block |
| `rds/` | 5 files from the kit's continuous R chain |

The 3M whitelist is not staged; learners extract it in episode 2. Sizes and md5 stay `-` until frozen on Negishi.

| Task | Command |
|---|---|
| Check the staged copy (read only) | `bash check_staged.sh [DIR] [--quick] [--perm group\|other]` |
| Fetch what is missing from 10x into a staging tree on scratch (52 GB for the FASTQs; use `--submit`) | `bash stage_fetch.sh [--only fastqs\|refdata] [--submit] [--dry-run]` |
| Publish a staging tree | `bash publish_staged.sh --from $KIT_SCRATCH_BASE/staging/scrna_workshop --part fastq\|reference [--dry-run]` |
| Publish the filtered matrix from a run whose `02-cellranger` and `02-compare` passed | `bash publish_staged.sh --matrix-from-run RUN_ID [--dry-run]` |
| Publish an existing matrix (outside `$STAGED`) whose counts equal the episode 3 output block | `bash publish_staged.sh --matrix-from-path DIR [--dry-run]` |
| Publish the top-level README (refused while it has `FILL-IN`) | `bash publish_staged.sh --part readme [--dry-run]` |
| Publish the five RDS files from a run's continuous pass, with a side-by-side summary | `bash publish_rds.sh --run RUN_ID [--dry-run]` |
| Freeze the manifest from a verified copy; bring the output back for round 3 | `python3 make_staged_manifest.py /depot/workshop/data/scrna_workshop` |

### Building a new staged directory

Order used for the 2026-10-29 delivery, after the old directory was moved aside. Arun runs the Depot `mv` and `mkdir`; everything after that goes through the publish scripts.

1. `bash stage_fetch.sh --submit` (FASTQs and refdata from 10x, verified, into scratch).
2. `bash publish_staged.sh --matrix-from-path <old copy>/.../filtered_feature_bc_matrix` (does not wait for step 1).
3. When the fetch job ends: `publish_staged.sh --from $KIT_SCRATCH_BASE/staging/scrna_workshop --part fastq`, then `--part reference`.
4. Kit run: `01_learner_setup.sh`, `submit_all.sh --only shell` and `--only r`.
5. `publish_rds.sh --run RUN_ID` once r03 to r08 pass, then `submit_all.sh --only cold`.
6. Fill the `FILL-IN` lines of `staged-top/README.txt`, rsync the kit, `publish_staged.sh --part readme`.
7. `check_staged.sh` (full, with md5), then `make_staged_manifest.py` to freeze the manifest.

How the fetch and publish scripts behave:

- **`stage_fetch.sh`** verifies every archive against the byte count and S3 multipart ETag recorded in `sources.tsv`. 10x publishes no md5; `py/s3etag.py` reproduces their ETag. It records the md5 of every extracted file in `staging/provenance.tsv`.
- **Both publish scripts:**
  - list every file as add, replace or same, and print the exact `cp` and `chmod` commands;
  - then ask you to type the destination directory; anything else cancels;
  - back up replaced files to `$KIT_SCRATCH_BASE/publish-backups/<time>/` and never delete;
  - log to `$KIT_SCRATCH_BASE/publish.log`, then run `check_staged.sh --quick`.
- **`publish_rds.sh`** reads the handoff copies the run kept right after each save chunk, so episode 8 overwriting `ifnb_annotated.rds` cannot reach Depot.

Learners need read access through the `workshop` Unix group; preflight checks the permissions.

## What is not stored in git

`generated/` (rebuild with `python3 negishi-run/build_kit.py`), `staged-data/`, `records/` (static reports and returned tarballs), and `*.tar.gz`. `updates-patches/` at the repo root is tracked.

## What could not be verified before the first Negishi run

On the laptop:

- `build_kit.py` builds deterministically (`--check`).
- shellcheck 0.10.0 passes on every kit script and generated job.
- `tests/run_all.sh` passes: mutations `--check` must catch, the bind-variable guard, and the summarizer on synthetic records.
- `py/s3etag.py` reproduces a real 10x ETag.
- The R chain ran end to end in Docker (`arnstrm2/r-bioconductor:3.23-R-4.6.1`, Seurat 5.5.1) on the 10x public matrix. That is a Cell Ranger 3.0 matrix, so its numbers differ from the lesson; errors are findings, warnings are provisional. See maintenance/log/updates.md, Runtime record.

What depends on Negishi, and what covers it:

| Not verifiable off-cluster | What covers it |
|---|---|
| Module defaults (cellranger 10.0.0, star), `BIOC_IMAGE_DIR` and the 10.0.0 SIF used for the whitelist | preflight |
| `rcac-rnaseq` account membership, partition limits, memory per core | preflight: sacctmgr, scontrol, sinfo, `sbatch --test-only` on every script |
| Staged data complete and group-readable | preflight runs `check_staged.sh --quick`; full check with `check_staged.sh` |
| OOD image (or module), R and Bioconductor versions, every package the code uses | preflight, one row per package from `generated/r_packages.tsv` (the episode code plus `package_requirements.tsv`: presto blocks episode 6, metap and multtest block episode 7, sceasy 0.0.7 or older is a WARN for Seurat v5), and step env-egress |
| Compute-node egress per call (ExperimentHub, SeuratData repository, GitHub, KEGG, Bioconductor, CRAN) | preflight (curl) and env-egress (each episode call in R, with proxy and cache paths) |
| The OOD form, launch time, the interactive session | `docs/ood_manual_check.md` |
| Real runtimes, memory, queue waits | the run itself (sacct and VmHWM per chunk in SUMMARY.md) |
