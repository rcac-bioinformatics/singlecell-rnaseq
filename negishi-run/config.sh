# negishi-run/config.sh: settings shared by every kit script. Sourced, not executed.
# shellcheck shell=bash disable=SC2034  # variables are used by the scripts that source this file
# Edit here, not in the scripts; any value can also be overridden from the environment.
# Values marked FILL-IN are placeholders: preflight FAILs while they still contain "FILL-IN".

# SLURM settings used by the episode 2 scripts (keep identical to the lesson; build_kit.py checks)
KIT_ACCOUNT=${KIT_ACCOUNT:-workshop}
KIT_PARTITION=${KIT_PARTITION:-cpu}
# QoS for kit jobs: empty means none is passed (the lesson's scripts set none)
KIT_QOS=${KIT_QOS:-}

# Staged workshop data on Depot. Read only for every kit script except publish_staged.sh and
# publish_rds.sh, which write here only after a dry-run listing and typed confirmation.
# The same path appears in learners/setup.md and episode 2; build_kit.py checks they agree.
STAGED=${STAGED:-/depot/workshop/data/scrna_workshop}
# Permission model of the staged copy: group (Depot group-readable) or other (world-readable)
STAGED_PERM=${STAGED_PERM:-group}

# R for episodes 3 to 8: image (primary) or module (fallback). run_R_chain.sh also accepts
# --r current, which uses whatever Rscript is on PATH (an OOD RStudio terminal, or an srun shell
# after module load).
KIT_R_MODE=${KIT_R_MODE:-image}
# image mode: the Open OnDemand RStudio (bioconductor) image, run with --cleanenv as the app does
OOD_SIF=${OOD_SIF:-/apps/biocontainers/images/r-bioconductor_3.23-R-4.6.1.sif}
# host package library the OOD app binds into the image (empty if the image is self-contained)
OOD_HOST_LIB=${OOD_HOST_LIB:-}
OOD_HOST_LIB_MOUNT=${OOD_HOST_LIB_MOUNT:-/opt/R/host-site-library}
# module mode: the module string that gives the same R as the OOD app (loaded after biocontainers)
R_MODULE=${R_MODULE:-FILL-IN/r-bioconductor-module}
# Account the OOD form uses (learners/setup.md says "provided by instructor"); R jobs run under it
OOD_ACCOUNT=${OOD_ACCOUNT:-rcac-rnaseq}
# OOD form resources (learners/setup.md): cores and wall time for the R jobs
OOD_CORES=${OOD_CORES:-16}
OOD_TIME=${OOD_TIME:-8:00:00}

# Where test runs live. Must be under /scratch/negishi/: the R episodes build their working
# directory as /scratch/negishi/$USER/scrna_workshop (see docs/ADAPTATIONS.md).
KIT_SCRATCH_BASE=${KIT_SCRATCH_BASE:-/scratch/negishi/$(id -un)/scrna-kit}

# Expected module defaults. The episodes load the default version (ml cellranger, ml star);
# episode 2 also hardcodes cellranger 10.0.0 in the whitelist extraction. Empty version = record only.
EXPECTED_MODULES="cellranger/10.0.0 star/"
# R and Bioconductor the OOD app should provide (empty = record only)
EXPECTED_R_VERSION=${EXPECTED_R_VERSION:-}
EXPECTED_BIOC_VERSION=${EXPECTED_BIOC_VERSION:-}

# Output of build_kit.py (do not edit by hand)
KIT_GEN_DIRNAME=generated
