#!/usr/bin/env python3
"""build_kit.py: generate the Negishi test kit from the lesson episodes.

Run from anywhere:  python3 negishi-run/build_kit.py [--check]
Writes negishi-run/generated/ (replaced on every run). --check builds in memory and fails if
generated/ differs from a fresh build (run it before sending the kit).

Code is never hand-copied: every learner command the kit runs is extracted from the episodes
(R chunks by label, shell blocks by a content anchor), and the build fails loudly when an
anchor matches zero or several blocks, so an episode that changes shape cannot be silently
mis-tested. Output is deterministic (no timestamps, sorted listings), so rerunning on unchanged
episodes gives byte-identical files. Python 3.6 compatible (the Negishi system python).
"""
import hashlib
import pathlib
import re
import shutil
import sys

KIT = pathlib.Path(__file__).resolve().parent
ROOT = KIT.parent
GEN = KIT / "generated"
EPISODES = ROOT / "episodes"
SETUP = ROOT / "learners" / "setup.md"

sys.dont_write_bytecode = True  # keep the kit directory free of __pycache__
sys.path.insert(0, str(KIT / "py"))
from lessonparse import parse_rmd  # noqa: E402

ORDER = ["introduction", "raw-data-processing", "quality-control", "normalization-and-features",
         "dimensionality-reduction-clustering", "cell-type-annotation", "multi-sample-integration",
         "differential-expression"]
cfg = (ROOT / "config.yaml").read_text()
m_eps = re.search(r'^episodes:\s*\n((?:-\s+\S+\s*\n)+)', cfg, re.M)
listed = re.findall(r'^-\s+(\S+)\.Rmd\s*$', m_eps.group(1), re.M) if m_eps else []
if not listed:
    sys.exit("build_kit: could not read the episodes: list from config.yaml")
if listed and listed != ORDER:
    sys.exit("build_kit: config.yaml episode order changed: %s; update ORDER" % listed)
KEY = {name: "%02d" % (i + 1) for i, name in enumerate(ORDER)}
EPS = {KEY[name]: parse_rmd(EPISODES / (name + ".Rmd")) for name in ORDER}
SET = parse_rmd(SETUP)
STAGED_PATH = "/depot/workshop/data/scrna_workshop"


def die(msg):
    sys.exit("build_kit: " + msg)


def find(ep, anchor, lang="bash"):
    """The single block of episode ep (or 'setup') whose text contains anchor."""
    e = SET if ep == "setup" else EPS[ep]
    hits = [b for b in e.chunks if anchor in b.text and (lang is None or b.lang == lang)]
    if len(hits) != 1:
        die("anchor %r in %s matched %d %s blocks (need exactly 1)" % (anchor, ep, len(hits), lang))
    return hits[0]


def chunk(ep, label):
    hits = [c for c in EPS[ep].r_chunks() if c.label == label]
    if len(hits) != 1:
        die("R chunk %r in episode %s found %d times" % (label, ep, len(hits)))
    return hits[0]


def rel(ep):
    return "episodes/%s.Rmd" % ORDER[int(ep) - 1]


# ---------------------------------------------------------------- episode 2 shell blocks
HEREDOC = re.compile(r"^cat << 'EOF' > (\$\{RCAC_SCRATCH\}/scrna_workshop/([\w.]+))\n(.*)\nEOF\n$", re.S)


def heredoc(anchor):
    b = find("02", anchor)
    m = HEREDOC.match(b.text)
    if not m:
        die("episode 2 block at line %d is no longer a single cat << 'EOF' > ... heredoc" % b.start)
    return b, m.group(2), m.group(3) + "\n"


B2 = dict(
    mkdir=find("02", "mkdir -p ${RCAC_SCRATCH}/scrna_workshop/{fastq,reference,cellranger_output,starsolo_output}"),
    rsync=find("02", "rsync -avP %s/ ${RCAC_SCRATCH}/scrna_workshop/" % STAGED_PATH),
    whitelist=find("02", "3M-february-2018_TRU.txt.gz"),
    ls_fastq=find("02", "ls -1 ${RCAC_SCRATCH}/scrna_workshop/fastq/"),
    ls_ref=find("02", "ls ${RCAC_SCRATCH}/scrna_workshop/reference/\n"),
    head_r1=find("02", "zcat pbmc_10k_v3_S1_L001_R1_001.fastq.gz | head -8"),
    head_r2=find("02", "zcat pbmc_10k_v3_S1_L001_R2_001.fastq.gz | head -8"),
    ls_refdata=find("02", "ls ${RCAC_SCRATCH}/scrna_workshop/reference/refdata-gex-GRCh38-2024-A/"),
    sub_cr=find("02", "sbatch run_cellranger.sh"),
    ls_outs=find("02", "ls ${RCAC_SCRATCH}/scrna_workshop/cellranger_output/pbmc10k/outs/"),
    sub_idx=find("02", "sbatch ${RCAC_SCRATCH}/scrna_workshop/build_star_index.sh"),
    sub_solo=find("02", "sbatch ${RCAC_SCRATCH}/scrna_workshop/run_starsolo.sh"),
    ls_solo=find("02", "ls ${RCAC_SCRATCH}/scrna_workshop/starsolo_output/starsolo_pbmc10k_Solo.out/Gene/filtered/"),
    cmp_cr=find("02", 'echo "Cell Ranger cells:"'),
    cmp_solo=find("02", 'echo "STARsolo cells:"'),
)
H_CR = heredoc("cat << 'EOF' > ${RCAC_SCRATCH}/scrna_workshop/run_cellranger.sh")
H_IDX = heredoc("cat << 'EOF' > ${RCAC_SCRATCH}/scrna_workshop/build_star_index.sh")
H_SOLO = heredoc("cat << 'EOF' > ${RCAC_SCRATCH}/scrna_workshop/run_starsolo.sh")
B_ARRAY = find("02", "#SBATCH --array=")
SCRIPTS = {"run_cellranger.sh": H_CR, "build_star_index.sh": H_IDX, "run_starsolo.sh": H_SOLO}
for name, (b, fname, body) in SCRIPTS.items():
    if fname != name:
        die("episode 2 heredoc at line %d writes %s, expected %s" % (b.start, fname, name))
if not B_ARRAY.text.startswith("#!/bin/bash"):
    die("episode 2 array example (line %d) is not a SLURM script" % B_ARRAY.start)

# the staged path learners copy from must match config.sh STAGED default and setup.md
bS_copy = find("setup", "rsync -avP %s/ ${RCAC_SCRATCH}/scrna_workshop/" % STAGED_PATH)
bS_verify = find("setup", "ls ${RCAC_SCRATCH}/scrna_workshop/filtered_feature_bc_matrix/")
if "STAGED=${STAGED:-%s}" % STAGED_PATH not in (KIT / "config.sh").read_text():
    die("config.sh STAGED default differs from the path in learners/setup.md and episode 2 (%s)" % STAGED_PATH)


def header_of(text):
    h = {}
    for line in text.splitlines():
        m = re.match(r'#SBATCH\s+--([\w-]+)(?:=(\S+))?', line.strip())
        if m:
            h[m.group(1)] = m.group(2) or ""
    return h


# the account the episode 2 #SBATCH headers must name; the kit's own jobs use config.sh KIT_ACCOUNT
LESSON_ACCOUNT = "rcac-rnaseq"   # patch 09, 2026-10-21
_cfg_sh = (KIT / "config.sh").read_text()
KIT_ACCOUNT = re.search(r'^KIT_ACCOUNT=\$\{KIT_ACCOUNT:-([^}]+)\}', _cfg_sh, re.M).group(1)
OOD_ACCOUNT = re.search(r'^OOD_ACCOUNT=\$\{OOD_ACCOUNT:-([^}]+)\}', _cfg_sh, re.M).group(1)
HEADER_ROWS = ["script\tline\taccount\tpartition\tnodes\tntasks\ttime\tmem\tverdict"]
for name, (b, fname, body) in sorted(SCRIPTS.items()) + [("cellranger_array_template.sh", (B_ARRAY, "", B_ARRAY.text))]:
    h = header_of(body)
    problems = []
    for k, v in (("account", LESSON_ACCOUNT), ("partition", "cpu"), ("nodes", "1")):
        if h.get(k) != v:
            problems.append("--%s=%s (expected %s)" % (k, h.get(k, "MISSING"), v))
    if "ntasks" not in h:
        problems.append("no --ntasks")
    if "mem" in h or "mem-per-cpu" in h:
        problems.append("sets --mem")
    HEADER_ROWS.append("\t".join([name, str(b.start), h.get("account", "-"), h.get("partition", "-"), h.get("nodes", "-"),
                                  h.get("ntasks", "-"), h.get("time", "-"), h.get("mem", "-"), "; ".join(problems) or "ok"]))
    # the three scripts the kit submits must be valid as written; the array is a template (updates.md)
    if problems and name in SCRIPTS:
        die("%s (%s:%d): %s" % (name, rel("02"), b.start, "; ".join(problems)))
# the tool modules the episode 2 scripts load must be the ones config.sh EXPECTED_MODULES names
# (preflight and collect.sh load exactly those strings)
_ep_mods = set()
for name, (b, fname, body) in SCRIPTS.items():
    for line in body.splitlines():
        if re.match(r'\s*ml\s', line) and "purge" not in line:
            _ep_mods.update(t for t in line.split()[1:] if t != "biocontainers")
_cfg_mods = set(re.search(r'^EXPECTED_MODULES="([^"]*)"', _cfg_sh, re.M).group(1).split())
if _ep_mods != _cfg_mods:
    die("episode 2 loads modules %s but config.sh EXPECTED_MODULES is %s" % (" ".join(sorted(_ep_mods)), " ".join(sorted(_cfg_mods))))
H_TIME = {name: header_of(body)["time"] for name, (b, f, body) in SCRIPTS.items()}
H_NTASKS = {name: int(header_of(body)["ntasks"]) for name, (b, f, body) in SCRIPTS.items()}

# ---------------------------------------------------------------- sessions (shell)
S = {}
S["02-prep"] = [("block", "02", B2[k]) for k in ("mkdir", "rsync", "whitelist", "ls_fastq", "ls_ref", "head_r1", "head_r2",
                                                 "ls_refdata")]
S["02-prep"] += [("block", "02", SCRIPTS[n][0]) for n in ("run_cellranger.sh", "build_star_index.sh", "run_starsolo.sh")]
S["02-prep"] += [("extra", "the written job scripts are exactly the episode text",
                  "\n".join('kit_check_written "$RCAC_SCRATCH/scrna_workshop/%s" "$KIT_GEN/learner-scripts/%s"' % (n, n)
                            for n in sorted(SCRIPTS))),
                 ("extra", "listings of the learner copy",
                  'ls -la "$RCAC_SCRATCH/scrna_workshop" > "$RUN/records/listings/02-scrna_workshop.txt"; '
                  'ls -la "$RCAC_SCRATCH/scrna_workshop/fastq" "$RCAC_SCRATCH/scrna_workshop/reference" > "$RUN/records/listings/02-fastq-reference.txt"; '
                  'wc -l < "$RCAC_SCRATCH/scrna_workshop/reference/3M-february-2018.txt" > "$RUN/records/listings/02-whitelist-lines.txt"')]
S["02-compare"] = [("block", "02", B2[k]) for k in ("ls_outs", "ls_solo", "cmp_cr", "cmp_solo")]
S["02-compare"] += [("extra", "Cell Ranger and STARsolo metrics, barcode agreement, fresh vs staged matrix", "kit_ep2_metrics")]

# ---------------------------------------------------------------- R plans
R_EPS = ["03", "04", "05", "06", "07", "08"]
OPTIONAL_CTX = ("challenge", "solution", "spoiler", "callout", "discussion")
R_SKIP = {
    ("03", "challenge1"): "fill-in-the-blank template (contains ???); solution1 runs the filled-in code",
}
for (ep, lab), why in R_SKIP.items():
    chunk(ep, lab)  # must exist
if "???" not in chunk("03", "challenge1").text:
    die("03 challenge1 no longer contains ???; reconsider skipping it")

# KIT EXTRA R code run after a chunk in the same session; never prints to the learner console.
# kit_anchor(name, value) writes records/anchors/<step>.tsv; kit_handoff(file) keeps a copy of a
# saved handoff RDS (publish_rds.sh publishes these, so episode 8's overwrite of ifnb_annotated.rds
# cannot leak into the staged copy).
EXTRAS = {
    ("03", "check-dims"): 'kit_anchor("ep3_read_genes", nrow(pbmc.data)); kit_anchor("ep3_read_cells", ncol(pbmc.data))',
    ("03", "create-seurat"): 'kit_anchor("ep3_seurat_genes", nrow(pbmc)); kit_anchor("ep3_seurat_cells", ncol(pbmc))',
    ("03", "calc-mt"): 'kit_anchor("ep3_median_percent_mt_prefilter", round(median(pbmc$percent.mt), 2))',
    ("03", "filter"): 'kit_anchor("ep3_cells_after_qc", ncol(pbmc)); kit_anchor("ep3_genes_after_qc", nrow(pbmc))',
    ("03", "save"): 'kit_handoff("pbmc_filtered.rds")',
    ("04", "find-variable"): 'kit_anchor("ep4_n_variable", length(VariableFeatures(pbmc))); kit_anchor("ep4_top10_variable", paste(head(VariableFeatures(pbmc), 10), collapse = ","))',
    ("04", "save"): 'kit_handoff("pbmc_normalized.rds")',
    ("05", "find-clusters"): 'kit_anchor("ep5_clusters_res0.5", nlevels(pbmc$seurat_clusters)); kit_anchor("ep5_cells", ncol(pbmc))',
    ("05", "multi-resolution"): 'for (r in c("0.2", "0.5", "0.8", "1.2")) kit_anchor(paste0("ep5_clusters_res", r), length(unique(pbmc[[paste0("RNA_snn_res.", r)]][, 1])))',
    ("05", "save"): 'kit_handoff("pbmc_clustered.rds")',
    ("06", "store-labels"): 'kit_marker_check(pbmc, "manual_annotation", "ep6")',
    ("06", "compare"): 'kit_table(table(Manual = pbmc$manual_annotation, SingleR = pbmc$singler_labels), "ep6_singler_vs_manual"); kit_table(table(cluster = pbmc$RNA_snn_res.0.5, SingleR = pbmc$singler_labels), "ep6_cluster_vs_singler")',
    ("06", "challenge1"): 'kit_label_profile(pbmc, c("CD3D", "SELL", "CCR7", "CD69", "IL7R"), "ep6_challenge1")',
    ("06", "save"): 'kit_handoff("pbmc_annotated.rds"); kit_anchor("ep6_cells", ncol(pbmc)); kit_anchor("ep6_n_labels", nlevels(Idents(pbmc)))',
    ("07", "load-ifnb"): 'kit_anchor("ep7_ifnb_cells", ncol(ifnb)); kit_anchor("ep7_ifnb_genes", nrow(ifnb))',
    ("07", "check-stim"): 'kit_anchor("ep7_ctrl_cells", sum(ifnb$stim == "CTRL")); kit_anchor("ep7_stim_cells", sum(ifnb$stim == "STIM"))',
    ("07", "post-integrate"): 'kit_anchor("ep7_clusters_integrated", nlevels(ifnb$seurat_clusters))',
    ("07", "de-condition"): 'kit_top_genes(mono.de, "ep7_mono_top20")',
    ("07", "save.obj"): 'kit_handoff("ifnb_annotated.rds"); kit_anchor("ep7_n_celltypes", length(unique(ifnb$celltype)))',
    ("08", "wilcox-count"): 'kit_anchor("ep8_wilcox_sig", sum(mono.wilcox$p_val_adj < 0.05, na.rm = TRUE)); kit_top_genes(mono.wilcox, "ep8_wilcox_top20")',
    ("08", "pseudo-count"): 'kit_anchor("ep8_pseudobulk_sig", sum(mono.pseudo$p_val_adj < 0.05, na.rm = TRUE)); kit_anchor("ep8_aggregate_class", class(bulk)[1]); kit_top_genes(mono.pseudo, "ep8_pseudo_top20")',
    ("08", "enrichgo"): 'kit_anchor("ep8_go_top15", paste(head(as.data.frame(go_results)$Description, 15), collapse = "|"))',
    ("08", "kegg"): 'kit_anchor("ep8_kegg_top5", paste(head(as.data.frame(kegg_results)$Description, 5), collapse = "|"))',
}
for (ep, lab) in EXTRAS:
    chunk(ep, lab)

# Where episode 3 reads the count matrix (relative to scrna_workshop/), from its set-paths chunk.
# It decides whether r03 waits for 02-cellranger and where the cold pass seeds the staged matrix.
m3 = re.search(r'data_dir\s*<-\s*paste0\(.*?"/scrna_workshop/([^"]*)"', chunk("03", "set-paths").text, re.S)
if not m3:
    die("03 set-paths: cannot find data_dir <- paste0(..., \"/scrna_workshop/...\")")
EP3_MATRIX = m3.group(1)
EP3_NEEDS_CR = EP3_MATRIX.startswith("cellranger_output/")

PLANS = {}


def plan_rows(ep, labels=None, extra_before=None):
    rows = []
    if extra_before:
        rows.append((ep, "kit-pre", "extras/%s/kit-pre.R" % ep, "extra", "0", "0", "-"))
    for c in EPS[ep].r_chunks():
        if labels is not None and c.label not in labels:
            continue
        opt = "1" if any(x in c.context for x in OPTIONAL_CTX) or c.label.startswith(("challenge", "solution")) else "0"
        act = "run"
        if labels is None and (ep, c.label) in R_SKIP:
            act = "skip:" + R_SKIP[(ep, c.label)]
        rows.append((ep, c.label, "blocks/%s/%s.R" % (ep, c.label), act, opt, str(c.start), c.context))
        if (ep, c.label) in EXTRAS and labels is None:
            rows.append((ep, c.label + "+kit", "extras/%s/%s.R" % (ep, c.label), "extra", "0", str(c.start), "-"))
    return rows


for ep in R_EPS:
    PLANS[ep] = plan_rows(ep)

# ---------------------------------------------------------------- figures and expected outputs
FIG_RE_NAME = re.compile(r'fig/([\w.-]+\.png)$')


def figmap_rows():
    rows = []
    for ep in R_EPS:
        e = EPS[ep]
        code = [c for c in e.chunks if c.lang in ("r", "bash")]
        for i, c in enumerate(code):
            if c.lang != "r" or not c.fence.startswith("```{"):
                continue
            nxt = code[i + 1].start if i + 1 < len(code) else 10 ** 9
            # code-generated figures only (static diagrams such as seurat_object.png are not drawn by code)
            figs = [f for (ln, f, alt, cap) in e.figures if c.end < ln < nxt and re.search(r'fig/(ep\d_|challenge)', f)]
            for k, f in enumerate(figs, 1):
                rows.append((ep, c.label, k, FIG_RE_NAME.search(f).group(1) if FIG_RE_NAME.search(f) else f))
    return rows


def expected_outputs():
    out = {}
    for ep in R_EPS:
        e = EPS[ep]
        r = [c for c in e.chunks if c.lang == "r" and c.fence.startswith("```{")]
        for i, c in enumerate(r):
            nxt = r[i + 1].start if i + 1 < len(r) else 10 ** 9
            outs = [o for o in e.chunks if o.lang in ("output", "plain") and c.end < o.start < nxt]
            if outs:
                out[(ep, c.label)] = "".join("# %s:%d\n%s" % (rel(ep), o.start, o.text) for o in outs)
    return out


# ---------------------------------------------------------------- packages used by the R code
def r_packages(by_ep=None):
    pk = set()
    for ep in R_EPS:
        before = set(pk)
        for c in EPS[ep].r_chunks():
            code = "\n".join(l.split("#")[0] for l in c.text.splitlines())
            pk.update(re.findall(r'\blibrary\(([\w.]+)\)', code))
            pk.update(re.findall(r'\brequire(?:Namespace)?\("?([\w.]+)"?', code))
            pk.update(re.findall(r'\b([A-Za-z][\w.]*)::', code))
            pk.update(re.findall(r'install\("([\w.]+)"\)', code))
            if re.search(r'test\.use\s*=\s*"DESeq2"', code):
                pk.add("DESeq2")
            if re.search(r'test\.use\s*=\s*"MAST"', code):
                pk.add("MAST")
            if "FindConservedMarkers" in code:
                pk.update(["metap", "multtest"])  # Seurat suggests them; FindConservedMarkers stops without metap
            if re.search(r'OrgDb\s*=\s*org\.Hs\.eg\.db', code):
                pk.add("org.Hs.eg.db")
        if by_ep is not None:
            for p in pk:
                if p not in before or any(re.search(r'\b%s\b' % re.escape(p), c.text) for c in EPS[ep].r_chunks()):
                    by_ep.setdefault(p, set()).add(ep)
    pk.discard("BiocManager")  # used, but only to install
    return sorted(pk | {"BiocManager", "BiocVersion", "SeuratObject", "Matrix"})


def r_package_table():
    """Derived packages (blocking for the episodes that use them) merged with package_requirements.tsv."""
    by_ep = {}
    names = r_packages(by_ep)
    req = {}
    for line in (KIT / "package_requirements.tsv").read_text().splitlines():
        if line.startswith("#") or line.startswith("package\t") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) != 5:
            die("package_requirements.tsv: bad row %r" % line)
        req[f[0]] = f
    rows = [("package", "episodes", "blocking", "warn_if_version_le", "reason")]
    for p in sorted(set(names) | set(req)):
        if p in req:
            rows.append(tuple(req[p]))
        else:
            eps = ",".join(sorted(by_ep.get(p, []))) or "all"
            rows.append((p, eps, "yes", "-", "used by the episode code (library(), ::, install() or test.use)"))
    return rows


# ---------------------------------------------------------------- step table
# kind: login   run by submit_all.sh on the login node before submitting (learners run these there)
#       learner the learner's own SLURM script, as written by the episode's heredoc, submitted from $W
#       kitjob  batch job running generated learner-shell sessions
#       rjob    R episode (or kit R check) in the configured R environment, OOD form resources
STEPS = []


def step(name, kind, deps, desc, dep_type="afterok", script="-", ntasks="", time=""):
    STEPS.append(dict(name=name, kind=kind, deps=deps, desc=desc, dep_type=dep_type, script=script, ntasks=ntasks, time=time))


step("02-prep", "login", [], "02 directories, rsync, whitelist, listings, write the three job scripts (login node)")
step("02-cellranger", "learner", [], "02 Cell Ranger count (run_cellranger.sh)", script="run_cellranger.sh",
     ntasks=H_NTASKS["run_cellranger.sh"], time=H_TIME["run_cellranger.sh"])
step("02-starindex", "learner", [], "02 STAR index (build_star_index.sh)", script="build_star_index.sh",
     ntasks=H_NTASKS["build_star_index.sh"], time=H_TIME["build_star_index.sh"])
step("02-starsolo", "learner", ["02-starindex"], "02 STARsolo (run_starsolo.sh)", script="run_starsolo.sh",
     ntasks=H_NTASKS["run_starsolo.sh"], time=H_TIME["run_starsolo.sh"])
step("02-compare", "kitjob", ["02-cellranger", "02-starsolo"], "02 compare outputs, metrics, barcodes vs staged",
     dep_type="afterany", ntasks=1, time="1:00:00")
step("env-egress", "rjob", [], "R environment, package versions, per-call egress from a compute node")


def rds_io(ep):
    """(reads, saves): the .rds names episode ep's R chunks pass to readRDS() / saveRDS(), comments excluded."""
    reads, saves = set(), set()
    for c in EPS[ep].r_chunks():
        code = "\n".join(l for l in c.text.splitlines() if not l.lstrip().startswith("#"))
        reads |= set(re.findall(r'readRDS\([^;]*?"(?:[^"]*/)?([\w.-]+\.rds)"', code))
        saves |= set(re.findall(r'saveRDS\([^;]*?"(?:[^"]*/)?([\w.-]+\.rds)"', code))
    return reads, saves


# Each continuous step waits for the steps that save the RDS files its episode reads (episode 8
# reads episode 6's pbmc_annotated.rds as well as episode 7's ifnb_annotated.rds).
RDS_IO = {ep: rds_io(ep) for ep in R_EPS}
for ep in R_EPS:
    deps = ["02-cellranger"] if ep == "03" and EP3_NEEDS_CR else []
    for f in sorted(RDS_IO[ep][0]):
        src = [e for e in R_EPS if e < ep and f in RDS_IO[e][1]]
        if src:
            deps.append("r" + src[-1])
        elif f not in RDS_IO[ep][1]:
            die("%s reads %s, which no earlier episode saves" % (rel(ep), f))
    note = "; starts from SeuratData" if ep == "07" else ""
    step("r" + ep, "rjob", sorted(set(deps)), "%s continuous (%s%s)" % (ep, ORDER[int(ep) - 1], note))
for ep in ("03", "04", "05", "06", "08"):
    step("c" + ep, "rjob", [], "%s cold start from the staged copies (%s)" % (ep, ORDER[int(ep) - 1]))

# What a cold start copies from $STAGED/rds/: every RDS the episode reads and does not save itself
# (episode 8 needs pbmc_annotated.rds as well as ifnb_annotated.rds). Each must be a staged rds/ file.
STAGED_RDS = {l.split("\t")[0][4:] for l in (KIT / "staged_manifest.tsv").read_text().splitlines()
              if l.startswith("rds/") and "\tfile\t" in l}
COLD_SEED = []
for ep in R_EPS:
    for f in sorted(RDS_IO[ep][0] - RDS_IO[ep][1]):
        if f not in STAGED_RDS:
            die("%s reads %s, which staged_manifest.tsv does not list under rds/ (cold start)" % (rel(ep), f))
        COLD_SEED.append((ep, f))

# ---------------------------------------------------------------- writers
OUT = {}


def put(relpath, text):
    OUT[relpath] = text if text.endswith("\n") else text + "\n"


def sha256(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def render_session(name, pieces):
    """Learner-like shell: no -e/-u/pipefail; failing commands are trapped and logged, not fatal."""
    out = ["# shellcheck shell=bash disable=SC2164  # learner code verbatim: plain cd, as typed",
           "# GENERATED by build_kit.py; learner session for step %s." % name,
           "# Run by kit_run_session in a child bash without set -euo pipefail, like a learner's shell.",
           "# Lmod's module function, as in a learner's login shell (sbatch normally exports it)",
           'type module >/dev/null 2>&1 || source "${LMOD_PKG:-/opt/lmod/lmod}/init/bash" 2>/dev/null || source /etc/profile.d/lmod.sh',
           "trap 'echo \"### KIT-ERR rc=$? line=$LINENO cmd=$BASH_COMMAND\" >&2' ERR",
           'cd "$RCAC_SCRATCH"']
    for kind, a, b in pieces:
        if kind == "block":
            out.append('echo "### KIT-BLOCK %s:%d begin $(date +%%s)"' % (rel(a), b.start))
            out.append(b.text.rstrip("\n"))
            out.append('echo "### KIT-BLOCK %s:%d end $(date +%%s)"' % (rel(a), b.start))
        else:
            out.append("# KIT EXTRA: %s" % a)
            out.append("echo '### KIT-EXTRA begin' \"$(date +%%s)\" '%s'" % a.replace("'", ""))
            out.append(b)
            out.append('echo "### KIT-EXTRA end $(date +%s)"')
    return "\n".join(out)


def header(name, desc, ntasks, time, account_note, account=KIT_ACCOUNT):
    return ["#!/bin/bash",
            "# GENERATED by build_kit.py; do not edit. Step: %s" % name,
            "# %s" % desc,
            "# %s" % account_note,
            "#SBATCH --job-name=kit-%s" % name,
            "#SBATCH --account=%s" % account,
            "#SBATCH --partition=cpu",
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks=%s" % ntasks,
            "#SBATCH --time=%s" % time,
            "set -euo pipefail",
            'source "$KIT_DIR/lib.sh"']


def render_job(d):
    n = d["name"]
    if d["kind"] == "kitjob":
        lines = header(n, d["desc"], d["ntasks"], d["time"], "Account from config.sh KIT_ACCOUNT (passed by submit_all.sh).")
        lines += ['kit_job_begin "%s"' % n, 'kit_run_session "%s" "$KIT_GEN/sessions/%s.sh"' % (n, n), 'kit_job_end "%s"' % n]
    else:
        note = "Resources match the OOD form in learners/setup.md; submit_all.sh passes --account=$OOD_ACCOUNT, --ntasks and --time from config.sh."
        lines = header(n, d["desc"], 16, "8:00:00", note, OOD_ACCOUNT)
        lines += ['kit_job_begin "%s"' % n]
        if n == "env-egress":
            lines += ['log="$RUN/logs/R-env-egress.${SLURM_JOB_ID:-$$}.log"',
                      'rc=0; kit_r "$KIT_R_MODE" env-egress "$KIT_USER_OVERRIDE" "$log" "$KIT_DIR/R/env_egress.R" || rc=$?',
                      'cp "$log" "$RUN/logs/R-env-egress.log"; cp "$log" "$RUN/records/console/R-env-egress.log"',
                      'if [ $rc -eq 0 ]; then kit_job_end env-egress PASS; else kit_job_end env-egress WARN; fi',
                      '# egress failures are findings, not a failed step: dependents do not wait on this step']
        else:
            pass_ = "continuous" if n.startswith("r") else "cold"
            lines += ['kit_r_episode "$KIT_R_MODE" %s %s' % (pass_, n[1:])]
    return "\n".join(lines)


def tsv(rows):
    return "\n".join("\t".join(str(x) for x in r) for r in rows)


def build():
    sums = []
    for name in ORDER:
        p = EPISODES / (name + ".Rmd")
        sums.append("%s  %s" % (sha256(p), p.relative_to(ROOT).as_posix()))
    sums.append("%s  %s" % (sha256(SETUP), SETUP.relative_to(ROOT).as_posix()))
    put("episodes.sha256", "\n".join(sums))

    # every block, verbatim (R chunks by label, shell blocks by number), plus the index
    idx = [("episode", "nn", "lang", "label", "start", "end", "context", "file")]
    for ep in sorted(EPS):
        for c in EPS[ep].chunks:
            if c.lang == "r" and c.fence.startswith("```{"):
                f = "blocks/%s/%s.R" % (ep, c.label)
            elif c.lang == "bash":
                f = "blocks/%s/%02d.sh" % (ep, c.nn)
            else:
                f = "-"
            if f != "-":
                put(f, c.text)
            idx.append((ep, c.nn, c.lang, c.label or "-", c.start, c.end, c.context, f))
    put("index.tsv", tsv(idx))

    for name, (b, f, body) in SCRIPTS.items():
        put("learner-scripts/" + name, body)
    put("learner-scripts/cellranger_array_template.sh", B_ARRAY.text)
    put("learner-header-check.tsv", "\n".join(HEADER_ROWS))
    put("submit-blocks/02-cellranger.sh", B2["sub_cr"].text)
    put("submit-blocks/02-starindex.sh", B2["sub_idx"].text)
    put("submit-blocks/02-starsolo.sh", B2["sub_solo"].text)

    for name, pieces in S.items():
        put("sessions/%s.sh" % name, render_session(name, pieces))
    put("setup/01-copy.sh", bS_copy.text)
    put("setup/02-verify.sh", bS_verify.text)

    hdr = ("episode", "label", "file", "action", "optional", "line", "context")
    for name, rows in sorted(PLANS.items()):
        put("rplans/%s.tsv" % name, tsv([hdr] + rows))
    for (ep, lab), code in sorted(EXTRAS.items()):
        put("extras/%s/%s.R" % (ep, lab), code)

    for d in STEPS:
        if d["kind"] in ("kitjob", "rjob"):
            put("jobs/%s.sh" % d["name"], render_job(d))
    st = [("step", "kind", "deps", "dep_type", "ntasks", "time", "script", "description")]
    for d in STEPS:
        st.append((d["name"], d["kind"], ",".join(d["deps"]) or "-", d["dep_type"], d["ntasks"] or "-", d["time"] or "-",
                   d["script"], d["desc"]))
    put("steps.tsv", tsv(st))

    put("figmap.tsv", tsv([("episode", "label", "k", "figure")] + figmap_rows()))
    exp = expected_outputs()
    for (ep, lab), text in sorted(exp.items()):
        put("expected/%s/%s.txt" % (ep, lab), text)
    pt = r_package_table()
    put("r_packages.tsv", tsv(pt))
    put("r_packages.txt", "\n".join(r[0] for r in pt[1:]))
    put("ep3_matrix_path.txt", EP3_MATRIX)
    put("cold_seed.tsv", tsv([("episode", "rds")] + COLD_SEED))
    fm = [("episode", "file", "title", "teaching", "exercises")]
    for name in ORDER:
        e = EPS[KEY[name]]
        fm.append((KEY[name], rel(KEY[name]), e.front.get("title", "").strip("'\""), e.front.get("teaching", ""),
                   e.front.get("exercises", "")))
    put("frontmatter.tsv", tsv(fm))


def write_all(check=False):
    if check:
        diffs = []
        for r, text in sorted(OUT.items()):
            p = GEN / r
            if not p.exists() or p.read_text() != text:
                diffs.append("generated/" + r)
        have = {p.relative_to(GEN).as_posix() for p in GEN.rglob("*") if p.is_file()} if GEN.exists() else set()
        diffs += ["generated/%s (stale)" % r for r in sorted(have - set(OUT))]
        if diffs:
            print("build_kit --check: FAIL, differs from a fresh build:\n  " + "\n  ".join(diffs[:40]))
            sys.exit(1)
        print("build_kit --check: generated/ matches a fresh build (%d files)" % len(OUT))
        return
    if GEN.exists():
        shutil.rmtree(str(GEN))
    for r, text in sorted(OUT.items()):
        p = GEN / r
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
        if r.endswith(".sh"):
            p.chmod(0o755)
    print("build_kit: wrote %d files to %s" % (len(OUT), GEN.relative_to(ROOT)))


if __name__ == "__main__":
    build()
    write_all(check="--check" in sys.argv[1:])
