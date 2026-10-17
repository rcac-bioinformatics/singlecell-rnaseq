#!/usr/bin/env python3
"""static_checks.py: round 1 static checks of the lesson (laptop, no data, no cluster).

Usage: python3 negishi-run/tools/static_checks.py [--links] [--out DIR]
  --links   also request every external URL (HEAD, then GET on 4xx/5xx); needs network
  --out     where to write the report (default negishi-run/records/static/)

Checks, one section each in static-report.md:
  handoff      saveRDS/readRDS file names chain across episodes 3 to 8
  set-paths    set-paths chunks identical, work_dir ends with "/", setwd() placement
  names        object and metadata column names used per episode, and objects used
               before any assignment or readRDS in that episode (cold start)
  slurm        episode 2 SLURM headers: account, partition, nodes, ntasks, no --mem
  frontmatter  teaching/exercises in each episode vs instructors/instructor-notes.md
  figures      referenced figures exist, have alt text, are described in
               episodes/fig/FIGURE_DESCRIPTIONS.md; unreferenced figures
  code         smart quotes, typographic dashes, shell/R prompts, hard wraps in code;
               bash -n on every bash fence; R parse() is done by build_kit.py's R check
  outputs      output fence spellings (```output, ```{output}, bare ```)
  versions     versions asserted in prose, for comparison with the environment check
  links        (with --links) HTTP status of every external URL
Exit status 0 always; the report is the result.
"""
import pathlib
import re
import subprocess
import sys
import tempfile
import concurrent.futures
import urllib.request
import urllib.error

KIT = pathlib.Path(__file__).resolve().parent.parent
ROOT = KIT.parent
EP = ROOT / "episodes"
ORDER = ["introduction", "raw-data-processing", "quality-control", "normalization-and-features",
         "dimensionality-reduction-clustering", "cell-type-annotation", "multi-sample-integration",
         "differential-expression"]
R_EPS = ORDER[2:]

sys.dont_write_bytecode = True
sys.path.insert(0, str(KIT / "py"))
from lessonparse import parse_rmd  # noqa: E402

args = sys.argv[1:]
OUTDIR = KIT / "records" / "static"
if "--out" in args:
    OUTDIR = pathlib.Path(args[args.index("--out") + 1])
OUTDIR.mkdir(parents=True, exist_ok=True)

EPS = {name: parse_rmd(EP / (name + ".Rmd")) for name in ORDER}
REPORT = []


def section(title):
    REPORT.append("\n## %s\n" % title)


def say(line=""):
    REPORT.append(line)


def loc(ep, line):
    return "episodes/%s.Rmd:%d" % (ep, line)


# ---------------------------------------------------------------- handoff
section("RDS handoff chain")
SAVE_RE = re.compile(r'saveRDS\(\s*([\w.]+)\s*,\s*(?:file\s*=\s*)?(paste0\([^)]*\)|"[^"]+")')
READ_RE = re.compile(r'readRDS\(\s*(paste0\([^)]*\)|"[^"]+")')
say("| episode | loads | saves |")
say("|---|---|---|")
chain = {}
for ep in R_EPS:
    loads, saves = [], []
    for c in EPS[ep].chunks:
        if c.lang != "r":
            continue
        for i, l in enumerate(c.text.splitlines()):
            if l.lstrip().startswith("#"):
                continue
            for m in SAVE_RE.finditer(l):
                f = re.search(r'"([^"]+\.rds)"', m.group(2))
                saves.append((m.group(1), f.group(1) if f else m.group(2), c.label, c.start + i + 1, m.group(2)))
            for m in READ_RE.finditer(l):
                f = re.search(r'"([^"]+\.rds)"', m.group(1))
                loads.append((f.group(1) if f else m.group(1), c.label, c.start + i + 1, m.group(1)))
    chain[ep] = (loads, saves)
    say("| %s | %s | %s |" % (ep, "; ".join("%s (%s, :%d)" % (f, lab, ln) for f, lab, ln, _ in loads) or "-",
                              "; ".join("%s <- %s (%s, :%d)" % (f, o, lab, ln) for o, f, lab, ln, _ in saves) or "-"))
say()
expect = [("normalization-and-features", "pbmc_filtered.rds"), ("dimensionality-reduction-clustering", "pbmc_normalized.rds"),
          ("cell-type-annotation", "pbmc_clustered.rds"), ("differential-expression", "ifnb_annotated.rds")]
prev = {"normalization-and-features": "quality-control", "dimensionality-reduction-clustering": "normalization-and-features",
        "cell-type-annotation": "dimensionality-reduction-clustering", "differential-expression": "multi-sample-integration"}
for ep, f in expect:
    loaded = [x[0] for x in chain[ep][0]]
    saved = [x[1] for x in chain[prev[ep]][1]]
    say("- %s loads %s: %s; %s saves it: %s" % (ep, f, "yes" if f in loaded else "NO", prev[ep], "yes" if f in saved else "NO"))
for ep in R_EPS:
    for o, f, lab, ln, expr in chain[ep][1]:
        if "work_dir" not in expr:
            say("- %s saves %s relative to the working directory (depends on setwd(work_dir)): %s" % (loc(ep, ln), f, expr))
    for f, lab, ln, expr in chain[ep][0]:
        if "work_dir" not in expr:
            say("- %s reads %s without work_dir: %s" % (loc(ep, ln), f, expr))
allsaves = [(ep, x) for ep in R_EPS for x in chain[ep][1]]
seen = {}
for ep, (o, f, lab, ln, _) in allsaves:
    if f in seen:
        say("- %s is written twice: %s and %s (the second overwrites the handoff copy)" % (f, seen[f], loc(ep, ln)))
    seen[f] = loc(ep, ln)

# ---------------------------------------------------------------- set-paths
section("set-paths chunks")
sp = {}
for ep in R_EPS:
    c = [c for c in EPS[ep].chunks if c.label == "set-paths"]
    if not c:
        say("- %s: no set-paths chunk" % ep)
        continue
    sp[ep] = c[0]
texts = {}
for ep, c in sp.items():
    texts.setdefault(c.text, []).append(ep)
say("%d distinct set-paths texts across %d episodes" % (len(texts), len(sp)))
for t, eps in texts.items():
    say("\nUsed by: %s\n```r\n%s```" % (", ".join(eps), t))
for ep, c in sp.items():
    m = re.search(r'work_dir\s*<-\s*(paste0\((?:[^()]|\([^()]*\))*\))', c.text, re.S)
    if not m:
        say("- %s: no work_dir <- paste0(...) assignment" % loc(ep, c.start))
    elif not re.search(r'/"\s*\)$', m.group(1)):
        say("- %s: work_dir does not end with \"/\": %s" % (loc(ep, c.start), " ".join(m.group(1).split())))
    else:
        say("- %s: work_dir ends with \"/\"" % ep)
    where = [ch.label for ch in EPS[ep].chunks if ch.lang == "r" and re.search(r'^\s*setwd\(', ch.text, re.M)]
    say("- %s: setwd() in chunk(s): %s" % (ep, ", ".join(where) or "NONE"))
    if "RCAC_SCRATCH" not in c.text:
        say("- %s: set-paths does not use RCAC_SCRATCH (%s)" % (loc(ep, c.start), "uses $USER" if "USER" in c.text else "?"))

# ---------------------------------------------------------------- names
section("Object and metadata column names")
TOKENS = {
    "objects": ["pbmc", "ifnb", "sce", "bulk", "mono", "pbmc.ref", "ifnb.split", "ref", "pred", "sc_obj"],
    "columns": ["percent.mt", "percent_mt", "percent.mito", "pct.mt", "seurat_clusters", "seurat_cluster",
                "cell_type", "celltype", "celltype.stim", "cell_type.stim", "stim", "condition", "orig.ident",
                "singler_labels", "SingleR.labels", "seurat_annotations", "replicate", "sample"],
}
say("Counts of whole-word uses in R chunks (comments excluded).\n")
hdr = "| token | " + " | ".join(e.split("-")[0][:6] for e in R_EPS) + " |"
say(hdr)
say("|---" * (len(R_EPS) + 1) + "|")
for group in ("objects", "columns"):
    for tok in TOKENS[group]:
        pat = re.compile(r'(?<![\w.$])' + re.escape(tok) + r'(?![\w])')
        row = []
        for ep in R_EPS:
            n = 0
            for c in EPS[ep].chunks:
                if c.lang == "r":
                    code = "\n".join(l.split("#")[0] for l in c.text.splitlines())
                    n += len(pat.findall(code))
            row.append(str(n) if n else ".")
        if any(x != "." for x in row):
            say("| `%s` | %s |" % (tok, " | ".join(row)))
say("\nObjects used in an episode before that episode assigns or loads them (a cold start fails):\n")
OBJ = ["pbmc", "ifnb", "sce", "bulk", "mono", "pbmc.markers", "ifnb.markers", "ref", "pred", "counts_matrix"]
for ep in R_EPS:
    assigned = set()
    for c in EPS[ep].chunks:
        if c.lang != "r":
            continue
        for i, l in enumerate(c.text.splitlines()):
            code = re.sub(r'"[^"]*"', '""', l.split("#")[0])
            for o in OBJ:
                if re.search(r'(?<![\w.$])' + re.escape(o) + r'(?![\w.])', code) and o not in assigned:
                    lhs = re.match(r'\s*' + re.escape(o) + r'(\$\w+|\[\[.*?\]\])?\s*(<-|=)', code)
                    if lhs and not lhs.group(1):
                        assigned.add(o)
                        continue
                    if not re.search(r'(?<![\w.$])' + re.escape(o) + r'\s*(<-)', code):
                        say("- %s (%s): `%s` used before assignment: `%s`" % (loc(ep, c.start + i + 1), c.label, o, code.strip()[:100]))
                        assigned.add(o)  # report once
            for m in re.finditer(r'([\w.]+)\s*<-', code):
                assigned.add(m.group(1))
say("\nPipes and the packages that provide them:\n")
for ep in R_EPS:
    libs = set()
    for c in EPS[ep].chunks:
        if c.lang == "r":
            libs.update(re.findall(r'library\(([\w.]+)\)', c.text))
    for c in EPS[ep].chunks:
        if c.lang == "r":
            for i, l in enumerate(c.text.splitlines()):
                if "%>%" in l.split("#")[0] and not ({"dplyr", "magrittr", "tidyverse"} & libs):
                    say("- %s (%s): `%%>%%` but the episode loads only %s" % (loc(ep, c.start + i + 1), c.label, ", ".join(sorted(libs))))

# ---------------------------------------------------------------- slurm
section("Episode 2 SLURM headers")
ep = "raw-data-processing"
say("| line | job-name | account | partition | nodes | ntasks | time | array | mem | other |")
say("|---|---|---|---|---|---|---|---|---|---|")
for c in EPS[ep].chunks:
    if c.lang != "bash" or "#SBATCH" not in c.text:
        continue
    h = {}
    for l in c.text.splitlines():
        m = re.match(r'#SBATCH\s+--([\w-]+)(?:=(\S+))?', l.strip())
        if m:
            h[m.group(1)] = m.group(2) or ""
    other = {k: v for k, v in h.items() if k not in ("job-name", "account", "partition", "nodes", "ntasks", "time", "array",
                                                    "mem", "output", "error")}
    say("| %d | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
        c.start, h.get("job-name", "-"), h.get("account", "MISSING"), h.get("partition", "MISSING"), h.get("nodes", "MISSING"),
        h.get("ntasks", "MISSING"), h.get("time", "-"), h.get("array", "-"), h.get("mem", "-"),
        ", ".join("%s=%s" % kv for kv in other.items()) or "-"))

# ---------------------------------------------------------------- front matter
section("Front matter vs instructors/instructor-notes.md")
notes = (ROOT / "instructors" / "instructor-notes.md").read_text()
rows = re.findall(r'^\|\s*(\d)\s*\|\s*([^|]+?)\s*\|\s*(\d+) min\s*\|\s*(\d+) min\s*\|', notes, re.M)
tt = te = 0
for n, title, teach, ex in rows:
    e = EPS[ORDER[int(n) - 1]]
    t2, x2 = e.front.get("teaching"), e.front.get("exercises")
    ok = str(t2) == teach and str(x2) == ex and e.front.get("title", "").strip("'\"") == title
    tt += int(t2 or 0); te += int(x2 or 0)
    say("- %s %s: front matter %s/%s \"%s\", notes %s/%s \"%s\"%s" % (n, ORDER[int(n) - 1], t2, x2, e.front.get("title"), teach, ex,
                                                                title, "" if ok else "  MISMATCH"))
say("- front matter totals: %d teaching + %d exercises = %d min" % (tt, te, tt + te))

# ---------------------------------------------------------------- figures
section("Figures")
desc = (EP / "fig" / "FIGURE_DESCRIPTIONS.md").read_text()
referenced = set()
for ep in ORDER:
    for line, path, alt, caption in EPS[ep].figures:
        referenced.add(pathlib.Path(path).name)
        p = EP / path
        problems = []
        if not p.exists():
            problems.append("FILE MISSING")
        if not alt or len(alt) < 40:
            problems.append("alt text short or missing (%d chars)" % len(alt or ""))
        name = pathlib.Path(path).name
        mm = re.match(r'ep(\d)_(.+)\.png$', name)
        if mm:
            if not re.search(r'^#### %s\.\d+ %s \(line' % (mm.group(1), re.escape(mm.group(2))), desc, re.M):
                problems.append("no '#### %s.N %s' entry in FIGURE_DESCRIPTIONS.md" % (mm.group(1), mm.group(2)))
        elif name not in desc:
            problems.append("not in FIGURE_DESCRIPTIONS.md")
        if problems:
            say("- %s %s: %s" % (loc(ep, line), path, "; ".join(problems)))
for p in sorted((EP / "fig").glob("*.png")):
    if p.name not in referenced:
        say("- episodes/fig/%s is not referenced by any episode" % p.name)
for num, lab, ln in re.findall(r'^#### (\d)\.\d+ ([\w-]+) \(line (\d+)\)', desc, re.M):
    epn = ORDER[int(num) - 1]
    c = [c for c in EPS[epn].chunks if c.label == lab]
    if not c:
        say("- FIGURE_DESCRIPTIONS.md %s.%s: no chunk `%s` in %s" % (num, lab, lab, epn))
    elif abs(c[0].start - int(ln)) > 2:
        say("- FIGURE_DESCRIPTIONS.md `%s` (ep %s) cites line %s; the chunk is at %s:%d" % (lab, num, ln, epn, c[0].start))
for name in sorted(set(re.findall(r'([\w.-]+\.png)', desc))):
    if not (EP / "fig" / name).exists():
        say("- FIGURE_DESCRIPTIONS.md describes %s, which does not exist" % name)
# code-generated figures and the chunk that draws them
say("\nCode-generated figures (episodes/fig/ep*_*.png and challenge*_*.png) and the chunk before each reference:\n")
for ep in R_EPS:
    for line, path, alt, caption in EPS[ep].figures:
        name = pathlib.Path(path).name
        if re.match(r'(ep\d|challenge)', name):
            before = [c for c in EPS[ep].chunks if c.lang == "r" and c.end < line]
            say("- %s <- %s chunk `%s` (:%d)" % (name, ep, before[-1].label if before else "?", before[-1].start if before else 0))

# ---------------------------------------------------------------- code hazards
section("Code hazards (smart quotes, dashes, prompts, hard wraps, bash -n)")
SMART = re.compile('[‘’“”–— ]')
for ep in ORDER:
    for c in EPS[ep].chunks:
        if c.lang not in ("r", "bash"):
            continue
        lines = c.text.splitlines()
        for i, l in enumerate(lines):
            ln = c.start + i + 1
            if SMART.search(l):
                say("- %s: typographic character %r in code: `%s`" % (loc(ep, ln), SMART.search(l).group(0), l.strip()[:90]))
            if c.lang == "bash" and re.match(r'\s*\$\s', l):
                say("- %s: shell prompt in code: `%s`" % (loc(ep, ln), l.strip()[:90]))
            if c.lang == "r" and re.match(r'\s*>\s', l):
                say("- %s: R prompt in code: `%s`" % (loc(ep, ln), l.strip()[:90]))
            if c.lang == "bash" and i > 0 and re.match(r'\s+--?\w', l) and not lines[i - 1].rstrip().endswith("\\") \
                    and not lines[i - 1].lstrip().startswith("#"):
                say("- %s: continuation line without a trailing backslash on the line before: `%s`" % (loc(ep, ln), l.strip()[:90]))
            if len(l) > 120:
                say("- %s: code line of %d characters" % (loc(ep, ln), len(l)))
        if c.lang == "bash":
            with tempfile.NamedTemporaryFile("w", suffix=".sh", delete=False) as fh:
                fh.write(c.text)
            r = subprocess.run(["bash", "-n", fh.name], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               universal_newlines=True)
            pathlib.Path(fh.name).unlink()
            if r.returncode:
                say("- %s: bash -n: %s" % (loc(ep, c.start), r.stderr.strip().replace(fh.name, "block")))

# ---------------------------------------------------------------- outputs
section("Output fences")
for ep in ORDER:
    kinds = {}
    for c in EPS[ep].chunks:
        if c.lang in ("output", "{output}", "plain"):
            kinds.setdefault(c.fence, []).append(c.start)
    if kinds:
        say("- %s: %s" % (ep, "; ".join("%s x%d (first :%d)" % (k, len(v), v[0]) for k, v in sorted(kinds.items()))))

# ---------------------------------------------------------------- versions
section("Versions asserted in prose and code")
VER = re.compile(r'(Cell ?Ranger|cellranger|STAR(?:solo)?|Seurat(?:Object)?|R |Bioconductor|BiocManager|SingleR|celldex|'
                 r'clusterProfiler|DESeq2|sceasy|refdata-gex|GRCh38|Chromium|v[23] chemistry|r-bioconductor|Apptainer|singularity)'
                 r'[^\n]{0,40}?\b(v?\d+\.\d+(?:\.\d+)?[a-z]?|20\d\d-[A-Z]|v\d)\b')
files = [EP / (e + ".Rmd") for e in ORDER] + [ROOT / "learners" / "setup.md", ROOT / "instructors" / "instructor-notes.md"]
for f in files:
    code_lines = set()
    if f.suffix == ".Rmd":
        for c in EPS[f.stem].chunks:
            if c.lang not in ("bash", "r"):
                code_lines.update(range(c.start, c.end + 1))
    for i, l in enumerate(f.read_text().splitlines(), 1):
        if i in code_lines:
            continue
        for m in VER.finditer(l):
            say("- %s:%d: %s" % (f.relative_to(ROOT), i, m.group(0).strip()))

# ---------------------------------------------------------------- links
if "--links" in args:
    section("External links")
    URL = re.compile(r'https?://[^\s)>\]"\'`]+')
    urls = {}
    for f in [EP / (e + ".Rmd") for e in ORDER] + sorted((ROOT / "learners").glob("*.md")) + \
            sorted((ROOT / "instructors").glob("*.md")) + [ROOT / "index.md", ROOT / "links.md", ROOT / "config.yaml"]:
        for i, l in enumerate(f.read_text().splitlines(), 1):
            for u in URL.findall(l):
                u = u.rstrip(".,;:*")
                urls.setdefault(u, []).append("%s:%d" % (f.relative_to(ROOT), i))

    def probe(u):
        if re.search(r'negishi\.rcac|gateway\.negishi|example\.|localhost', u):
            return u, "skip (cluster or placeholder)"
        for method in ("HEAD", "GET"):
            try:
                req = urllib.request.Request(u, method=method, headers={"User-Agent": "Mozilla/5.0 lesson-link-check"})
                with urllib.request.urlopen(req, timeout=25) as r:
                    return u, str(r.status)
            except urllib.error.HTTPError as e:
                code = str(e.code)
                if method == "GET":
                    return u, code
            except Exception as e:  # noqa: BLE001
                if method == "GET":
                    return u, "ERR %s" % type(e).__name__
        return u, "?"
    with concurrent.futures.ThreadPoolExecutor(12) as ex:
        res = dict(ex.map(probe, sorted(urls)))
    bad = [(u, s) for u, s in res.items() if not s.startswith(("2", "3", "skip"))]
    say("%d URLs, %d not OK\n" % (len(res), len(bad)))
    for u, s in sorted(res.items()):
        mark = "" if s.startswith(("2", "3", "skip")) else "  <-- check"
        say("- %s %s (%s)%s" % (s, u, ", ".join(urls[u][:3]), mark))

out = OUTDIR / "static-report.md"
out.write_text("# Static checks\n" + "\n".join(REPORT) + "\n")
print("static_checks: wrote %s" % out)
