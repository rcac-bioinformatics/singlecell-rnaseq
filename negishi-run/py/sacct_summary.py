#!/usr/bin/env python3
"""sacct_summary.py SACCT_TXT JOBS_TSV: queue wait, elapsed, CPU efficiency and peak MaxRSS per job.

SACCT_TXT is `sacct -P -n --format=JobID,JobName%40,Elapsed,MaxRSS,ReqMem,State,Submit,Start,End,
AllocCPUS,Timelimit,ExitCode,TotalCPU` output (pipe-separated, no header); JOBS_TSV is $RUN/jobs.tsv
(step, jobid, submitted, deps). Queue wait is measured from submission or, for a job that waited on
dependencies, from the end of its last dependency, whichever is later (the variant-calling kit
counted dependency time as queue time). MaxRSS is the maximum over the job's steps.
Writes TSV to stdout. Python 3.6 compatible.
"""
import sys
from datetime import datetime

COLS = ["JobID", "JobName", "Elapsed", "MaxRSS", "ReqMem", "State", "Submit", "Start", "End", "AllocCPUS",
        "Timelimit", "ExitCode", "TotalCPU"]


def seconds(t):
    if not t or t in ("Unknown", "INVALID"):
        return None
    days = 0
    if "-" in t:
        d, t = t.split("-", 1)
        days = int(d)
    parts = [float(p) for p in t.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0.0)
    h, m, s = parts
    return days * 86400 + h * 3600 + m * 60 + s


def gib(rss):
    if not rss:
        return None
    units = {"K": 1 / 1048576.0, "M": 1 / 1024.0, "G": 1.0, "T": 1024.0}
    return float(rss[:-1]) * units[rss[-1]] if rss[-1] in units else float(rss) / 1073741824.0


def when(t):
    try:
        return datetime.strptime(t, "%Y-%m-%dT%H:%M:%S")
    except (TypeError, ValueError):
        return None


def main(sacct_txt, jobs_tsv):
    jobs = {}
    for line in open(sacct_txt):
        f = line.rstrip("\n").split("|")
        if len(f) < len(COLS) - 1:
            continue
        r = dict(zip(COLS, f))
        jid = r["JobID"]
        base = jid.split(".")[0]
        if "." not in jid:
            jobs[base] = dict(r)
            jobs[base]["_rss"] = gib(r.get("MaxRSS", ""))
        elif base in jobs:
            rss = gib(r.get("MaxRSS", ""))
            if rss is not None and (jobs[base]["_rss"] is None or rss > jobs[base]["_rss"]):
                jobs[base]["_rss"] = rss
    deps = {}
    step_of = {}
    for line in open(jobs_tsv):
        f = line.rstrip("\n").split("\t")
        if f[0] == "step" or len(f) < 4:
            continue
        step_of[f[1]] = f[0]
        deps[f[1]] = [d for d in f[3].split(":") if d]
    print("\t".join(["jobid", "step", "name", "state", "cpus", "wait_s", "elapsed_s", "cpu_efficiency", "maxrss_gib", "timelimit", "exitcode"]))
    for jid, r in sorted(jobs.items()):
        sub, start = when(r.get("Submit")), when(r.get("Start"))
        ready = sub
        for d in deps.get(jid, []):
            e = when(jobs.get(d, {}).get("End"))
            if e and (ready is None or e > ready):
                ready = e
        wait = (start - ready).total_seconds() if ready and start else None
        el = seconds(r.get("Elapsed"))
        cpu = seconds(r.get("TotalCPU"))
        n = int(r.get("AllocCPUS") or 0)
        eff = cpu / (el * n) if el and n and cpu is not None else None
        rss = r["_rss"]
        print("\t".join(str(x) for x in [jid, step_of.get(jid, "-"), r.get("JobName", ""), r.get("State", ""), n,
                                         "" if wait is None else int(max(wait, 0)), "" if el is None else int(el),
                                         "" if eff is None else "%.2f" % eff, "" if rss is None else "%.2f" % rss,
                                         r.get("Timelimit", ""), r.get("ExitCode", "")]))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
