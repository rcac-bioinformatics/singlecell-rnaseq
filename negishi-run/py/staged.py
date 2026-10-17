#!/usr/bin/env python3
"""staged.py: the staged learner copy and its manifest (used by check_staged.sh, publish_staged.sh,
make_staged_manifest.py). Python 3.6 compatible.

  staged.py check DIR MANIFEST [--quick] [--perm group|other]
      every manifest row present, readable, expected size and md5 (when frozen); directories
      traversable; nothing extra except inside "tree" rows. Exit 1 on any FAIL.
  staged.py plan SRC DEST
      what copying SRC over DEST would do, one row per file: add | replace | same (md5 compared
      when sizes agree). Never writes.
  staged.py record DIR MANIFEST OUT [NOTE]
      measure DIR against MANIFEST's rows and write a frozen manifest (bytes and md5 filled in).

Manifest columns: path type bytes md5 need source; type file | dir | tree (a directory whose
contents are not listed one by one); bytes and md5 "-" until frozen; need required | optional.
"""
import hashlib
import os
import stat
import sys

# files that may be present without being a FAIL, and why (reported as WARN)
TOLERATED = {"reference/3M-february-2018.txt": "learners extract it in episode 2 (the staged copy is overwritten); remove it from the staged copy"}


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def read_manifest(path):
    rows = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if f[0] == "path":
            continue
        rows.append(dict(path=f[0], type=f[1], bytes=f[2], md5=f[3], need=f[4], source=f[5] if len(f) > 5 else ""))
    return rows


def readable(st, perm):
    return bool(st.st_mode & (stat.S_IRGRP if perm == "group" else stat.S_IROTH))


def traversable(st, perm):
    return bool(st.st_mode & (stat.S_IXGRP if perm == "group" else stat.S_IXOTH))


def check(d, manifest, quick=False, perm="group"):
    bad = 0
    out = []

    def fail(msg):
        nonlocal bad
        bad += 1
        out.append("FAIL  " + msg)

    rows = read_manifest(manifest)
    trees = [r["path"] for r in rows if r["type"] == "tree"]
    listed = set()
    for r in rows:
        p = os.path.join(d, r["path"])
        listed.add(r["path"])
        if r["type"] in ("dir", "tree"):
            if os.path.isdir(p):
                out.append("ok    %s/" % r["path"])
            elif r["need"] == "optional":
                out.append("WARN  %s/ absent (optional; %s)" % (r["path"], r["source"]))
            else:
                fail("%s/ missing (%s)" % (r["path"], r["source"]))
            continue
        if not os.path.isfile(p):
            if r["need"] == "optional":
                out.append("WARN  %s absent (optional; %s)" % (r["path"], r["source"]))
            else:
                fail("%s missing (%s)" % (r["path"], r["source"]))
            continue
        st = os.stat(p)
        if not readable(st, perm):
            fail("%s not readable by %s (%s)" % (r["path"], perm, stat.filemode(st.st_mode)))
            continue
        if r["bytes"] != "-" and str(st.st_size) != r["bytes"]:
            fail("%s size %d, expected %s" % (r["path"], st.st_size, r["bytes"]))
            continue
        if r["md5"] == "-":
            out.append("ok    %s (%d bytes; size and md5 not frozen yet)" % (r["path"], st.st_size))
            continue
        if quick and st.st_size > 50000000:
            out.append("ok    %s (size only)" % r["path"])
            continue
        m = md5(p)
        if m == r["md5"]:
            out.append("ok    %s" % r["path"])
        else:
            fail("%s md5 %s, expected %s" % (r["path"], m, r["md5"]))
    # directories traversable, nothing extra, nothing unreadable inside trees
    extra = []
    for root, dirs, files in os.walk(d):
        rel = os.path.relpath(root, d)
        st = os.stat(root)
        if not traversable(st, perm):
            fail("directory %s not traversable by %s (%s)" % (root, perm, stat.filemode(st.st_mode)))
        for f in files:
            rp = os.path.normpath(os.path.join(rel, f))
            in_tree = any(rp.startswith(t + "/") for t in trees)
            if in_tree:
                if not readable(os.stat(os.path.join(root, f)), perm):
                    fail("%s not readable by %s" % (rp, perm))
                continue
            if rp not in listed:
                extra.append(rp)
    for rp in sorted(extra):
        if rp in TOLERATED:
            out.append("WARN  %s present: %s" % (rp, TOLERATED[rp]))
        else:
            fail("%s is not in the manifest (learners would receive it)" % rp)
    # path to the copy
    p = os.path.realpath(d)
    while p != "/":
        st = os.stat(p)
        if not traversable(st, perm):
            fail("directory %s on the path is not traversable by %s (%s)" % (p, perm, stat.filemode(st.st_mode)))
        p = os.path.dirname(p)
    out.append("")
    out.append("%s: %d problem(s)" % (d, bad))
    print("\n".join(out))
    return bad


def walk_files(d):
    for root, dirs, files in os.walk(d):
        dirs.sort()
        for f in sorted(files):
            p = os.path.join(root, f)
            yield os.path.relpath(p, d), p


def plan(src, dest):
    print("action\tpath\tsrc_bytes\tdest_bytes\tsrc_md5\tdest_md5")
    for rel, p in walk_files(src):
        q = os.path.join(dest, rel)
        sb = os.path.getsize(p)
        if not os.path.exists(q):
            print("add\t%s\t%d\t-\t-\t-" % (rel, sb))
            continue
        db = os.path.getsize(q)
        if sb != db:
            print("replace\t%s\t%d\t%d\t-\t-" % (rel, sb, db))
            continue
        ms, md = md5(p), md5(q)
        print("%s\t%s\t%d\t%d\t%s\t%s" % ("same" if ms == md else "replace", rel, sb, db, ms, md))


def record(d, manifest, outp, note=""):
    rows = read_manifest(manifest)
    lines = ["# Expected content of the staged learner copy ($STAGED). Generated by make_staged_manifest.py from",
             "# %s%s; do not edit by hand. Anything else in $STAGED (episode 2 outputs, job logs) is an error." % (d, ("; " + note) if note else ""),
             "path\ttype\tbytes\tmd5\tneed\tsource"]
    for r in rows:
        p = os.path.join(d, r["path"])
        if r["type"] == "file" and os.path.isfile(p):
            r["bytes"], r["md5"] = str(os.path.getsize(p)), md5(p)
            print("%-70s %s" % (r["path"], r["md5"]))
        lines.append("\t".join([r["path"], r["type"], r["bytes"], r["md5"], r["need"], r["source"]]))
    open(outp, "w").write("\n".join(lines) + "\n")
    print("wrote %s" % outp)


if __name__ == "__main__":
    a = sys.argv[1:]
    if not a:
        sys.exit(__doc__)
    if a[0] == "check":
        perm = a[a.index("--perm") + 1] if "--perm" in a else "group"
        sys.exit(1 if check(a[1], a[2], quick="--quick" in a, perm=perm) else 0)
    if a[0] == "plan":
        plan(a[1], a[2])
    elif a[0] == "record":
        record(a[1], a[2], a[3], a[4] if len(a) > 4 else "")
    else:
        sys.exit(__doc__)
