#!/usr/bin/env python3
"""make_staged_manifest.py DIR [OUT] [--note TEXT]: freeze staged_manifest.tsv from a verified copy.

Measures DIR (normally $STAGED on Negishi, after check_staged.sh passes and the staged files are
the ones the lesson's outputs were produced from) against the rows of the current
staged_manifest.tsv and writes a manifest with every file's size and md5 (default OUT:
staged_manifest.measured.tsv next to this script; bring it back and rename it to
staged_manifest.tsv in round 3). Python 3.6 compatible.
"""
import pathlib
import sys

KIT = pathlib.Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(KIT / "py"))
import staged  # noqa: E402

a = sys.argv[1:]
note = ""
if "--note" in a:
    i = a.index("--note"); note = a[i + 1]; del a[i:i + 2]
if not a:
    sys.exit(__doc__)
out = a[1] if len(a) > 1 else str(KIT / "staged_manifest.measured.tsv")
staged.record(a[0], str(KIT / "staged_manifest.tsv"), out, note)
