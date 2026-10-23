#!/usr/bin/env python3
"""s3etag.py FILE PART_MIB: print the S3 multipart ETag of FILE ("<md5 of part md5s>-<parts>") and its md5.

cf.10xgenomics.com serves the 10x archives from S3 with multipart ETags; sources.tsv records the
ETag and the part size (8 MiB), so a download can be verified without a published md5.
Output: etag<TAB>md5. Python 3.6 compatible.
"""
import hashlib
import sys

path, part = sys.argv[1], int(sys.argv[2]) * 1024 * 1024
whole, parts = hashlib.md5(), []
with open(path, "rb") as fh:
    while True:
        b = fh.read(part)
        if not b:
            break
        whole.update(b)
        parts.append(hashlib.md5(b).digest())
etag = hashlib.md5(b"".join(parts)).hexdigest() + "-%d" % len(parts) if len(parts) > 1 else whole.hexdigest()
print("%s\t%s" % (etag, whole.hexdigest()))
