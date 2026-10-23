"""lessonparse.py: parse a Workbench episode (.Rmd/.md) into code blocks, figures and front matter.

Shared by build_kit.py and tools/static_checks.py. Python 3.6 compatible (Negishi system python).

Every fence is kept, in file order:
  knitr chunk    ```{r label, opts}      lang "r", label from the header
  display fence  ```bash / ```r / ```output / ```{output}
  bare fence     ```                     lang "plain" (the episodes use these for output too)
Indented fences (inside list items) are recognised; the indent is stripped from the body.
Each block records the fenced-div context it sits in (e.g. "challenge/solution"), so the
kit can tell spoiler and solution code from the main line.
"""
import re

CHUNK_RE = re.compile(r'^(\s*)```\{(\w+)\s*([^,}\s]*)\s*,?\s*(.*)\}\s*$')
OPEN_RE = re.compile(r'^(\s*)```\s*(\{?[A-Za-z]*\}?)\s*$')
CLOSE_RE = re.compile(r'^\s*```\s*$')
DIV_OPEN = re.compile(r'^\s*:{3,}\s*\{?\.?\s*([\w-]+)')
DIV_CLOSE = re.compile(r'^\s*:{3,}\s*$')
MDFIG_RE = re.compile(r'!\[([^\]]*)\]\(([^)\s]+)\)(\{[^}]*\})?')
ALT_RE = re.compile(r'alt\s*=\s*"([^"]*)"')
IMG_RE = re.compile(r'<img\s[^>]*src="([^"]+)"[^>]*>')


class Block(object):
    def __init__(self, lang, label, start, end, context, text, fence, opts=""):
        self.lang, self.label, self.start, self.end = lang, label, start, end
        self.context, self.text, self.fence, self.opts = context, text, fence, opts
        self.nn = 0

    def __repr__(self):
        return "Block(%s %s :%d)" % (self.lang, self.label, self.start)


class Episode(object):
    def __init__(self, path, front, chunks, figures):
        self.path, self.front, self.chunks, self.figures = path, front, chunks, figures

    def r_chunks(self):
        return [c for c in self.chunks if c.lang == "r" and c.fence.startswith("```{")]


def parse_rmd(path):
    lines = path.read_text().splitlines()
    front = {}
    i = 0
    if lines and lines[0].strip() == "---":
        i = 1
        while i < len(lines) and lines[i].strip() != "---":
            m = re.match(r'^(\w+):\s*(.*)$', lines[i])
            if m:
                front[m.group(1)] = m.group(2).strip()
            i += 1
        i += 1
    chunks, figs, divs = [], [], []
    while i < len(lines):
        L = lines[i]
        if DIV_CLOSE.match(L):
            if divs:
                divs.pop()
            i += 1
            continue
        if DIV_OPEN.match(L):
            divs.append(DIV_OPEN.match(L).group(1))
            i += 1
            continue
        m = CHUNK_RE.match(L)
        o = OPEN_RE.match(L) if not m else None
        if m or o:
            indent = (m or o).group(1)
            if m:
                lang, label, opts, fence = m.group(2).lower(), m.group(3), m.group(4), L.strip()
            else:
                tag = o.group(2)
                lang = tag.strip("{}").lower() or "plain"
                label, opts, fence = "", "", "```" + tag
            start, body = i + 1, []
            i += 1
            while i < len(lines) and not CLOSE_RE.match(lines[i]):
                body.append(lines[i][len(indent):] if lines[i].startswith(indent) else lines[i])
                i += 1
            chunks.append(Block(lang, label, start, i + 1, "/".join(divs) or "-",
                                "\n".join(body) + "\n" if body else "", fence, opts))
            i += 1
            continue
        for f in MDFIG_RE.finditer(L):
            alt = ALT_RE.search(f.group(3) or "")
            figs.append((i + 1, f.group(2), alt.group(1) if alt else "", f.group(1)))
        for f in IMG_RE.finditer(L):
            alt = ALT_RE.search(f.group(0))
            figs.append((i + 1, f.group(1), alt.group(1) if alt else "", ""))
        i += 1
    n = 0
    for c in chunks:
        n += 1
        c.nn = n
    return Episode(path, front, chunks, figs)
