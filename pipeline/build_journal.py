"""TIC 402710280.01: the submission files built from the filled manuscript (paper/tic402710280.tex).

  release/journal/tic402710280_AJ_source.zip   flat AASTeX source for the AAS portal (line numbers on)
  release/journal/tic402710280_AJ.pdf          the PDF that source compiles to
  release/arxiv/tic402710280_arxiv.tar.gz      the same source without line numbers, for arXiv
  release/arxiv/tic402710280.pdf               the PDF arXiv will produce; also the Zenodo preprint file
  release/arxiv/abstract.txt                   the abstract for the arXiv form

Each package is compiled from an empty directory with pdflatex alone, as the AAS portal and arXiv do.

  python build_journal.py
"""
import re
import shutil
import subprocess
import tarfile
import tempfile
import zipfile
from pathlib import Path

import env  # noqa: F401

from env import RESULTS

ROOT = RESULTS.parent
PAPER = ROOT / "paper"
FILES = ["tic402710280.tex", "tic402710280.bbl", "refs.bib", "aastex701.cls", "aasjournalv7.bst"]
FIGS = sorted((PAPER / "figures").glob("fig?.pdf"))


def compile_in(d):
    for _ in range(2):
        subprocess.run(["pdflatex", "-interaction=nonstopmode", "tic402710280.tex"], cwd=d, capture_output=True)
    log = (d / "tic402710280.log").read_text(errors="replace")
    assert not re.search(r"^!", log, re.M) and "undefined" not in log, f"LaTeX errors in {d}"
    return d / "tic402710280.pdf"


def stage(d, tex):
    for f in FILES[1:]:
        shutil.copy(PAPER / f, d / f)
    for f in FIGS:
        shutil.copy(f, d / f.name)
    (d / "tic402710280.tex").write_text(tex)


def arxiv_abstract(tex):
    a = tex[tex.find("\\begin{abstract}") + 16:tex.find("\\end{abstract}")].strip()
    for x, y in (("\\tic", "TIC 402710280"), ("\\tess\\ ", "TESS "), ("\\tess", "TESS"), ("\\gaia\\ ", "Gaia "),
                 ("\\gaia", "Gaia"), ("\\texttt{TRICERATOPS}", "TRICERATOPS"), ("\\,\\rearth", "\\,R_\\oplus"),
                 ("~", " ")):
        a = a.replace(x, y)
    left = set(re.findall(r"\\[A-Za-z]+", a)) - {"\\pm", "\\geq", "\\Delta", "\\chi", "\\oplus", "\\sigma"}
    assert not left, f"unconverted macros in the abstract: {left}"
    assert len(a) <= 1920, "arXiv abstract limit"
    return a + "\n"


def main():
    tex = (PAPER / "tic402710280.tex").read_text()
    assert tex.startswith("\\documentclass[twocolumn,linenumbers]{aastex701}")
    jdir, adir = ROOT / "release" / "journal", ROOT / "release" / "arxiv"
    adir.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as t:
        d = Path(t)
        stage(d, tex)
        with zipfile.ZipFile(jdir / "tic402710280_AJ_source.zip", "w", zipfile.ZIP_DEFLATED) as z:
            for f in FILES + [f.name for f in FIGS]:
                z.write(d / f, f)
        shutil.copy(compile_in(d), jdir / "tic402710280_AJ.pdf")
    arx = tex.replace("[twocolumn,linenumbers]", "[twocolumn]", 1)
    with tempfile.TemporaryDirectory() as t:
        d = Path(t)
        stage(d, arx)
        with tarfile.open(adir / "tic402710280_arxiv.tar.gz", "w:gz") as tz:
            for f in ["tic402710280.tex", "tic402710280.bbl", "aastex701.cls"] + [f.name for f in FIGS]:
                tz.add(d / f, f)
        shutil.copy(compile_in(d), adir / "tic402710280.pdf")
    (adir / "abstract.txt").write_text(arxiv_abstract(tex))
    print("AJ source, AJ PDF, arXiv tarball, arXiv PDF and abstract written")


if __name__ == "__main__":
    main()
