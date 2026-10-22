# Replacing arXiv:2610.26056 with the journal version (v2)

The files here are built (`python build_journal.py` in `pipeline/`) from the final manuscript, which cites data record
version 4 (doi:10.5281/zenodo.23553448). Upload them as they are.

## Files

- `tic402710280_arxiv.tar.gz`: `tic402710280.tex`, `tic402710280.bbl`, `aastex701.cls` and `fig1-6.pdf`, the AJ
  source without line numbers. It compiles from an empty directory with `pdflatex` alone, as arXiv does.
- `tic402710280.pdf`: the PDF it produces, 11 pages. The same file goes to Zenodo as the preprint version 3.
- `abstract.txt`: the abstract for the form.

## Steps

1. Log in at arxiv.org → your article 2610.26056 → **Replace**.
2. Upload `tic402710280_arxiv.tar.gz`; check that the processed PDF matches `tic402710280.pdf`.
3. Abstract: paste `abstract.txt`. Title, author, category and licence stay as they are.
4. Comments: `11 pages, 6 figures, 4 tables. Submitted to AJ. Data and code: doi:10.5281/zenodo.23390273.
   Pre-registered test of the 2027 transits: doi:10.5281/zenodo.23390325`
5. Submit. The replacement is announced as v2 the next working day; v1 stays available.

## Zenodo preprint record (doi:10.5281/zenodo.23390316)

New version → remove the old PDF, upload `tic402710280.pdf` → Version `3.0` → append to the description:

> Version 3 is the manuscript submitted to the Astronomical Journal, identical to arXiv:2610.26056v2. The predicted
> transit times now carry the covariance of T0 and P, the alias test covers every integer fraction of the period
> down to one day, and the acknowledgments thank Keshav Aggarwal and Pascal S. D'Andrea.

Publish.
