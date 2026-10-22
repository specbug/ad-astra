# Replacing arXiv:2610.26056 with the journal version (v2)

Submit to AJ first: the comment below says the paper is submitted. The files here are built
(`python build_journal.py` in `pipeline/`) from the final manuscript, which cites data record version 4
(doi:10.5281/zenodo.23553448). Upload them as they are.

## Files

- `tic402710280_arxiv.tar.gz`: `tic402710280.tex`, `tic402710280.bbl`, `aastex701.cls` and `fig1-6.pdf`, the AJ
  source without line numbers. It compiles from an empty directory with `pdflatex` alone, as arXiv does.
- `tic402710280.pdf`: the PDF it produces, 11 pages. The same file goes to Zenodo as the preprint version 3.
- `abstract.txt`: the abstract for the form.

## Steps

Your replacement is already unsubmitted: it sits under **Article Submissions** as `submit/8266785`, "Replacement of
2610.26056", status *incomplete* (it expires on 2026-11-05). Finish that one; do **not** click *Start new submission*
or *Replace* again. The account allows two submissions per calendar month and one is left, so check the preview
carefully before the final click.

1. Submit to AJ first (STEPS.md §2).
2. arxiv.org → your account → in the `submit/8266785` row, click the **pencil icon** (Update).
3. **License / policy pages:** leave everything as it is and continue.
4. **Add Files:** delete every file already there (*Delete all files*, or the bin icon on each row), then upload
   `tic402710280_arxiv.tar.gz` from this folder. Check that the list shows `tic402710280.tex`, `tic402710280.bbl`,
   `aastex701.cls` and `fig1.pdf`–`fig6.pdf`, nothing else. Continue.
5. **Process / Review Files:** the scan reports three missing files, `aasjournalv7`, `CoFirstAuthor` and
   `orcid-ID.png`. They are optional (see below); uncheck their rows. Compiler pdfLaTeX, top-level file
   `tic402710280.tex`. Continue.
6. **Metadata:** title, author and category unchanged. Replace the abstract with `abstract.txt`. Comments:
   `11 pages, 6 figures, 4 tables. Submitted to AJ. Data and code: doi:10.5281/zenodo.23390273. Pre-registered test
   of the 2027 transits: doi:10.5281/zenodo.23390325`
7. **Preview:** open the PDF arXiv built and compare it with `tic402710280.pdf` in this folder. Check:
   - 11 pages, no line numbers;
   - acknowledgments begin "The author thanks Keshav Aggarwal (student visitor, Space Physics Laboratory, Vikram
     Sarabhai Space Centre, India)";
   - §6.4 says "a circle of 2 pixels radius in Sector 93", and its footnote ends with "NFPP = 6 × 10⁻⁷." (no
     sentence about stacking on the linear ephemeris);
   - references read "(Gaia Collaboration et al. 2016, 2023)" with no space after the parenthesis, and no ",,".
8. **Submit.** The status changes to *submitted*; v2 is announced on the next announcement day and v1 stays
   available.

Why the three files are reported missing: `orcid-ID.png` and `CoFirstAuthor.png` are icons that `aastex701.cls`
adds next to author names only if the files exist; `aasjournalv7` is the bibliography style, which arXiv never needs
because the `.bbl` is included.

## Zenodo preprint record (doi:10.5281/zenodo.23390316)

New version → remove the old PDF, upload `tic402710280.pdf` → Version `3.0` → append to the description:

> Version 3 is the manuscript submitted to the Astronomical Journal, identical to arXiv:2610.26056v2. The predicted
> transit times now carry the covariance of T0 and P, the alias test covers every integer fraction of the period
> down to one day, and the acknowledgments thank Keshav Aggarwal and Pascal S. D'Andrea.

Publish.
