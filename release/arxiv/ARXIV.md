# Submitting the preprint to arXiv

Submit **version 2** of the preprint, the same PDF that goes to Zenodo as version 2 of record B.

## Files

- `tic402710280_arxiv_v2.tar.gz`: `tic402710280.tex`, `tic402710280.bbl` and `figures/fig1-6.pdf`. It was tested by
  compiling it in an empty directory with `pdflatex` alone, which is what arXiv does (arXiv does not run BibTeX, so
  the `.bbl` is included). The output is 10 pages, and its text is identical to `paper/tic402710280.pdf`.
- `abstract.txt`: the abstract for the arXiv form, 1900 characters (arXiv's limit is 1920).

## Form fields

- **Archive and category:** Astrophysics, primary `astro-ph.EP` (Earth and Planetary Astrophysics). No
  cross-lists are needed.
- **Licence:** CC BY 4.0, the same as the Zenodo records.
- **Title:** A 36.34-day Sub-Neptune Candidate Transiting the Nearby K Dwarf TIC 402710280
- **Authors:** Rishit Vora
- **Abstract:** paste `abstract.txt`.
- **Comments:** 10 pages, 6 figures, 4 tables. Data and code: doi:10.5281/zenodo.23390273. Pre-registered test of
  the 2027 transits: doi:10.5281/zenodo.23390325.

## Endorsement

arXiv asks a first-time submitter to astro-ph for an endorsement. Starting the submission produces a six-character
endorsement code. Send that code to the endorser, who enters it at <https://arxiv.org/auth/endorse>. arXiv itself
checks whether the endorser is eligible for astro-ph, so a refusal there means asking someone else, not resubmitting.
