# TIC 402710280.01: what to do, in order

## 1. Zenodo data record version 4: published; add one note

Published as 10.5281/zenodo.23553448 (md5 2d64942bb2af548adef4e9c786248d22). Its manuscript copy is a
pre-submission snapshot, so open the record → **Edit** → append the note at the end of
`1_zenodo_data_v4/zenodo_v4_data.md` to the description → **Publish**. This edits metadata only and creates no new
version; do not upload files. The pre-registration (10.5281/zenodo.23390325) is never changed.

## 2. AJ submission at https://aas.msubmit.net

Do this before step 3: the arXiv comment and the Zenodo preprint description say "submitted to AJ".
Field by field in `4_aj_form/SUBMISSION_DETAILS.md`. In short:

1. Account: link ORCID; enter your full postal address in the profile (not shown in the paper).
2. Journal **AJ**; corridor **The Solar System, Exoplanets, and Astrobiology**.
3. Upload every file in `3_aj_upload_files/` and `tic402710280_AJ.pdf` (the compiled merged PDF). File Types: the
   `.tex` and `tic402710280_AJ.pdf` as Manuscript file, `fig1.pdf`–`fig6.pdf` as figures, the rest as supporting
   LaTeX files. Compare the portal's own merged PDF with `tic402710280_AJ.pdf` (11 pages).
4. Abstract: paste `4_aj_form/abstract_plain.txt`.
5. Dual anonymous review: opt out.
6. Datasets: MAST yes (five DOIs) and IPAC yes (10.26133/NEA13, 10.26134/ExoFOP5), as listed in the details file. LLM use: yes, paste the disclosure paragraph from the acknowledgments.
7. Financial: request the **AAS Publication Support Fund** for USD 925; you pay USD 500. Paste
   `4_aj_form/SUPPORT_FUND.txt`. Answer its questions (career stage, funding) truthfully.
8. Cover letter: paste `4_aj_form/COVER_LETTER.txt`. Submit.

## 3. arXiv v2 and Zenodo preprint version 3 (`2_arxiv_and_preprint/`)

The replacement already sent to arXiv carries the earlier files. If it is not yet announced, unsubmit it, replace the
files with `tic402710280_arxiv.tar.gz` from this bundle, and resubmit. Then publish `tic402710280.pdf` from the same
folder as version 3 of the Zenodo preprint record. Steps in `ARXIV.md`.

## 4. Indexing (any time)

- NASA ADS: create an account and link ORCID.
- ORCID: Works → Add → Search & link → DataCite → add doi:10.48550/arXiv.2610.26056.
- Google Scholar: create a profile; add the paper when it appears (usually 1–2 weeks).

## 5. While under review

- Revision tests prepared in parallel: the Sector 100/102 image-clock check and the Sector 102 systematics spread.
- Leave arXiv and Zenodo unchanged until acceptance; then update both to the accepted text, add the journal reference
  on arXiv, and submit the community TOI at ExoFOP.
