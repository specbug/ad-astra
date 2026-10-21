# TIC 402710280.01: what to do, in order

## 1. Zenodo: data record version 4 (`1_zenodo_data_v4/`)

1. Open doi:10.5281/zenodo.23390273 → **New version** → **Get a DOI now**. Send me the reserved DOI; I put it in the
   manuscript and rebuild the zip and the AJ files.
2. Files: delete the version 3 zip; upload the rebuilt `tic402710280_zenodo.zip`.
3. Version: `4.0`. Publication date: today.
4. Description: append the paragraph in `zenodo_v4_data.md`. Related works and the AAS Journals request carry over.
5. Publish.

Do not touch the pre-registration.

## 2. Zenodo: preprint record, metadata only (skip if already done)

Edit the latest version of doi:10.5281/zenodo.23390316 → Related works: add `arXiv:2610.26056`, relation
**Is identical to**. Save.

## 3. AJ submission at https://aas.msubmit.net

Field by field in `3_aj_form/SUBMISSION_DETAILS.md`. In short:

1. Account: link ORCID; enter your full postal address in the profile (not shown in the paper).
2. Journal **AJ**; corridor **The Solar System, Exoplanets, and Astrobiology**.
3. Upload every file in `2_aj_upload_files/` individually: `tic402710280.tex` as the manuscript, `fig1.pdf`–`fig6.pdf`
   as figures, the rest as supporting LaTeX files. Compare the built PDF with `tic402710280_AJ.pdf` (11 pages).
4. Abstract: paste `3_aj_form/abstract_plain.txt`.
5. Dual anonymous review: opt out.
6. Datasets: MAST DOIs listed in the details file. LLM use: yes, paste the disclosure paragraph from the acknowledgments.
7. Financial: request the **AAS Publication Support Fund**, 100% to you; paste `3_aj_form/SUPPORT_FUND.txt`. Answer its
   questions (career stage, funding) truthfully.
8. Cover letter: paste `3_aj_form/COVER_LETTER.txt`. Submit.

## 4. Indexing (any time)

- NASA ADS: create an account and link ORCID.
- ORCID: Works → Add → Search & link → DataCite → add doi:10.48550/arXiv.2610.26056.
- Google Scholar: create a profile; add the paper when it appears (usually 1–2 weeks).

## 5. While under review

- Acknowledgments for Keshav and Pascal go in at revision, once they say how they want to appear.
- Revision tests prepared in parallel: the Sector 100/102 image-clock check and the Sector 102 systematics spread.
- The journal version differs from arXiv v1 in the prediction uncertainties (34 instead of 30 min) and the extended
  alias test. Replace the arXiv version at acceptance, or earlier if you want the preprint to match the submission.
- Leave Zenodo unchanged until acceptance; then update both to the accepted text, add the journal reference
  on arXiv, and submit the community TOI at ExoFOP.
