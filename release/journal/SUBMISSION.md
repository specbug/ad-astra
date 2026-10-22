# Submitting to the Astronomical Journal

AJ and ApJ share one submission portal, one set of author guidelines and one fee schedule. The AAS assigns
observational exoplanet papers to AJ: in the corridor *The Solar System, Exoplanets, and Astrobiology*, "the division
between AJ and ApJ primarily depends on whether the content is primarily observational or theoretical". The Lead
Editor makes the final assignment.

## Files (folders of the submission bundle)

| File | Use |
|---|---|
| `3_aj_upload_files/` | Flat LaTeX source: `tic402710280.tex`, `.bbl`, `refs.bib`, `aastex701.cls`, `aasjournalv7.bst`, `fig1.pdf`–`fig6.pdf`. Compiles with pdflatex alone. The portal cannot read folders, so upload the files individually. |
| `tic402710280_AJ.pdf` (top level) | The compiled manuscript (AASTeX 7.0.1, line numbers on), for checking against the PDF the portal builds. |
| `4_aj_form/abstract_plain.txt` | The abstract in plain Unicode for the form (226 words; the limit is 250). |
| `4_aj_form/COVER_LETTER.txt` | Paste into the cover-letter field. |
| `4_aj_form/SUPPORT_FUND.txt` | Justification for the AAS Publication Support Fund request. |

The data record version 4 (10.5281/zenodo.23553448), which the manuscript cites, is published.

## The portal: https://aas.msubmit.net

1. Create an account, link your ORCID, and enter your postal address.
2. **Submit a manuscript** → journal **AJ**, corridor **The Solar System, Exoplanets, and Astrobiology**.
3. **Files:** upload `tic402710280.tex` as the manuscript; `fig1.pdf`–`fig6.pdf` as figures; `refs.bib`,
   `tic402710280.bbl`, `aastex701.cls`, `aasjournalv7.bst` as supporting LaTeX files. Check the PDF the portal builds
   against `tic402710280_AJ.pdf`: 11 pages, six figures, four tables, line numbers.
4. **Title, abstract, authors:** the title as in the manuscript; paste `abstract_plain.txt`; one author, Rishit Vora,
   corresponding.
5. **Dual anonymous review:** opt out. The preprint is already public on arXiv and Zenodo.
6. **Ethics:** read and accept the Professional and Ethical Standards.
7. **Digital assets:** none in the article. Give the Zenodo data record (version 4, 10.5281/zenodo.23553448) as data stored in a
   repository.
8. **Datasets in astronomy archives:** yes, MAST.
   - TESS-SPOC light curves: 10.17909/t9-wpz1-8s54
   - QLP light curves: 10.17909/t9-r086-e880
   - TESS calibrated full-frame images: 10.17909/0cp4-2j79
   - TESS two-minute light curves (Sectors 99–107, used in the search): 10.17909/t9-nmc8-f686
   - TESS Input Catalog: 10.17909/fwdt-2x66
9. **LLM use:** answer yes, and paste the disclosure paragraph from the acknowledgments.
10. **Financial:** request support from the **AAS Publication Support Fund**, with 100% of charges allocated to you.
    Answer the questions it asks (career stage, grant funding) truthfully, and paste `SUPPORT_FUND.txt` where it
    asks for a justification. At 2026 rates the article falls in Tier 1 (≤ 30 quanta: about 12 for text,
    6 figures, 4 tables), $1,425 before support.
11. **Funding:** none.
12. **UAT concepts:** Exoplanets (498); Transit photometry (1709); Mini Neptunes (1063); K dwarf stars (876).
13. **Cover letter:** paste `COVER_LETTER.txt`.
14. Submit. The charge confirmation arrives within a few business days.

## After submission

- Expect the referee report within several weeks. Revisions use the `trackchanges` class option and
  `\added{}`/`\deleted{}` markup, with a response letter.
- On acceptance: replace the arXiv version with the accepted text (it becomes v3) and add the journal reference and
  DOI; publish the accepted text as a new version of the Zenodo preprint record (version 4) and of the data record
  (version 5); then submit the community TOI at ExoFOP, which requires the refereed paper.
