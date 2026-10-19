# Publishing TIC 402710280.01 on Zenodo

## Version 2 (2026 October 19)

Version 2 updates two of the three records: A (data and code) and B (preprint). **Record C, the pre-registration,
is not touched.** Its value comes from being fixed before the 2027 data exist.

What changed and why: QLP light curves of Sectors 100-103 exist at MAST but not on the AWS mirror the pipeline reads,
so version 1 stated wrongly that none existed. Version 2 corrects that statement and uses the QLP light curves to
check our timing. It also refits the Sector 66 window with the final transit shape. The light curves, transit times,
fits, ephemeris and predictions are unchanged; `release/tic402710280_zenodo/SHA256SUMS.txt` shows that every data
file is byte-identical to version 1.

Zenodo keeps version 1 online and citable. The concept DOIs (A: 10.5281/zenodo.23390273, B:
10.5281/zenodo.23390316) always resolve to the latest version; the preprint now cites A by its concept DOI.

For each of A and B, open the published record on Zenodo and press **New version**. That opens a draft with a new DOI
and the old files. Then:

**A, data and code** (opened from https://doi.org/10.5281/zenodo.23390274)
1. Delete the old file and upload `release/tic402710280_zenodo.zip`.
2. Version: `2`. Publication date: the day you publish.
3. Append to the description: *Version 2 adds a comparison with the QLP light curves of Sectors 100 and 102
   (`data/results/qlp_comparison.json`, `qlp_recheck_402710280.py`) and the version 2 manuscript. Data products,
   transit times, fits and predictions are unchanged from version 1.*
4. Related identifiers: change "is supplement to" 10.5281/zenodo.23390317 to 10.5281/zenodo.23390316 (the preprint, all
   versions), keep the rest, and publish.

**B, preprint** (opened from https://doi.org/10.5281/zenodo.23390317)
1. Delete the old PDF and upload `paper/tic402710280.pdf`.
2. Version: `2`. Publication date: the day you publish.
3. Description: paste the version 2 abstract (`release/arxiv/abstract.txt`), then add: *Version 2 corrects a
   statement of version 1 (QLP light curves of Sectors 100-103 do exist at MAST), adds a check of the transit times
   against them, and refits the Sector 66 window. The ephemeris, parameters and predicted transits are unchanged.*
4. Related identifiers: "is supplemented by" 10.5281/zenodo.23390273, "is referenced by" 10.5281/zenodo.23390325.
   Publish.

After publishing, check that each new version's files match this repository byte for byte (`md5sum`). Then submit
the same PDF to arXiv (`release/arxiv/ARXIV.md`).

---

## Version 1 (2026 October 15)

**Status: published 2026-10-15.** A: https://doi.org/10.5281/zenodo.23390274 (data and code), B:
https://doi.org/10.5281/zenodo.23390317 (preprint), C: https://doi.org/10.5281/zenodo.23390325 (pre-registration).
Every published file was checked byte for byte (MD5) against this repository. The steps below are kept as the record
of how they were made; the next action on C is a new version with the outcome of the 2027 transits.

Three records, uploaded in this order so that each can cite the others. Zenodo can **reserve a DOI** for a draft
before it is published ("Reserve DOI" in the upload form); use it to get all three DOIs first.

## Before uploading

1. Author details are in place (rishitv@protonmail.com, ORCID 0009-0003-8778-9883); on Zenodo, add the same ORCID
   to the creator field of all three records.
2. The three DOIs are reserved and already written into every file:
   - A, data and code record: `10.5281/zenodo.23390274`
   - B, preprint: `10.5281/zenodo.23390317`
   - C, pre-registration: `10.5281/zenodo.23390325`

   Upload each file to the draft whose DOI it carries; a mismatch would be permanent.
3. Re-check ExoFOP for TIC 402710280 on the day (target page and the full TOI/CTOI tables).

## A. Data and code record (DOI 10.5281/zenodo.23390274)

- File: `release/tic402710280_zenodo.zip` (upload the zip, or its contents).
- Resource type: Dataset. Licence: Creative Commons Attribution 4.0 (the code inside is MIT; see LICENSE-CODE).
- Title: *TIC 402710280.01: data products and analysis code*
- Description: Light curves of every transit window, TESS full-frame-image photometry of Sectors 100-103, transit
  times, posterior summaries and chains, alias, secondary-eclipse, timing, blend and TRICERATOPS results, the
  characterisation of the 3.78-day signal, the analysis pipeline, and the manuscript source for the preprint "A
  36.34-day Sub-Neptune Candidate Transiting the Nearby K Dwarf TIC 402710280".
- Keywords: exoplanets; transit photometry; TESS; sub-Neptune; TIC 402710280
- Related identifiers: "is supplement to" the preprint DOI (B).

## B. Preprint (DOI 10.5281/zenodo.23390317)

- File: `paper/tic402710280.pdf`.
- Resource type: Publication → Preprint. Licence: CC BY 4.0.
- Title: *A 36.34-day Sub-Neptune Candidate Transiting the Nearby K Dwarf TIC 402710280*
- Description: paste the abstract from the PDF.
- Keywords: as above.
- Related identifiers: "is supplemented by" the data record (A); "is referenced by" the pre-registration (C).

## C. Pre-registration (DOI 10.5281/zenodo.23390325)

- File: `release/prereg_tic402710280.zip` (or its six files).
- Resource type: Other (or Publication → Other). Licence: CC BY 4.0 (code MIT).
- Title: *A Pre-registered Test of the Transiting Planet Candidate TIC 402710280.01 at Its Predicted Transits of 2027*
- Description: This record fixes, before the data exist, the predicted transit times of TIC 402710280.01 for 2027
  (TESS Sector 114 on 2027 February 19; ground-based photometry on 2027 May 3), the data to be used, the processing,
  and the rules that assign each tested transit one of four outcomes. Start with PROTOCOL.pdf.
- Related identifiers: "references" the preprint (B) and the data record (A).
- Must be published before 2027 February 19, and before any Sector 114 data are released.

## After publishing

- arXiv: once you have an endorser for astro-ph.EP, submit the same LaTeX source (`paper/`), citing the Zenodo DOIs.
- ExoFOP: a new community candidate still requires refereed publication; the preprint is the record to cite when
  that paper is written.
