# Zenodo data record, version 4 (record A only): published as 10.5281/zenodo.23553448

Kept for the record. The only remaining step is the description note at the end of this file.

Upload `release/tic402710280_zenodo.zip` as a new version of doi:10.5281/zenodo.23390273. This step changes only the
data record; the pre-registration is never changed.

1. Open the draft of version 4 (reserved DOI 10.5281/zenodo.23553448).
2. **Files:** remove the version 3 zip; upload `tic402710280_zenodo.zip` (md5 2d64942bb2af548adef4e9c786248d22).
3. **Version:** 4.0. **Publication date:** the day you publish.
4. **Description:** append this paragraph.

> Version 4 accompanies the manuscript submitted to the Astronomical Journal. The uncertainties of the predicted
> transit times now include the covariance of T0 and P, which are correlated (0.47) because the reference epoch is not
> the weighted mean epoch; at the Sector 114 transit the uncertainty is 34 rather than 30 minutes. The pre-registered
> test (doi:10.5281/zenodo.23390325) is unchanged. The alias test now covers every integer fraction of the period
> down to one day. This version adds the Exoplanet Archive comparison sample quoted in the introduction and the
> intermediate results, so that paper_numbers.py regenerates every number in the manuscript without network access.
> The light curves, transit times, fits and ephemeris are unchanged.

5. Related works and communities carry over from version 3; check that `arXiv:2610.26056` (*Is supplement to*) and the
   AAS Journals community request are still there. Publish.

## Note to add to the published version 4 (Edit → Description; metadata only, no new version)

> Note added 2026 October 22. The manuscript in this version (manuscript/ and code/paper/) is a snapshot taken
> shortly before submission. The text submitted to the Astronomical Journal, also posted as arXiv:2610.26056v2,
> differs in wording only: it gives the Sector 93 aperture used in the TRICERATOPS calculation; it removes a
> sentence that gave FPP = 0.015 for stacking on the linear ephemeris, a value that came from an earlier
> configuration of the calculation; it counts 27 rather than 30 additional alias windows at P/4 and P/5; it calls
> periods below one day strongly disfavored rather than excluded; and it gives the acknowledgments in full. The
> data, code and numerical results in this version are unchanged.
