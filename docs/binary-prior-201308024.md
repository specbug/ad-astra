# Scoped binary-catalogue review of TIC 201308024

The target's Gaia radial-velocity and astrometric diagnostics motivate checking prior binary reports; they do not themselves establish a binary or an orbit.

The authoritative CDS descriptions for SB9 (`B/sb9`) and Traven et al. (2020), *The GALAH survey: Multiple stars and our Galaxy* (`J/A+A/638/A145`), were retrieved before choosing query columns. Exact GALAH spectrum `170909003201311` and Gaia DR2 `4906011953900670976` queries each returned the requested Traven binary table with zero rows. A target-centered 10-arcsecond SB9 positional query also returned its requested table with zero rows.

These are scoped nonmatches in two catalogues, not an exhaustive binary-literature clearance. They do not exclude a single-lined, unresolved, unmodeled or previously reported binary, nor establish new variability. The nearby bright star TIC 201308023 already has a Gaia spectroscopic trend, as documented in the stellar-context report.

Initial strict parsing of the GALAH table rejected a legacy VOTable unit `log(cm.s**-2)`. The authoritative ReadMe identifies these columns as logarithmic surface gravity in `[cm/s2]`. A separately declared local reparse accepts that specific metadata warning and generated XML-ID warnings; it performs no physical unit conversion and rejects unexpected warnings or service errors. The shared parser and frozen science pipeline were unchanged. Original failed parses, ReadMes, raw XML, query provenance and the amendment remain in `data/binary-prior-201308024`.

The reproducible helpers are `scripts/review_sap_binary_catalogs.py` and `scripts/parse_sap_sb2_catalogs.py`. No light curves, epoch velocities or protected TESS data were opened.
