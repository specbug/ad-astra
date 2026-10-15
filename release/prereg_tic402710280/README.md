# Pre-registered test of TIC 402710280.01 (2027 transits)

Start with `PROTOCOL.pdf` (two pages). It fixes, before the data exist, the predicted transit times of TIC 402710280.01
for 2027 (TESS Sector 114 on 2027 February 19; ground-based on 2027 May 3), the data to be used, the processing, and
the rules that assign each tested transit one of four outcomes: recovered, inconsistent, excluded, insufficient data.

| File | Content |
|---|---|
| `PROTOCOL.pdf`, `PROTOCOL.tex` | The protocol. |
| `prereg_test.py` | Frozen analysis code (Python 3, numpy, batman). Usage: `python prereg_test.py lightcurve.csv EPOCH [exposure_s]`. |
| `validation.json` | Output of the frozen code on the archival transits and on transit-free data (Section 7 of the protocol). |
| `predictions.json` | The registered ephemeris and predicted mid-times with uncertainties for epochs 57-66. |
| `SHA256SUMS.txt` | Checksums of every file in this record. |

Companion records: the preprint *A 36.34-day Sub-Neptune Candidate Transiting the Nearby K Dwarf TIC 402710280, Identified in TESS Full-Frame Images* (doi:10.5281/zenodo.23390317) and its data/code record (doi:10.5281/zenodo.23390274). This record: doi:10.5281/zenodo.23390325.

Licences: text CC BY 4.0; code MIT.
