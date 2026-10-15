# Pre-registered test of TIC 402710280.01 (2027 transits)

Start with `PROTOCOL.pdf` (two pages). It fixes, before the data exist, the predicted transit times of TIC 402710280.01
for 2027 (TESS Sector 114 on 2027 February 19; ground-based on 2027 May 3), the data to be used, the processing, and
the rules that assign each tested transit one of four outcomes: recovered, inconsistent, excluded, insufficient data.

| File | Content |
|---|---|
| `PROTOCOL.pdf`, `PROTOCOL.tex` | The protocol. |
| `prereg_test.py` | Frozen analysis code (Python 3, numpy, batman). Usage: `python prereg_test.py lightcurve.csv EPOCH [exposure_s]`. |
| `predictions.json` | The registered ephemeris and predicted mid-times with uncertainties for epochs 57-66. |
| `SHA256SUMS.txt` | Checksums of every file in this record. |

Companion records: the preprint and its data/code record (R. Vora 2026, Zenodo).

Licences: text CC BY 4.0; code MIT.
