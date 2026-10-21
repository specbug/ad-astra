"""TIC 402710280.01: the comparison sample quoted in the paper's introduction, from the NASA Exoplanet Archive.

Confirmed transiting planets (Planetary Systems Composite table, doi:10.26133/NEA13) of 3-4 Earth radii around stars
of 3900-5300 K brighter than V = 10.5, with the archive's own parameter values. The query and the date it was run are
written into the file header, so the table can be regenerated and compared.

  python archive_comparison.py
Output: results/candidates/archive_comparison.csv
"""
import csv
import datetime
import io

import env  # noqa: F401
import requests

from env import RESULTS

TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
QUERY = ("select pl_name, pl_orbper, pl_rade, sy_vmag, st_teff, disc_refname from pscomppars "
         "where tran_flag = 1 and pl_rade >= 3 and pl_rade <= 4 and sy_vmag < 10.5 "
         "and st_teff >= 3900 and st_teff <= 5300 order by pl_orbper")
OUT = RESULTS / "candidates" / "archive_comparison.csv"


def main():
    r = requests.get(TAP, params={"query": QUERY, "format": "csv"}, timeout=120)
    r.raise_for_status()
    rows = list(csv.DictReader(io.StringIO(r.text)))
    with open(OUT, "w", newline="") as fh:
        fh.write(f"# NASA Exoplanet Archive TAP, {datetime.date.today().isoformat()}: {QUERY}\n")
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(len(rows), "planets;", sum(float(x["pl_orbper"]) > 30 for x in rows), "with P > 30 d")


def load():
    """(query date, rows) of the stored comparison table."""
    with open(OUT) as fh:
        head = fh.readline()
        return head.split(":")[0].split(",")[-1].strip(), list(csv.DictReader(fh))


if __name__ == "__main__":
    main()
