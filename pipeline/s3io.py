"""Read TESS light curves straight from the MAST public S3 bucket (stpubdata).

Products:
  spoc   SPOC 2-minute target light curves (PDCSAP), tess/public/tid/sNNNN/...
  tspoc  TESS-SPOC full-frame-image light curves (PDCSAP), mast/hlsp/tess-spoc/...
  qlp    MIT Quick-Look Pipeline full-frame-image light curves, mast/hlsp/qlp/...

Every reader returns the same dict of numpy arrays (time in BTJD, flux normalised to ~1) plus header
metadata, so downstream code never cares which product it got.
"""
import io
import json
import re
import threading
import time

import numpy as np
import requests
from astropy.io import fits
from requests.adapters import HTTPAdapter

from env import CATALOGS

BUCKET = "https://stpubdata.s3.amazonaws.com"
_local = threading.local()

# lightkurve's default TESS quality bitmask: attitude tweak, safe mode, coarse/earth point, desaturation,
# manual exclude, impulsive outlier, argabrightening, straylight-in-aperture flags are dropped.
SPOC_BAD_BITS = 1 | 2 | 4 | 8 | 32 | 128 | 1024 | 2048 | 4096 | 16384


def session():
    s = getattr(_local, "s", None)
    if s is None:
        s = requests.Session()
        s.mount("https://", HTTPAdapter(pool_connections=8, pool_maxsize=64, max_retries=3))
        _local.s = s
    return s


def fetch(key, retries=4):
    """GET one object; None if it does not exist."""
    url = key if key.startswith("http") else f"{BUCKET}/{key}"
    for attempt in range(retries):
        try:
            r = session().get(url, timeout=60)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            return r.content
        except requests.RequestException:
            if attempt == retries - 1:
                raise
            time.sleep(1.5 * 2 ** attempt)


def list_keys(prefix, delimiter=None, max_keys=None):
    """Yield (key, size) under a prefix, or common prefixes when a delimiter is given."""
    token = None
    n = 0
    while True:
        params = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
        if delimiter:
            params["delimiter"] = delimiter
        if token:
            params["continuation-token"] = token
        for attempt in range(5):
            try:
                r = session().get(BUCKET + "/", params=params, timeout=60)
                r.raise_for_status()
                break
            except requests.RequestException:
                if attempt == 4:
                    raise
                time.sleep(2 ** attempt)
        txt = r.text
        if delimiter:
            for p in re.findall(r"<CommonPrefixes><Prefix>(.*?)</Prefix></CommonPrefixes>", txt):
                yield p, 0
        for k, s in re.findall(r"<Key>(.*?)</Key>.*?<Size>(\d+)</Size>", txt):
            yield k, int(s)
            n += 1
            if max_keys and n >= max_keys:
                return
        m = re.search(r"<NextContinuationToken>(.*?)</NextContinuationToken>", txt)
        if not m:
            return
        token = m.group(1)


def tic_dir(tic):
    s = f"{int(tic):016d}"
    return f"{s[0:4]}/{s[4:8]}/{s[8:12]}/{s[12:16]}"


def qlp_key(tic, sector):
    return (f"mast/hlsp/qlp/s{sector:04d}/{tic_dir(tic)}/"
            f"hlsp_qlp_tess_ffi_s{sector:04d}-{int(tic):016d}_tess_v01_llc.fits")


def tspoc_key(tic, sector):
    return (f"mast/hlsp/tess-spoc/s{sector:04d}/target/{tic_dir(tic)}/"
            f"hlsp_tess-spoc_tess_phot_{int(tic):016d}-s{sector:04d}_tess_v1_lc.fits")


_spoc_tmpl = None
_spoc_lock = threading.Lock()


def spoc_templates():
    """Per-sector filename stem of SPOC 2-min light curves, e.g. ('tess2018206045859', '0120')."""
    global _spoc_tmpl
    with _spoc_lock:
        if _spoc_tmpl is None:
            path = CATALOGS / "spoc_templates.json"
            _spoc_tmpl = json.loads(path.read_text()) if path.exists() else {}
        return _spoc_tmpl


def build_spoc_templates(sectors):
    tm = spoc_templates()
    pat = re.compile(r"/(tess\d{13})-s(\d{4})-(\d{16})-(\d{4})-s_lc\.fits$")
    for s in sectors:
        if str(s) in tm:
            continue
        for k, _ in list_keys(f"tess/public/tid/s{s:04d}/", max_keys=4000):
            m = pat.search(k)
            if m and int(m.group(2)) == s:
                tm[str(s)] = [m.group(1), m.group(4)]
                break
    (CATALOGS / "spoc_templates.json").write_text(json.dumps(tm, indent=0))
    return tm


def spoc_key(tic, sector):
    tm = spoc_templates().get(str(sector))
    if tm is None:
        return None
    return (f"tess/public/tid/s{sector:04d}/{tic_dir(tic)}/"
            f"{tm[0]}-s{sector:04d}-{int(tic):016d}-{tm[1]}-s_lc.fits")


def key_for(product, tic, sector):
    return {"qlp": qlp_key, "tspoc": tspoc_key, "spoc": spoc_key}[product](tic, sector)


def _col(d, name):
    return np.asarray(d[name], dtype=float) if name in d.columns.names else None


def parse_lc(raw, product):
    """FITS bytes -> standard light-curve dict (bad cadences removed, flux normalised)."""
    with fits.open(io.BytesIO(raw), memmap=False) as h:
        hd0 = h[0].header
        d = h[1].data
        hd1 = h[1].header
        t = _col(d, "TIME")
        q = np.asarray(d["QUALITY"], dtype=np.int64)
        meta = {"tic": int(hd0.get("TICID", 0)), "sector": int(hd0.get("SECTOR", 0)),
                "camera": int(hd0.get("CAMERA", 0)), "ccd": int(hd0.get("CCD", 0)),
                "tmag": float(hd0.get("TESSMAG") or np.nan), "product": product,
                "cadence_s": float(hd1.get("TIMEDEL", np.nan)) * 86400,
                "crowdsap": float(hd1.get("CROWDSAP") or np.nan),
                "flfrcsap": float(hd1.get("FLFRCSAP") or np.nan)}
        if product == "qlp":
            f = _col(d, "SYS_RM_FLUX")
            if f is None or not np.isfinite(f).any():
                f = _col(d, "SAP_FLUX")
            good = (q == 0)
            extra = {"sap_flux": _col(d, "SAP_FLUX"), "det_flux": _col(d, "DET_FLUX"), "flux_sml": _col(d, "DET_FLUX_SML"),
                     "flux_lag": _col(d, "DET_FLUX_LAG"), "bkg": _col(d, "SAP_BKG"),
                     "cx": _col(d, "SAP_X"), "cy": _col(d, "SAP_Y"), "orbit": _col(d, "ORBITID")}
            ferr = _col(d, "DET_FLUX_ERR")
        else:
            f = _col(d, "PDCSAP_FLUX")
            ferr = _col(d, "PDCSAP_FLUX_ERR")
            good = (q & SPOC_BAD_BITS) == 0
            extra = {"sap_flux": _col(d, "SAP_FLUX"), "bkg": _col(d, "SAP_BKG"),
                     "cx": _col(d, "MOM_CENTR1"), "cy": _col(d, "MOM_CENTR2")}
    good &= np.isfinite(t) & np.isfinite(f) & (f > 0)
    if ferr is None:
        ferr = np.full_like(f, np.nan)
    out = {"time": t[good], "flux": f[good], "flux_err": ferr[good]}
    for k, v in extra.items():
        if v is not None:
            out[k] = v[good]
    med = np.nanmedian(out["flux"]) if good.any() else np.nan
    for k in ("flux", "flux_err", "sap_flux"):
        if k in out and np.isfinite(med):
            out[k] = out[k] / (med if k != "sap_flux" else np.nanmedian(out[k]))
    out["meta"] = meta
    return out


def load(product, tic, sector):
    key = key_for(product, tic, sector)
    if key is None:
        return None
    raw = fetch(key)
    return None if raw is None else parse_lc(raw, product)


def read_target_list(product, sector):
    """TIC IDs in one sector's QLP or TESS-SPOC delivery (from the HLSP target_lists CSVs)."""
    import pandas as pd
    base = {"qlp": "mast/hlsp/qlp", "tspoc": "mast/hlsp/tess-spoc"}[product]
    raw = fetch(f"{base}/target_lists/s{sector:04d}.csv")
    if raw is None:
        return None
    df = pd.read_csv(io.BytesIO(raw), comment=None)
    df.columns = [c.strip("#").strip() for c in df.columns]
    return df.rename(columns={"TIC_ID": "tic", "RA": "ra", "DEC": "dec"})
