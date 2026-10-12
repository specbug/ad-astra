"""Strict native TGLC-v1 reader; flux units intentionally not inferred from BUNIT."""
import numpy as np
from adastra.archive import validated_fits

FLUX_COLUMNS = ('psf_flux', 'aperture_flux', 'cal_psf_flux', 'cal_aper_flux')

def read_tglc(path, tic, gaia, sector):
    if sector not in (1, 2):
        raise ValueError('Only declared sectors 1 and 2 are allowed')
    with validated_fits(path) as hdus:
        if len(hdus) != 2 or hdus[1].name != 'LIGHTCURVE':
            raise ValueError('Unexpected native TGLC table layout')
        primary, header = hdus[0].header, hdus[1].header
        if primary.get('ORIGIN') != 'UCSB/TGLC' or primary.get('TELESCOP') != 'TESS':
            raise ValueError('Unexpected TGLC origin')
        if int(primary.get('TICID', -1)) != tic or primary.get('SECTOR') != sector:
            raise ValueError('TIC/sector identity mismatch')
        for h in (primary, header):
            if int(h.get('GAIADR3', -1)) != gaia or h.get('TIMESYS') != 'TDB':
                raise ValueError('Gaia identity or time scale mismatch')
        if header.get('TIMEUNIT') != 'd' or header.get('TIMEREF') != 'SOLARSYSTEM':
            raise ValueError('Missing explicit barycentric day time convention')
        if 'BJDREFI' not in header or 'BJDREFR' not in header:
            raise ValueError('Missing BJD reference')
        reference = float(header['BJDREFI']) + float(header['BJDREFR'])
        if reference != 2457000:
            raise ValueError('Unexpected declared native-v1 BJD reference')
        names = hdus[1].data.names
        required = ('time', 'TESS_flags', 'TGLC_flags', 'background') + FLUX_COLUMNS
        if any(name not in names for name in required):
            raise ValueError('Missing native TGLC columns')
        columns = {name: np.array(hdus[1].data[name], copy=True) for name in names}
        if not all(np.issubdtype(columns[name].dtype, np.integer) for name in ('TESS_flags', 'TGLC_flags')):
            raise ValueError('Quality flags must be integer bit fields')
        metadata = {k: primary.get(k) for k in ('SECTOR', 'CAMERA', 'CCD', 'GAIADR3', 'TICID')}
        metadata.update({k: header.get(k) for k in ('WOTAN_WL', 'WOTAN_MT', 'PSF_ERR', 'APER_ERR', 'CPSF_ERR', 'CAPE_ERR', 'LOC_BG')})
        metadata['time_convention'] = 'Native time = BJD_TDB - 2457000, days'
        metadata['calibrated_unit_caution'] = 'v1 headers describe calibrated scatter as e-/s, but calibrated curves are normalized/detrended; no cadence uncertainties borrowed'
        metadata['checksum_cards'] = [{'hdu':i,'CHECKSUM':'CHECKSUM' in h.header,'DATASUM':'DATASUM' in h.header} for i,h in enumerate(hdus)]
    return columns, metadata

def series(columns, name):
    t, y = np.asarray(columns['time'], float), np.asarray(columns[name], float)
    good = (columns['TESS_flags'] == 0) & (columns['TGLC_flags'] == 0) & np.isfinite(t) & np.isfinite(y)
    # Do not clip or discard negative fluxes; failed normalization is explicit.
    t, y = t[good], y[good]
    order = np.argsort(t, kind='stable'); t, y = t[order], y[order]
    unique = np.r_[True, np.diff(t) > 0]; t, y = t[unique], y[unique]
    if len(t) < 100:
        raise ValueError('Too few usable native TGLC samples')
    if name in FLUX_COLUMNS:
        median = float(np.median(y))
        if not np.isfinite(median) or median <= 0:
            raise ValueError('Nonpositive native flux normalization')
        y = (y/median-1)*1e6; unit = 'relative ppm'
    else:
        y = y - np.median(y); unit = 'native background e-/s'
    return t, y, unit
