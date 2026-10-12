import numpy as np
import pytest
from astropy.io import fits

from adastra.qlp import read_qlp


def make_file(tmp_path, *, sector=1, fraction_key="BJDREFR", fraction=0.25,
              conflict=None, timesys="TDB", omit_error=False):
    primary = fits.PrimaryHDU()
    primary.header["SECTOR"] = sector
    primary.header["TICID"] = 261136679
    cols = [fits.Column(name="TIME", array=[100., 101., 102.], format="D"),
            fits.Column(name="SAP_FLUX", array=[1., .999, 1.], format="D"),
            fits.Column(name="KSPSAP_FLUX", array=[1., .999, 1.], format="D"),
            fits.Column(name="QUALITY", array=[0, 0, 1], format="J")]
    if not omit_error:
        cols.append(fits.Column(name="KSPSAP_FLUX_ERR", array=[.0001]*3, format="D"))
    table = fits.BinTableHDU.from_columns(cols)
    table.header["TIMESYS"] = timesys
    table.header["TIMEUNIT"] = "d"
    table.header["BJDREFI"] = 2457000
    table.header[fraction_key] = fraction
    if conflict is not None:
        table.header["BJDREFF"] = conflict
    path = tmp_path / "qlp.fits"
    fits.HDUList([primary, table]).writeto(path)
    return path


@pytest.mark.parametrize("fraction_key", ["BJDREFR", "BJDREFF"])
def test_explicit_time_reference_and_quality_preserved(tmp_path, fraction_key):
    lc = read_qlp(make_file(tmp_path, fraction_key=fraction_key))
    np.testing.assert_allclose(lc.time_bjd, [2457100.25, 2457101.25, 2457102.25], rtol=0, atol=1e-8)
    np.testing.assert_array_equal(lc.quality, [0, 0, 1])
    np.testing.assert_array_equal(lc.good_mask(), [True, True, False])
    assert len(lc.flux) == 3


def test_conflicting_fraction_rejected(tmp_path):
    with pytest.raises(ValueError, match="Conflicting BJDREF"):
        read_qlp(make_file(tmp_path, conflict=.5))


def test_sap_does_not_borrow_corrected_flux_errors(tmp_path):
    lc = read_qlp(make_file(tmp_path), "SAP_FLUX")
    assert lc.flux_error is None
    assert lc.metadata["error_column"] is None


@pytest.mark.parametrize("sector", [0, 74, 79, 94, 99, 104])
def test_out_of_scope_release_rejected(tmp_path, sector):
    with pytest.raises(ValueError, match="outside the validated legacy scope"):
        read_qlp(make_file(tmp_path, sector=sector))


def test_missing_corrected_errors_rejected(tmp_path):
    with pytest.raises(ValueError, match="Missing QLP columns"):
        read_qlp(make_file(tmp_path, omit_error=True))


def test_wrong_time_scale_rejected(tmp_path):
    with pytest.raises(ValueError, match="TIMESYS=TDB"):
        read_qlp(make_file(tmp_path, timesys="UTC"))
