"""Native timing identity/reference failures must be explicit."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest
from astropy.io import fits
spec=importlib.util.spec_from_file_location('native_spoc',Path(__file__).resolve().parents[1]/'scripts/native_spoc_201308023.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

def hdus(unit='BJD - 2457000, days'):
    p=fits.PrimaryHDU();p.header['TICID']=201308023;p.header['SECTOR']=1
    t=fits.BinTableHDU.from_columns([fits.Column(name='TIME',format='D',unit=unit,array=np.array([1325.,1326.]))])
    for k,v in {'TIMESYS':'TDB','TIMEUNIT':'d','BJDREFI':2457000,'BJDREFF':0.}.items():t.header[k]=v
    return fits.HDUList([p,t])

def test_native_bjd_unit_reference():
    np.testing.assert_array_equal(module.native_time(hdus(),201308023,1),[2458325.,2458326.])

def test_conflicting_embedded_reference_rejected():
    with pytest.raises(ValueError,match='reference conflict'):module.native_time(hdus('BJD - 2400000, days'),201308023,1)

@pytest.mark.parametrize('key,value',[('TIMESYS','UTC'),('TIMEUNIT','s')])
def test_wrong_native_time_system_rejected(key,value):
    h=hdus();h[1].header[key]=value
    with pytest.raises(ValueError):module.native_time(h,201308023,1)

def test_wrong_target_or_sector_rejected():
    with pytest.raises(ValueError,match='identity'):module.native_time(hdus(),201308024,1)
    with pytest.raises(ValueError,match='identity'):module.native_time(hdus(),201308023,68)

def test_primary_extension_reference_disagreement_rejected():
    h=hdus();h[0].header['BJDREFI']=2400000
    with pytest.raises(ValueError,match='Conflicting'):module.native_time(h,201308023,1)
