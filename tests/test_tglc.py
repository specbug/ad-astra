"""Guard independent quality fields and explicit native time/identity."""
import numpy as np
import pytest
from astropy.io import fits
from adastra.tglc import read_tglc,series,FLUX_COLUMNS

def fixture(tmp_path):
    primary=fits.PrimaryHDU();h=primary.header
    h.update({'ORIGIN':'UCSB/TGLC','TELESCOP':'TESS','TICID':'201308024','GAIADR3':4906011953900670976,'SECTOR':1,'TIMESYS':'TDB'})
    columns=[fits.Column(name='time',format='D',array=np.arange(110)/48+1325)]
    for name in FLUX_COLUMNS+('background',):columns.append(fits.Column(name=name,format='E',array=np.ones(110)))
    for name in ('TESS_flags','TGLC_flags'):columns.append(fits.Column(name=name,format='J',array=np.zeros(110,dtype=int)))
    table=fits.BinTableHDU.from_columns(columns,name='LIGHTCURVE');table.header.update({'GAIADR3':4906011953900670976,'TIMESYS':'TDB','TIMEUNIT':'d','TIMEREF':'SOLARSYSTEM','BJDREFI':2457000,'BJDREFR':0.0})
    p=tmp_path/'native.fits';fits.HDUList([primary,table]).writeto(p);return p

@pytest.mark.parametrize('key,value',[('TIMESYS','UTC'),('BJDREFI',2400000),('TIMEUNIT','s'),('GAIADR3',123)])
def test_native_time_and_identity_guard(tmp_path,key,value):
    p=fixture(tmp_path)
    with fits.open(p,mode='update') as h:h[1].header[key]=value
    with pytest.raises(ValueError):read_tglc(p,201308024,4906011953900670976,1)

def test_two_flags_are_elementwise_and(tmp_path):
    p=fixture(tmp_path);cols,_=read_tglc(p,201308024,4906011953900670976,1)
    cols['TESS_flags'][0]=1;cols['TGLC_flags'][1]=1
    t,y,u=series(cols,'psf_flux')
    assert len(t)==108 and t[0]==cols['time'][2]
    assert u=='relative ppm'

def test_no_later_sector(tmp_path):
    with pytest.raises(ValueError):read_tglc(fixture(tmp_path),201308024,4906011953900670976,27)
