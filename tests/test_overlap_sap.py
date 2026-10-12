import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from qlp_overlapfield_sap_followup import common_mode,residualize
from qlp_variability_followup import fit

def test_joint_commonmode_preserves_injected_sinusoid():
    t=np.arange(1340,1360,.02); p=3.2
    nuisance=800*np.sin(2*np.pi*(t-1340)/2.5)+400*np.cos(2*np.pi*(t-1340)/7)
    signal=200*np.sin(2*np.pi*(t-1340)/p)
    y=signal+1.7*nuisance+30*(t-t.mean())
    good,clean,info=residualize(t,y,p,nuisance)
    recovered,_=fit(t[good],clean,p)
    assert abs(recovered['amplitude']-200)<1e-7
    assert abs(info['common_mode_coefficient']-1.7)<1e-9

def test_common_mode_never_interpolates_missing_samples():
    t=np.arange(5.)
    mode,count=common_mode(t,[(t,np.arange(5.)),(t[::2],np.arange(5.)[::2])])
    assert np.isnan(mode[1]) and np.isnan(mode[3])
    assert np.all(count[[0,2,4]]==2)
