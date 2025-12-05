from hmf2smf import *
import numpy as np
import colossus.cosmology.cosmology as c
import pytest

test_HMF = np.logspace(7,15,3)
test_cosmo = c.setCosmology('planck18')

def test_snrate():
    assert snrate('SALPETER') == pytest.approx([0.007421779932978433],rel=1e-8)
    assert snrate('Kroupa') == pytest.approx([0.012845638907782389])
    with pytest.raises(TypeError):
        snrate()
    with pytest.raises(ValueError):
        snrate('topheavy')

def test_SFE_SNFeedback():
    with pytest.raises(TypeError):
        SFE_SNFeedback()
    assert SFE_SNFeedback(7,test_HMF,test_cosmo,0.5) == pytest.approx(np.array([3.74855551e-06, 1.73690976e-03, 4.46781739e-01]),rel = 1e-8)
    assert SFE_SNFeedback(7,test_HMF,test_cosmo,snrate('salpeter')) == pytest.approx(np.array([2.52474688e-04, 1.04919512e-01, 9.81951986e-01]),rel=1e-8)
    assert SFE_SNFeedback(7,1e10,test_cosmo,snrate('salpeter')) == pytest.approx([0.024631797131135728],rel=1e-8)
    assert SFE_SNFeedback(7,1e10,test_cosmo,snrate('salpeter'),0.5,1e7) == pytest.approx([0.20258833460543701],rel=1e-8)

def test_stellarmass_func():
    test_sfe_2 = [50,10]
    with pytest.raises(ValueError):
        stellarmass_func(test_HMF,test_sfe_2,test_cosmo)
    test_sfe_3 = [10,20,30]
    assert stellarmass_func(test_HMF,test_sfe_3,test_cosmo) == pytest.approx(np.array([1.57505625e+07, 3.15011250e+11, 4.72516876e+15]),rel=1e-8)
    assert stellarmass_func(test_HMF5,test_cosmo) == pytest.approx(np.array([7.87528126e+06, 7.87528126e+10, 7.87528126e+14]),rel=1e-8)
    assert stellarmass_func(1e7,test_sfe_3,test_cosmo) == pytest.approx(np.array([15750562.52009001, 31501125.04018001, 47251687.56027002]),rel=1e-8)
    