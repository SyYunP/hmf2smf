import numpy as np
from scipy.integrate import quad

def stellarmass_func(hmf,Omega_b = 0.0412,Omega_m=0.2725):
    '''
    Sets up a steller mass function from a halo mass function including supernova feedback.
    
    Parameters
    ----------
    hmf : array-like object
        Halo mass function
    Omega_b : `float`, optional
        Cosmic baryon mass density (default value: 0.0412)
    Omega_m : `float`, optional
        Cosmic mass density (default value: 0.2725)
    sfe : `float`, optional
        Star formation efficiency (default value: 1.0)

    Returns
    -------
    smf : `array-like object`, as the input halo mass function
        Stellar mass function
    '''
    sfe = SFE_SNFeedback()
    smf = hmf*(Omega_b/Omega_m)*sfe
    return smf

def snrate(imf,SNthreshold=8,mlowerbound=0.1,mupperbound=100):
    '''
    Solve the supernova rate.

    Parameters
    ----------
    imf : `str` or `function`
        Initial mass function. Users can directly request for `salpeter` (Salpeter 1955) or 
        `kroupa` (Kroupa 2002) initial mass function with `str` or input custom function for 
        integration.
    SNthreshold : `float`, optional
        The lower mass bound of supernova progenitors (default value: 8)
    mlowerbound : `float`, optional
        The lower mass bound of initial mass function (default value: 0.1)
    mupperbound : `float`, optional
        The upper mass bound of initial mass function (default value: 100)

    Returns
    -------
    snrate : `float`
        Supernova rate (number of supernova per solar mass)
    '''
    if isinstance(imf,str):
        imf = _setup_IMF_fromname(imf)
        pass
    snrate = quad(imf,SNthreshold,mupperbound)[0]/quad(_mass_integral,mlowerbound,mupperbound,args=imf)[0]
    return snrate

def _setup_IMF_fromname(imf):
    '''
    Setup common initial mass function from name.

    Supporting models: Salpeter (1955), Kroupa (2002).

    Parameters
    ----------
    imf : `str`
        Name of the initial mass function model. 
    
    Returns
    ----------
    Function of initial mass function.
    '''
    if imf.lower() == 'salpeter':
        return lambda m: m**-2.35
    elif imf.lower() == 'kroupa':
        def Kroupa(m):
            if m > 0.5:
                return m**-2.3
            if m < 0.08:
                return m**-0.3
            else:
                return m**-1.3
        return Kroupa
    else:
        return "This initial mass function does not exist in the package."

def _mass_integral(m,func):
    '''
    Set up for mass integration for any initial mass function.
    
    Parameters
    ----------
    m : array-like object
        mass
    func : `function`
        initial mass function that takes m as input.
    '''
    return m*func(m)

def SFE_SNFeedback(vc,f_gas,SNrate,SNEnergy):
    '''
    Solve star formation efficiency when gas binding energy equals to supernovae feedback energy.

    Parameters
    ----------
    vc : `float`
        Circular velocity of halo
    f_gas : `float`
        Fraction of energy contribution to ambient gas
    SNrate : `float`
        Supernova rate (number/solar mass)
    SNEnergy : `float`
        Released energy from every single supernova

    Returns
    ----------
    sfe : `float`
        Star formation efficiency when binding energy equals to supernova energy
    '''
    sfe = vc*vc/(vc*vc+f_gas*SNrate*SNEnergy)
    return sfe