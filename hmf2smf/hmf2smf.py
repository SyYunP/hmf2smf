import numpy as np
from scipy.integrate import quad

__all__ = ["stellarmass_func", "snrate", "SFE_SNFeedback"]

def stellarmass_func(hmf,sfe,cosmo):
    '''
    Sets up a steller mass function from a halo mass function including supernova feedback.
    
    Parameters
    ----------
    hmf : `float` or array-like object
        Halo mass function
    sfe : `float` or array-like object
        Star formation efficiency
        Note that the shape of the sfe must be the same as the hmf
    cosmo : `colossus.cosmology.cosmology`
        Cosmology object from colossus. 
    
    Returns
    -------
    smf : `array-like object`, as the input halo mass function
        Stellar mass function
    '''
    if isinstance(sfe, (list, tuple, np.ndarray)) and isinstance(hmf, (list, tuple, np.ndarray)) and (np.array(sfe).shape != np.array(hmf).shape):
        raise ValueError("The shape of star formation efficiency must be the same as the halo mass function.")
    
    smf = np.array(hmf)*(cosmo.Ob0/cosmo.Om0)*np.array(sfe)

    return smf

def snrate(imf,SNthreshold=8,mlowerbound=0.1,mupperbound=100):
    '''
    Solve the supernova rate.

    Parameters
    ----------
    imf : `str` or `func`
        Initial mass function. Users can directly request for "salpeter" (Salpeter 1955) or 
        "kroupa" (Kroupa 2002) initial mass function with `str` or input custom function for 
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
        return ValueError("This initial mass function does not exist in the package.")

def _mass_integral(m,func):
    '''
    Set up for mass integration for any initial mass function.
    
    Parameters
    ----------
    m : array-like object
        mass
    func : `func`
        initial mass function that takes m as input.
    '''
    return m*func(m)

def SFE_SNFeedback(z,hmf,cosmo,SNrate,f_gas = 1.0,SNEnergy = 50300737):
    '''
    Solve star formation efficiency when gas binding energy equals to supernovae feedback energy.

    Parameters
    ----------
    z : `float`
        Redshift
    hmf : `float` or array-like object
        Halo mass function
    cosmo : `colossus.cosmology.cosmology`
        Cosmology object from colossus. 
    SNrate : `float`
        Supernova rate (number/solar mass)
    f_gas : `float`, optional
        Fraction of energy contribution to ambient gas (default value: 1.0)
    SNEnergy : `float`, optional
        Released energy from every single supernova (solar mass*km^2/s^2 default value: 50300737)

    Returns
    ----------
    sfe : `float`
        Star formation efficiency when binding energy equals to supernova energy
    '''
    vc = _vvir(hmf,z,cosmo)
    vcsq = vc*vc
    sfe = vcsq/(vcsq+f_gas*SNrate*SNEnergy)
    return sfe

def _vir_overdens(z,cosmo):
    return 18*np.pi*np.pi+82*(cosmo.Om(z)-1)-39*(cosmo.Om(z)-1)*(cosmo.Om(z)-1)

def _vvir(hmf,z,cosmo):
    return 23.4*(hmf/1e8)**(1/3)*(cosmo.Om(z)*_vir_overdens(z,cosmo)/cosmo.Om(z)/18/np.pi/np.pi)**(1/6)*((1+z)/10)**0.5