hmf2smf Documentation
---------------------

`hmf2smf` calculates stellar mass function from halo mass function. It relies on 
`NumPy <https://numpy.org/>`_, `SciPy <https://scipy.org/>`_, and 
`Colossus <https://bdiemer.bitbucket.io/colossus/index.html>`_. The package is 
still under construction so all comments, suggestions or contribution are welcomed! 
Please open a new issue for bugs, feedback or feature requests on the `package github 
<https://github.com/SyYunP/hmf2smf>`_.

Background
---------------------

`hmf2smf` is build under the stellar mass-halo mass relation,

.. math::

   M_\star = f_\star^\text{eff} \frac{\Omega_b}{\Omega_m} M_\text{dm},

which is redshift-dependent as :math:`\Omega_b` and :math:`\Omega_m` are the cosmic baryon and matter density at the desired redshift and :math:`f_\star^\text{eff}` is the *efficiency of star formation*.
More specifically, `hmf2smf` is here to solve :math:`f_\star^\text{eff}` under different assumptions.

That being said, current version of `hmf2smf` only provide :math:`f_\star^\text{eff}` under the assumptions that the gas mass within each halo

.. math::

   M_g = \frac{\Omega_b}{\Omega_m}M_h

and the gas binding energy equals to the supernova energy

.. math::

   E_\text{binding} = E_\text{SN}

Therefore the star formation efficiency

.. math::

   f_{\star,\text{upper}}^\text{eff} = \frac{v_c^2}{v_c^2+\gamma f_g E_{51}}

where :math:`v_c` is the circular velocity (in the package we assume this equals to 
the virial velocity of halo); :math:`\gamma`, :math:`f_g` and :math:`E_{51}` is the 
supernova rate, fraction of energy contributed to the ambient gas and individual 
supernova energy, respectively.

The package is then designed to solve this star formation efficiency and transform the halo mass to the stellar mass.

Quickstart
---------------------
To install the code, simply

.. code::

   git clone https://github.com/SyYunP/hmf2smf.git

Before start using the code, users are expected to create their own cosmology with 
`Colossus <https://bdiemer.bitbucket.io/colossus/index.html>`_. For instance, one can 
set the cosmology as one of the example cosmology

   >>> from colossus.cosmology import cosmology as c
   >>> cosmo = c.setCosmology('planck18')

Next, we calculates the supernova rate of any given initial mass function. One can define 
their own initial mass function or use common initial mass functions contained in the package.

   >>> snr = snrate('salpeter')

or 

   >>> def my_salpeter(m):
   >>>     return m**-2.35
   >>> snr = snrate(my_salpeter)

For now, we provide solutions of `Salpeter (1955) <https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S/abstract>`_ 
and `Kroupa (2002) <https://ui.adsabs.harvard.edu/abs/2002Sci...295...82K/abstract>`_.

Next, one can calculate the star formation efficiency for any given 
redshift, halo mass, cosmology and the supernova rate. 

   >>> z = 5
   >>> halo_mass = 1e10 # Msun
   >>> sfe = SFE_SNFeedback(z,halo_mass,cosmo,snr)

Finally, we can calculate the stellar mass.

   >>> stellarmass_func(halo_mass,sfe,cosmo)

.. toctree::
   :maxdepth: 2
   :caption: Contents:


.. automodapi:: hmf2smf