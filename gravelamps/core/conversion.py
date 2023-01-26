"""Gravelamps Conversion Functions

Following are functions handling physical unit conversions

Written by Mick Wright 2022

Routines
--------
lens_mass_to_redshifted_lens_mass
    Converts from source frame lens mass to observer frame redshfited lens msas
solar_mass_to_natural_mass
    Converts mass from units of solar mass to natural units such that c=G=1
frequency_to_dimensionless_frequency
    Converts regular frequencies in Hz to dimensionless frequencies

"""

from astropy.constants import c, G, M_sun
import numpy as np

from bilby.gw.conversion import luminosity_distance_to_redshift

def lens_mass_to_redshifted_lens_mass(lens_mass, lens_fractional_distance, luminosity_distance):
    """
    Converts from source frame lens mass to observer frame redshifted lens mass

    Parameters
    ----------
    lens_mass : float
        Source frame mass of lensing object in Solar Masses
    lens_fractional_distance : float
        Distance to lensing object as a fraction of the luminosity distance to the source
    luminosity_distance : float
        Distance to the source object in Mpc

    Returns
    -------
    redshifted_lens_mass : float
        Observer frame mass of lensing object in natural units
    """

    lens_distance = lens_fractional_distance * luminosity_distance
    lens_redshift = luminosity_distance_to_redshift(lens_distance)
    redshifted_lens_mass = solar_mass_to_natural_mass(lens_mass * (1 + lens_redshift))

    return redshifted_lens_mass

def solar_mass_to_natural_mass(solar_mass):
    """
    Converts mass from units of solar mass to natural units such that c=G=1.

    Parameters
    ----------
    solar_mass : float
        Mass in solar mass units

    Returns
    -------
    natural_mass : float
        The same mass in natural units
    """

    kg_mass = solar_mass * M_sun
    natural_mass = kg_mass * (G/c**3)

    return natural_mass.value

def frequency_to_dimensionless_frequency(frequency_array,
                                         redshifted_lens_mass=None):
    """
    Converts regular frequencies in Hz to dimensionless frequencies.

    Can be done in two methods. If the redshifted lens mass is supplied this will be done for that
    specific mass object. If it is not supplied, it will be done using a generic agnostic form

    Parameters
    ----------
    frequency_array : ArrayLike
        Array of frequencies to convert
    redshifted_lens_mass : float, optional
        Mass of lensing object to use in conversion

    Returns
    -------
    dimensionless_frequency_array : ArrayLike
        Array of converted dimensionless frequencies
    """

    if redshifted_lens_mass is not None:
        dimensionless_frequency_array = 8 * np.pi * redshifted_lens_mass * frequency_array
    else:
        dimensionless_frequency_array = 1j * 2 * np.pi * frequency_array

    return dimensionless_frequency_array
