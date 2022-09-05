'''
Gravelamps Conversion Functions

Functions within handle conversion of quantities from one unit to another

Written by Mick Wright 2022
'''

from astropy.constants import c, G, M_sun

import numpy as np

from bilby.gw.conversion import luminosity_distance_to_redshift

def lens_mass_to_redshifted_lens_mass(lens_mass, lens_fractional_distance, luminosity_distance):
    '''
    Input:
        lens_mass - Source frame mass of lensing object in Solar Masses
        lens_fractional_distance - Distance to lensing object as fraction of luminosty distance to
                                   source
        luminosity_distance - Distance to the source object in Mpc

    Output:
        redshifted_lens_mass - Redshifted mass of the lensing object in natural units
    '''

    lens_distance = lens_fractional_distance * luminosity_distance
    lens_redshift = luminosity_distance_to_redshift(lens_distance)
    redshifted_lens_mass = solar_mass_to_natural_mass(lens_mass * (1 + lens_redshift))

    return redshifted_lens_mass

def solar_mass_to_natural_mass(solar_mass):
    '''
    Input:
        solar_mass - Mass in units of solar mass

    Output:
        natural_mnass - Mass in natural units

    Converts a mass from solar masses to natural units.
    '''

    kg_mass = solar_mass * M_sun
    natural_mass = kg_mass * (G/c**3)

    return natural_mass.value

def frequency_to_dimensionless_frequency(frequency_array,
                                         redshifted_lens_mass=None):
    '''
    Input:
        frequency_array - Array of frequencies for gravitational wave data
        redshifted_lens_mass - Observer frame mass of the lensing object in natural units

    Output;
        dimensionless_frequency_array - Corresponding array of dimensionless frequencies
                                        associated with the frequency array and the lensing object

    Function converts an arary of frequencies to the corresponding dimensionless frequencies
    with either the specified object mass or to a generic agnostic form if object mass is not
    specified
    '''

    if redshifted_lens_mass is not None:
        dimensionless_frequency_array = 8 * np.pi * redshifted_lens_mass * frequency_array
    else:
        dimensionless_frequency_array = 1j * 2 * np.pi * frequency_array

    return dimensionless_frequency_array
