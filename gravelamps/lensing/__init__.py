'''
Gravelamps Lensing Module
Mick Wright 2021

Module contains the functions necessary to construct lensed waveforms as well as generate the data
necessary to perform analysis runs for both local machines and for cluster machines utilising the
HTCondor scheduler platform
'''

import importlib

import numpy as np
import astropy.constants as const

import bilby
from bilby.gw.waveform_generator import WaveformGenerator

from . import utils

class LensedWaveformGenerator(WaveformGenerator):
    '''
    Lensed Waveform Generator Class

    Based upon the general bilby Waveform Generator class, but with the additional requirements for
    files containing the dimensionless frequency and source position arrays as well as the
    amplification factor matrices. Using these, it will then generate the interpolator for use
    by the frequency domain source model
    '''

    def __init__(self, duration=None, sampling_frequency=None, start_time=0,
                 frequency_domain_source_model=None, time_domain_source_model=None,
                 parameters=None, parameter_conversion=None, waveform_arguments=None):

        #Initialise the parent class
        super().__init__(duration, sampling_frequency, start_time, frequency_domain_source_model,
                         time_domain_source_model, parameters, parameter_conversion,
                         waveform_arguments)

        if waveform_arguments["methodology"] == "interpolate":
            #Extract the files necessary to generate the interpolator
            dimensionless_frequency_file = waveform_arguments["dim_freq_file"]
            source_position_file = waveform_arguments["sour_pos_file"]
            amplification_factor_real_file = waveform_arguments["amp_fac_real_file"]
            amplification_factor_imag_file = waveform_arguments["amp_fac_imag_file"]

            #Generate the interpolator and add it to the waveform_arguments dictionary
            waveform_arguments["interpolator"] = utils.generate_interpolator(
                dimensionless_frequency_file, source_position_file, amplification_factor_real_file,
                amplification_factor_imag_file)

        elif waveform_arguments["methodology"] == "direct":
            #Extract the lens model from the waveform arguments
            lens_model = waveform_arguments["lens_model"]

            lens_module_spec = importlib.util.find_spec(lens_model)

            if lens_module_spec is None:
                if lens_model.endswith("lens"):
                    lens_module_spec = importlib.util.find_spec(lens_model[:-4])
            if lens_module_spec is None:
                raise ModuleNotFoundError(f"No module corresponding to {lens_model} found!")

            try:
                waveform_arguments["amplification_factor_function"] = getattr(
                    importlib.import_module(lens_module_spec.name), "amplification_factor")
            except NameError:
                print(f"amplification_factor not defined within {lens_module_spec.name}")

def BBH_lensed_waveform(frequency_array, mass_1, mass_2, a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl,
                        luminosity_distance, theta_jn, phase, ra, dec, geocent_time, psi,
                        lens_mass, source_position, lens_fractional_distance, **kwargs):
    '''
    Inputs:
        frequency_array - frequencies over which to generate waveform
        mass_1 - non-redshifted primary mass in solar masses
        mass_2 - non-redshifted secondary mass in solar masses
        a_1 - dimensionless spin magnitude of the primary
        a_2 - dimensionless spin magnitude of the secondary
        tilt_1 - polar angle between primary spin and the orbital angular momentum in radians
        tilt_2 - polar angle between secondary spin and the orbital angular momentum in radians
        phi_12 - azimuthal angle bbetween primary and secondary spins in radians
        phi_jl - azimuthal angle between total angular momentum and the orbital angular momentum
                 in radians
        luminosity_distance - luminosity distance to the source in Mpc
        theta_jn - inclination angle between line-of-sight and the orbital angular momentum in
                   radians
        phase - phase in radians
        ra - right ascension of the source in radians
        dec - declination of the source in radians
        geocent_time - time of coalescene or peak amplitude in GPS seconds
        psi - gravitational wave polarisation angle in radians
        lens_mass - non-redshifted mass of the lensing object in solar masses
        source_position - distance from observer plane of the source
        lens_fractional_distance - fractional position of the lens compared to the luminosity
                                   distance

        **kwargs:
            waveform_approximant - the waveform approximant used to generate teh base waveform
            reference_frequency - the waveform reference frequency in Hz
            minimum_frequency - the waveform minimum frequency in Hz
            maximum_frequency - the waveform maximum frequency in Hz
            methodology - Can be either 'interpolate' or 'direct' determining which method of
                          caluclation will be used to generate the amplification factor
            interpolator - interpolating function to generate the amplification factor used to lens
                           the base waveform when using the interpolate methodology
            amplification_factor_func - function to generate the amplification factor used to lens
                                        the base waveform when using the direct methodology
            scaling_cosntant - constant used for the Navarro Frenk White (NFW) lens model
            lens_model - which lens model is being used
            image_position_func - function solving the lens equation for image positions needed for
                                  the NFW direct method
            min_time_delay_phase_func - function calculating the phase needed for a minimum time
                                        delay necessary for the NFW direct method
    Outputs:
        lens_waveform - dictionary containing the plus and cross polarisation mode strain data for
                        the lensed waveform

    Function takes a base waveform generated by the lal_binary_black_hole function and then lenses
    it by a lens with given mass, at the given source position. The amplification factor is
    generated by m,eans of the given interpolator function
    '''

    #Generate the waveform kwargs dict and then update it using the given kwargs
    waveform_kwargs = dict(waveform_appoximant="IMRPhenomPv2", reference_frequency=50,
                           minimum_frequency=20, maximum_frequency=1024, interpolator=None,
                           amplification_factor_func=None, scaling_constant=None,
                           lens_model=None, image_position_func=None,
                           min_time_delay_phase_func=None, methodology="interpolate")
    waveform_kwargs.update(kwargs)

    #Extract components
    waveform_approximant = waveform_kwargs["waveform_approximant"]
    reference_frequency = waveform_kwargs["reference_frequency"]
    minimum_frequency = waveform_kwargs["minimum_frequency"]
    maximum_frequency = waveform_kwargs["maximum_frequency"]
    interpolator = waveform_kwargs["interpolator"]
    amplification_factor_function = waveform_kwargs["amplification_factor_function"]
    methodology = waveform_kwargs["methodology"]

    #Construct the base waveform using lal_binary_black_hole
    base_waveform = bilby.gw.source.lal_binary_black_hole(frequency_array, mass_1=mass_1,
        mass_2=mass_2, a_1=a_1, a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2, phi_12=phi_12, phi_jl=phi_jl,
        luminosity_distance=luminosity_distance, theta_jn=theta_jn, phase=phase,
        waveform_approximant=waveform_approximant, reference_frequency=reference_frequency,
        minimum_frequency=minimum_frequency, ra=ra, dec=dec, geocent_time=geocent_time, psi=psi,
        maximum_frequency=maximum_frequency)

    if base_waveform is None:
        return None

    #For the point estimate of luminosity distance calculate the redshifted lens mass
    lens_distance = lens_fractional_distance * luminosity_distance
    lens_redshift = bilby.gw.conversion.luminosity_distance_to_redshift(lens_distance)
    redshifted_lens_mass = natural_mass(lens_mass * (1 + lens_redshift))

    #Use this to construct the dimensionless frequency array
    dimensionless_frequency_array = dimensionless_frequency(frequency_array, redshifted_lens_mass)

    #Get the amplification factor function and hence generate the amplification factor
    if methodology == "interpolate":
        if interpolator is None:
            raise ValueError("To use interpolate method, interpolator must be given!")
        amplification_factor_array = interpolator(dimensionless_frequency_array, source_position)
    else:
        if amplification_factor_function is None:
            raise ValueError("To use direct method, direct calculation function must be given!")
        amplification_factor_array = amplification_factor_function(dimensionless_frequency_array,
                                                                   source_position)

    #Create the lens waveform by multiplying the base waveform by the amplification factor array
    lens_waveform = {}

    for polarisation in base_waveform:
        lens_waveform[polarisation] = np.multiply(base_waveform[polarisation],
                                                  amplification_factor_array)

    #Return the lensed waveform
    return lens_waveform

def dimensionless_frequency(frequency, redshifted_lens_mass):
    '''
    Inputs:
        frequency - the frequencies to be converted
        redsfhited_lens_mass - the redshifted mass of the lensing object in natural units

    Outputs:
        dim_freq - array containing the dimensionless frequencies calculated from the input
                   frequencies

    Function takes an array of frequencies and converts to the dimensionless frequency array
    corresponding to a lens with the given redshifted mass. This is given by 8*pi*f*M_lz
    '''

    dim_freq = 8 * np.pi * redshifted_lens_mass * frequency

    return dim_freq

def natural_mass(mass, mode="solar"):
    '''
    Inputs:
        mass - the mass to be converted
        mode - units of the mass given. Can be either "solar" for Solar Masses or "kg"

    Outputs:
        m_nat - the mass in natural units

    Function converts a given mass from either solar masses or kg to natural units
    '''

    #If the mass is in solar masses, first convert to kg
    if mode == "solar":
        m_kg = mass * const.M_sun
    else:
        m_kg = mass

    #Convert from kg to natural units
    m_nat = m_kg * (const.G/const.c**3)

    return m_nat.value
