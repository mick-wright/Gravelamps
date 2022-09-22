'''
Millilensing Model Agnostic Functions

These functions perform calculations for a millilensing system using a phenomenological approach
that is lens model agnostic and reliant only on the lensing observables.

Written by Ania Liu 2022
'''

import numpy as np

def gather_parameter_lists(lens_parameters, parameters):
    '''
    Input:
        lens_parameters - List of ther lens parameters needed for the model
        parameters - Dictionary containing each parameter value

    Output:
        image_times - List of image time delay values
        luminosity_distances - List of luminosity distances to the millisignals
        phases - List of phases of the images.

    Function gathers the lens parameters into the lists for the amplification factor function.
    '''

    luminosity_distances = [parameters["luminosity_distance"]]
    image_times = []
    phases = []

    for parameter in lens_parameters:
        if parameter.startswith("dL"):
            luminosity_distances.append(parameters[parameter])
        elif parameter.startswith("dt"):
            image_times.append(parameters[parameter])
        else:
            phases.append(parameters[parameter])

    return image_times, luminosity_distances, phases

def get_lens_parameters(waveform_arguments):
    '''
    Input:
        waveform_arguments - Dictionary of arguments to the waveform

    Output:
        lens_parameters - List of the lens parameters needed for this model

    Function generates a list of the lens parameters that will be needed based on the waveform
    arguments
    '''

    max_number_of_images = waveform_arguments["millilensing_kmax"]

    luminosity_distance_list = np.array(range(1, max_number_of_images)).astype(str)
    time_delay_list = np.array(range(1, max_number_of_images-1)).astype(str)
    phase_list = np.array(range(max_number_of_images)).astype(str)

    luminosity_distance_list = np.char.add("dL", luminosity_distance_list)
    time_delay_list = np.char.add("dt", time_delay_list)
    phase_list = np.char.add("n", phase_list)

    lens_parameters = list(["k", *luminosity_distance_list, *time_delay_list, *phase_list])

    return lens_parameters

def frequency_to_dimensionless_frequency(frequency_array):
    '''
    Input:
        frequency_array - array of frequencies over which to generate amplification factor

    Output:
        dimensionless_frequency_array - Corresponding dimensionless frequency array

    Function converts an array of frequencies to the dimensionless equivalent geometrically.
    '''

    dimensionless_frequency_array = 1j * 2 * np.pi * frequency_array

    return dimensionless_frequency_array

def amplification_factor(frequency_array,
                         number_of_images,
                         image_times,
                         luminosity_distances,
                         morse_phases):
    '''
    Inputs:
        frequency_array - array of frequencies over which to generate amplification factor
        number_of_images - integer number of millilensing images
        image_times - Array of the time delays between millisignals
        luminosity_distances - Array containing luminosity distance, followed by the distances to
                               the next millisignal from the previous
        morse_phases - morse factor of each millisignal

    Outputs:
        amplification_factor_array - array of values for the amplification factor over the
                                     frequencies specified
    '''

    amplification_factor_value = np.exp(-1j * morse_phases[0] * np.pi)
    amplification_factor_array = np.full(len(frequency_array),
                                         amplification_factor_value,
                                         dtype=complex)

    if number_of_images > 1:
        dimensionless_frequency_array = frequency_to_dimensionless_frequency(frequency_array)
        time_delays = np.cumsum(image_times)

        for idx in range(1, number_of_images):
            luminosity_distance_term = luminosity_distances[0]/luminosity_distances[idx]
            time_delay_term = dimensionless_frequency_array * time_delays[idx-1]
            morse_phase_term = 1j * morse_phases[idx] * np.pi

            amplification_factor_array += luminosity_distance_term\
                                          * np.exp(time_delay_term - morse_phase_term)

    return amplification_factor_array
