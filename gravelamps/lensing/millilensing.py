"""Millilensing Model Agnostic Functions

Following are functions performing calculations for a mililensing system using a phenomenological
model-agnostic approach reliant only on the lensing observables.

Written by Ania Liu 2022

Routines
--------
gather_parameter_lists
    Gathers lens parameters into lists
get_lens_parameters
    Generates required parameters based on waveform arguments
frequency_to_dimensionless_frequency
    Conversion from frequencies to dimensionless frequencies in model agnostic fashion
amplification_factor
    Computs the amplification factor for the given parameters
"""

import numpy as np

def gather_parameter_lists(lens_parameters, parameters):
    """Gathers lens parameters into lists

    This gathers the image times, lumionsity distances, and morse phases of the individual
    millisignals into lists of each of these for management. This is done in time ordering

    Parameters
    ----------
    lens_parameters : list of strings
        Every lensing observable that is needed for the total millilensing signal
    parameters : dict
        Contains values for each of the specified parameters above

    Returns
    -------
    image_times : list of floats
        Image time delay values as compared to first image
    luminosity_distances : list of floats
        Luminosity distances to millisignals
    phases : list of floats
        Phases of each millisignal
    """

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
    """
    Generates required parameters based on waveform arguments

    The waveform arguments specify a number of signals, this will generate the required
    parameter names for each millisignal --- each needs a luminosty distance, time delay,
    and phase. The first signal requires only the phase. It's luminosity distance is in the
    main parameter list, and time delays are measured relative to the first image.

    Parameters
    ----------
    waveform_arugments : dict
        Contains arguments for the waveform generation

    Returns
    -------
    lens_parameters : list of strings
        Every lensing observable that is needed for the total millilensing signal
    """

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
    """
    Conversion of frequencies to dimensionless frequencies in model agnostic fashion

    Parameters
    ----------
    frequency_array : Array of floats
        Frequencies of interest

    Returns
    -------
    dimensionless_frequency_array : Array of floats
        Dimensionless form of the input frequencies from model agnostic conversion
    """

    dimensionless_frequency_array = 1j * 2 * np.pi * frequency_array

    return dimensionless_frequency_array

def amplification_factor(frequency_array,
                         number_of_images,
                         image_times,
                         luminosity_distances,
                         morse_phases):
    """
    Computes the amplification factor for the given parameters.

    The total amplification of a millilensing signal is given as the individual amplifications
    from each image (or millisignal) contributing. The number of these must be specified and
    lensing observables given for each of the images in question.

    Parameters
    ----------
    frequency_array : Array of floats
        Frequencies to be amplified
    number_of_images : float
        The total number of millisignals that comprise the total amplification
    image_times : Array of floats
        Time delays between millisignals relative to the first image
    luminosity_distances : Array of floats
        Luminosity distance of the first image, followed by the distance to the next millisignal
        from the previous
    morse_phases : Array of floats
        Morse factor of each millisignal

    Returns
    -------
    amplification_factor_array : Array of floats
        Values of the amplification factor over the specified frequencies
    """

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
