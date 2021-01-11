import numpy as np
import mpmath as mp
import scipy.special as scisp

##Lensing Array Calculation Functions
def impact_normalisation_ratio_min(y):
    '''
    Inputs:
        y - variable y corresponding to (distance_to_lens/(length_normalisation_constant*distance_to_source))*source_position

    Outputs:
        xm - value of impact_parameter/length_normalisation_constant for the phase constant corresponding to the minimum time delay i.e. that of the image that travels the shortest path to the observer

    Function takes variable y, and calculates the corresponding value of xm
    '''

    xm = (y + np.sqrt(y**2 + 4))/2
    return xm

def phase_const(y):
    '''
    Inputs:
        y - variable y corresponding to (distance_to_lens/(length_normalisation_constant*distance_to_source))*source_position

    Outputs:
        phi - value of the phase constant used to obtain the minimum time delay produced by the lensing

    Function takes the variable y and calculates the value of the phase constant used to obtain the minimum time delay induced by the lensing.
    '''

    phi = ((impact_normalisation_ratio_min(y)-y)**2)/2 - np.log(impact_normalisation_ratio_min(y))

    return phi

def amplification_factor(w,y):
    '''
    Inputs:
        w - dimensionless frequency array over which the amplification factor is to be calculated
        y - value of the variable y which corresponds to the gravitational wave source position modified by a factor of (distance_to_lens/length_normalisation_constant*distance_to_source)

    Outputs:
        f - value of the amplification factor for the given values of w and y

    Function takes given values of w and y and computes the amplification factor for a lensing object that can be considered a point mass and that is axially-symmetric
    '''

    exponent = (np.pi*w)/4 + 1j*(w/2)*(np.log(w/2)-2*phase_const(y))

    exp_array = np.vectorize(mp.exp)
    exponential_component = exp_array(exponent)

    gamma_component = scisp.gamma(1-(1j/2)*w)

    hyper_arg_a = (1j/2)*w
    hyper_arg_b = 1
    hyper_arg_z = hyper_arg_a*(y**2)

    hyper_array = np.vectorize(mp.hyp1f1, excluded="maxterms")
    hyper_component = hyper_array(hyper_arg_a, hyper_arg_b, hyper_arg_z, maxterms=10**20)

    f = exponential_component*gamma_component*hyper_component

    return f

def generate_amplification_factor_matrix(dim_freq_array, impact_array, lens_model="Point"):
    '''Inputs:
            dim_freq_array - Array of dimensionless_frequency to generate the amplification factor over
            impact_array - Array of impact parameter variable y, to generate the amplification factor over
            lens_model - Which lens model to use, defaults to point mass.

        Outputs:
            amp_fac_matrix - Matrix containing all values of the amplification factor calculated

        Function takes arrays of dimensionless frequency and impact parameter and calculates the amplification factor for all values using the specified lens model.'''

    amp_fac_func = np.vectorize(amplification_factor, excluded="w")

    amp_fac_matrix = np.stack(amp_fac_func(w=dim_freq_array, y=impact_array), axis=1).astype(complex)

    return(amp_fac_matrix)
