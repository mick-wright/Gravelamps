import numpy as np
import astropy.constants as const
import scipy.interpolate as scint

def dimensionless_frequency(f, Mlz):
    '''Inputs:
            f - frequency array to be converted
            Mlz - redshifted lens mass of lens

       Outputs:
            w - dimensionless frequency array

        Function takes a frequency array and converts to the dimensionless frequency corresponding to the input frequency array and a lens with redshifted lens mass (Mlz = M_lens(1+z))'''
    return(8*np.pi*Mlz*f)

def natural_mass(mass, mode="solar"):
    '''Inputs:
            mass - the mass of the object to be converted
            mode - units of the mass given in mass, this can be either solar or kg, default being solar

       Outputs:
            m_nat - the mass in natural units corresponding to input mass

    Function takes a given value of mass in either solar masses or kg and converts it to the equivalent natural units. This is done by multiplying the value of the mass in kg by a factor of (G/c**3)'''

    if mode == "solar":
        m_kg = mass*const.M_sun
    else:
        m_kg = mass

    m_nat = m_kg*(const.G/const.c**3)

    return(m_nat.value)

def generate_dimensionless_frequency_array(highest_freq, Mlz, lowest_freq=0, changeover_frequency=10, below_change_npoints=2000, above_change_npoints=8000):
    '''Inputs:
            highest_freq - the upper bound of the frequencies that you are interested in
            Mlz - the redshifted lens mass of the object doing the lensing
            lowest_freq - the lower bound of the frequencies that you are interested in, default is 0, though, if it is 0, will instead do the first point after 0
            changeover_frequency - the point at which the points may be investigated less frequently, defaults to 10
            below_change_npoints - the number of points between the value of lowest_freq and changeover_frequency to evaluate, defaults to 2000
            above_change_npoints - the number of points between the value of changeover_frequency and highest_freq to evaluate, defaults to 8000

       Outputs:
            dim_freq_array - dimensionless frequency array from lowest freq to highest freq

        Function takes at minimum a highest frequency and generates a number of points to evaluate the dimensionless frequency of. It generates two linear spaces to do this, due to the need to investigate the lower band in more detail than the higher. User is able to customise the point that is considered to be the point below which greater resolution is needed, and the amount of points generated for each side.'''

    #Amplification Factor maths cannot process a value of 0 for w, avoiding by shifting to the next data point
    if lowest_freq == 0:
        lowest_freq += 1/below_change_npoints

    lower_frequency_array = np.linspace(lowest_freq, changeover_frequency, below_change_npoints)
    #Avoiding calculating the same point twice
    higher_frequency_array = np.linspace(changeover_frequency+1/above_change_npoints, highest_freq, above_change_npoints)

    complete_frequency_array = np.concatenate((lower_frequency_array, higher_frequency_array))

    dim_freq_array = dimensionless_frequency(complete_frequency_array, Mlz)

    return(dim_freq_array)

def generate_amplification_factor_matrix(dim_freq_array, impact_array, lens_model="Point"):
    '''Inputs:
            dim_freq_array - Array of dimensionless_frequency to generate the amplification factor over
            impact_array - Array of impact parameter variable y, to generate the amplification factor over
            lens_model - Which lens model to use, defaults to point mass.

        Outputs:
            amp_fac_matrix - Matrix containing all values of the amplification factor calculated

        Function takes arrays of dimensionless frequency and impact parameter and calculates the amplification factor for all values using the specified lens model.'''

    if lens_model == "Point":
        amp_fac_func = np.vectorize(pl.amplification_factor, excluded="w")
    else:
        raise ValueError("Lens Model not recognised!")

    amp_fac_matrix = np.stack(amp_fac_func(w=dim_freq_array, y=impact_array), axis=1).astype(complex)

    return(amp_fac_matrix)

def generate_interpolator(dim_freq_array, impact_array, amp_fac_matrix):
    '''Inputs:
            dim_freq_array - Array of the dimensionless frequencies used to generate the rows of the amplification factor matrix
            impact_array - Array of the impact parameters used to generate the columns of the amplification factor matrix
            amp_fac_matrix - Matrix of amplification factor values from which to generate an interpolator

       Outputs:
            interpolator_func - interpolation function based on values given, will be called with syntax interpolator_func(w,y) for given values of dimensionless frequency and impact parameter

        Function generates a two dimensional complex interpolator for the value of the amplification factor outside of the generated values.'''

    amp_fac_real = np.real(amp_fac_matrix)
    amp_fac_imag = np.imag(amp_fac_matrix)

    amp_interp_real = scint.RectBivariateSpline(dim_freq_array, impact_array, amp_fac_real, kx=1, ky=1, s=0)
    amp_interp_imag = scint.RectBivariateSpline(dim_freq_array, impact_array, amp_fac_imag, kx=1, ky=1, s=0)

    interpolator_func = lambda w, y: (amp_interp_real(w,y) + 1j*amp_interp_imag(w,y)).flatten()

    return(interpolator_func)
