// Program definition files for calculating the amplification factor for the
// Singular Isothermal Sphere (SIS) lens model in both the wave and geometric
// optic styles
//
// Mick Wright 2021

#include "src/sis.h"

// Function computes the amplification factor for an axially symmetric singular
// isothermal sphere (SIS) style lensing using full wave optics for given
// values of dimensionless frequency and source position. It does this using a
// summation method with the infinite sum approximated up to a given threshold
// and with a given arithmetic precision
void AmplificationFactorCalculation(acb_t amplification_factor,
                                    double dimensionless_frequency,
                                    double source_position,
                                    slong summation_upper_limit,
                                    slong precision) {
    // For the given value of impact parameter, calculate the phi, the phase
    // constant used to obtain a minimum time delay of zero
    double phi = source_position + 1./2.;

    // Calculate the prefactor term for the summation. This is given by
    // exp(i(w/2 * (y^2 + 2phi(y))))
    acb_t prefactor;
    acb_init(prefactor);

    double prefactor_dim_freq_term = dimensionless_frequency/2.;
    double prefactor_sour_pos_term = source_position * source_position;
    double prefactor_phi_term = 2. * phi;

    double prefactor_exponent_imag = prefactor_dim_freq_term
                                     * (prefactor_sour_pos_term
                                     + prefactor_phi_term);

    acb_set_d_d(prefactor, 0, prefactor_exponent_imag);
    acb_exp(prefactor, prefactor, precision);

    // Initialise the summation term
    acb_t summation_term;
    acb_init(summation_term);
    acb_zero(summation_term);

    // Loop through the summation from 0 up to the upper summation limit
    for (int n=0; n <= summation_upper_limit; n++) {
        // Calculate the gamma term, gamma(1+n/2)/n!
        double gamma_term_val = 1 + n/2.;

        acb_t gamma_term;
        acb_init(gamma_term);
        acb_set_d(gamma_term, gamma_term_val);
        acb_gamma(gamma_term, gamma_term, precision);

        arb_t factorial;
        arb_init(factorial);
        arb_fac_ui(factorial, n, precision);

        acb_div_arb(gamma_term, gamma_term, factorial, precision);

        // Calculate the power term (2w e^i3pi/2)^n/2
        double power_term_prefactor_val = 2 * dimensionless_frequency;
        double power_term_exponent_val = 3 * M_PI/2.;
        double power_term_power_val = n/2.;

        acb_t power_term;
        acb_t power_term_exponent;
        arb_t power_term_power;

        acb_init(power_term);
        acb_init(power_term_exponent);
        arb_init(power_term_power);

        acb_set_d(power_term, power_term_prefactor_val);
        acb_set_d_d(power_term_exponent, 0, power_term_exponent_val);
        arb_set_d(power_term_power, power_term_power_val);

        acb_exp(power_term_exponent, power_term_exponent, precision);
        acb_mul(power_term, power_term, power_term_exponent, precision);
        acb_pow_arb(power_term, power_term, power_term_power, precision);

        // Calculate the Hypergeometric function term:
        // 1f1(1 + n/2; 1; -i/2 wy^2)
        double hyper_arg_a_real = 1 + n/2.;
        double hyper_arg_z_imag = (-1./2.) * dimensionless_frequency
                                  * source_position * source_position;

        acb_t hyper_arg_a;
        acb_t hyper_arg_b;
        acb_t hyper_arg_z;

        acb_init(hyper_arg_a);
        acb_init(hyper_arg_b);
        acb_init(hyper_arg_z);

        acb_set_d(hyper_arg_a, hyper_arg_a_real);
        acb_one(hyper_arg_b);
        acb_set_d_d(hyper_arg_z, 0, hyper_arg_z_imag);

        acb_t hyper_term;
        acb_init(hyper_term);
        acb_hypgeom_1f1(
           hyper_term, hyper_arg_a, hyper_arg_b, hyper_arg_z, 0, precision);

        // Construct the value to add to the summation
        acb_t summation_temp;
        acb_init(summation_temp);
        acb_mul(summation_temp, gamma_term, power_term, precision);
        acb_mul(summation_temp, summation_temp, hyper_term, precision);

        acb_add(summation_term, summation_term, summation_temp, precision);

        // Memory Management - clear the declared acbs
        acb_clear(gamma_term);
        arb_clear(factorial);
        acb_clear(power_term);
        acb_clear(power_term_exponent);
        arb_clear(power_term_power);
        acb_clear(hyper_arg_a);
        acb_clear(hyper_arg_b);
        acb_clear(hyper_arg_z);
        acb_clear(hyper_term);
        acb_clear(summation_temp);
    }

    // Get the final result by multiplying the prefactor and the summation term
    acb_mul(amplification_factor, prefactor, summation_term, precision);

    // Memory Management - clear the remaining acbs
    acb_clear(prefactor);
    acb_clear(summation_term);
}

// Function computes the amplification factor for an axially symmetric singular
// isothermal sphere (SIS) style lens using the geometric optics approximation
// for given values of dimensionless frequency and source position
std::complex<double> AmplificationFactorGeometric(
    double dimensionless_frequency, double source_position) {
    // Using the complex_literals i operator
    using std::literals::complex_literals::operator""i;

    // Compute the positive magnification and square root
    double mag_plus = sqrt(abs(1. + 1./source_position));

    // If the source position is greater than 1, return the root of the
    // positive magnification
    if (source_position >= 1) {
        std::complex<double> geometric_factor = mag_plus;
        return geometric_factor;
    }

    // Compute the negative magnification and square root
    double mag_minus = sqrt(abs(1. - 1./source_position));

    // Calculate the time delay
    double time_delay = 2 * source_position;

    // Calculate the exponential term
    std::complex<double> exponent = 1i * dimensionless_frequency * time_delay;
    std::complex<double> exp_term = exp(exponent);

    // Construct the final amplification factor
    std::complex<double> geometric_factor =
        mag_plus - 1i * mag_minus * exp_term;

    return geometric_factor;
}

// Wrapper function converting the amplification factor result calculated by
// AmplificationFactorGeometric to a pair of floats for compatibility with
// ctypes for python
double * AFGRealOnly(
    double dimensionless_frequency,
    double source_position) {

    // Get the amplification factor value in the manner expected
    std::complex<double> afg = AmplificationFactorGeometric(
        dimensionless_frequency, source_position);

    // Split into real and imaginary components and place in double array
    double * afg_arr = new double[2];
    afg_arr[0] = std::real(afg);
    afg_arr[1] = std::imag(afg);

    return afg_arr;
}

// Function destroys the object placed within, used for deallocating the memory
// used by the above function
void destroyObj(double * object) {
    delete object;
    object = NULL;
}
