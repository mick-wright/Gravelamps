// Function definitions for caluclating the amplification factor for the
// isolated singular isothermal sphere (SIS) lens model in the wave optics
// regime
//
// Mick Wright 2021

#include "src/sis.h"

// Function computes the value of the phase required for the minimum time delay
// induced by the lensing to be zero. This value is given by:
// phi = y + 1/2
//
// Input:
//      double source_position : dimensionless displacement from optical axis
//
// Output:
//      double phase : phase as described above
double Phase(double source_position) {
    double phase = source_position + 1./2.;
    return phase;
}

// Function computes the prefactor component for the amplification factor for
// the wave optics case at the dimensionless frequency and source position
// specified with the exponentiation occuring with the specified precision.
// The resultant value is assigned to the given acb_t object
//
// Input:
//      acb_t prefactor_term : object for containing the complex value of the
//                             calculated prefactor
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//      slong precision : arithmetic precision to use in arbitrary precision
//                        calculations
void PrefactorCalculation(acb_t prefactor,
                          double dimensionless_frequency,
                          double source_position,
                          slong precision) {
    double dimensionless_frequency_term = dimensionless_frequency/2.;
    double source_position_sq = source_position * source_position;
    double phi_term = 2. * Phase(source_position);

    double prefactor_exponent_imag = dimensionless_frequency_term
                                     * (source_position_sq
                                     + phi_term);

    acb_set_d_d(prefactor, 0, prefactor_exponent_imag);
    acb_exp(prefactor, prefactor, precision);
}

// Function computes the gamma function component of the nth summation term for
// the given value of n with the desired precision. The resultant value is
// assigned to the given acb_t object.
//
// Input:
//      acb_t gamma_term : object for containing the complex value of the gamma
//                         function term
//      int n : which term in the summation is being calculated
//      slong precision : arithmetic precision to use in arbitrary precision
//                        calculations
void GammaComponent(acb_t gamma_term, int n, slong precision) {
    double gamma_argument = 1 + n/2.;
    acb_set_d(gamma_term, gamma_argument);
    acb_gamma(gamma_term, gamma_term, precision);

    arb_t n_factorial;
    arb_init(n_factorial);
    arb_fac_ui(n_factorial, n, precision);

    acb_div_arb(gamma_term, gamma_term, n_factorial, precision);

    arb_clear(n_factorial);
}

// Function computes the power compoonent of the nth summation term for the
// given value of dimensionless frequency and n with the desired precision. The
// resultant value is assigned to the given acb_t object
//
// Input:
//      acb_t power_term : object for containing the complex value of the power
//                         component
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      int n : which term in the summation is being calculated
//      slong precision : arithmetic precision to use in arbitrary precision
//                        calculations
void PowerComponent(acb_t power_term,
                    double dimensionless_frequency,
                    int n,
                    slong precision) {
    double dimensionless_frequency_term = 2 * dimensionless_frequency;
    acb_set_d(power_term, dimensionless_frequency_term);

    double exponential_term_imag = 3 * M_PI/2.;
    acb_t exponential_term;
    acb_init(exponential_term);
    acb_set_d_d(exponential_term, 0, exponential_term_imag);

    double power_value = n/2.;
    arb_t power;
    arb_init(power);
    arb_set_d(power, power_value);

    acb_exp(exponential_term, exponential_term, precision);
    acb_mul(power_term, power_term, exponential_term, precision);
    acb_pow_arb(power_term, power_term, power, precision);

    acb_clear(exponential_term);
    arb_clear(power);
}

// Function computes the confluent hypergeometric function component of the nth
// summation term for the given dimensionless frequency, source position, and n
// with the desired preicsion. The resultant value is assigned to the given
// acb_t object.
//
// Input:
//      acb_t hyper_term : object for containing the complex value of the
//                         confluent hypergeometric function component
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//      int n : which term in the summation is being calculated
//      slong precision : arithmetic precision to use in arbitrary precision
//                        calculations
void HyperComponent(acb_t hyper_term,
                    double dimensionless_frequency,
                    double source_position,
                    int n,
                    slong precision) {
    // Get the values of those arguments to the function that have non-trivial
    // value
    double hyper_arg_a_real = 1 + n/2.;
    double hyper_arg_z_imag = (-1./2.) * dimensionless_frequency\
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

    acb_hypgeom_1f1(hyper_term,
                    hyper_arg_a,
                    hyper_arg_b,
                    hyper_arg_z,
                    0,
                    precision);

    acb_clear(hyper_arg_a);
    acb_clear(hyper_arg_b);
    acb_clear(hyper_arg_z);
}

// Function computes the nth summation term for the amplification factor for
// the wave optics case at the dimensionless frequency and source position
// specified with the desired precision. The resultant value is assigned
// to the given acb_t object
//
// Input:
//      acb_t summation_term : object for containing the complex value of the
//                             nth summation term
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//      slong precision : arithmetic precision to use in arbitrary precision
//                        calculations
//      int n : which term in the summation is being calculated
void SummationNTerm(acb_t summation_term,
                   double dimensionless_frequency,
                   double source_position,
                   slong precision,
                   int n) {
    // Calculate the gamma function component of the summation term
    acb_t gamma_term;
    acb_init(gamma_term);
    GammaComponent(gamma_term,
                   n,
                   precision);

    // Calculate the power component of the summation term
    acb_t power_term;
    acb_init(power_term);
    PowerComponent(power_term,
                   dimensionless_frequency,
                   n,
                   precision);

    // Calculate the confluent hypergoemtric function component of the
    // summation term
    acb_t hyper_term;
    acb_init(hyper_term);
    HyperComponent(hyper_term,
                   dimensionless_frequency,
                   source_position,
                   n,
                   precision);

    // Construct the final summation term from the components
    acb_mul(summation_term, gamma_term, power_term, precision);
    acb_mul(summation_term, summation_term, hyper_term, precision);

    // Clear the declared acbs to free memory
    acb_clear(gamma_term);
    acb_clear(power_term);
    acb_clear(hyper_term);
}

// Function computes the amplification factor for an isolated axially symmetric
// SIS lens using full wave optics for given values of dimensionless frequency
// and source position with the given arithmetic precision. The infinite sum
// in the analytical expression is approximated up to a specified term. The
// resultant value is given to the provided amplification factor object.
//
// Input:
//      acb_t amplification_factor : object for containing the complex
//                                   amplification factor value
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//      slong summation_upper_limit : value to terminate calculation of the
//                                    infinite sum
//      slong precision : arithmetic precision to use in arbitrary precision
//                        calculations
void AmplificationFactorWave(acb_t amplification_factor,
                             double dimensionless_frequency,
                             double source_position,
                             slong summation_upper_limit,
                             slong precision) {
    // Calculate the prefactor term
    acb_t prefactor;
    acb_init(prefactor);
    PrefactorCalculation(prefactor,
                         dimensionless_frequency,
                         source_position,
                         precision);

    // Calculate the summation term
    acb_t summation_term;
    acb_init(summation_term);
    acb_zero(summation_term);

    acb_t summation_tmp;
    acb_init(summation_tmp);
    // Loop through the summation adding the value for that term to the total
    for (int n=0; n <= summation_upper_limit; n++) {
        acb_zero(summation_tmp);
        SummationNTerm(summation_tmp,
                       dimensionless_frequency,
                       source_position,
                       precision,
                       n);
        acb_add(summation_term, summation_term, summation_tmp, precision);
    }

    // Construct final result as multiplication of the two terms
    acb_mul(amplification_factor, prefactor, summation_term, precision);

    // Clear the declared internal abcs to clear memory
    acb_clear(prefactor);
    acb_clear(summation_term);
    acb_clear(summation_tmp);
}
