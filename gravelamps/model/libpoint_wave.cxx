// Function definitions for calculating the amplification factor for the
// isolated point mass lens model in the wave optics regime
//
// Mick Wright 2021

#include "src/point.h"

// Function computes the value of the displacement as a function of the
// dimensionless source position used in the calculation of the phase required
// for zero minimal time delay. Function is given by:
// xm = (y + sqrt(y^2 + 4))/2
//
// Input:
//      double source_position : dimensionless displacement from optical axis
//
// Output:
//      double stationary_point_min : value of 'xm' as described above
double StationaryPointMinimum(double source_position) {
    double source_position_sq = source_position * source_position;
    double stationary_point_min = (source_position +\
                                   sqrt(source_position_sq + 4))/2;
    return stationary_point_min;
}

// Function computes the value of the phase required for the minimum time delay
// induced by the lensing to be zero. This value is given by:
// phi = (xm(y) - y)^2/2 - log(xm(y))
//
// Input:
//      double source_position : dimensionless displacement from optical axis
//
// Output:
//      double phase : phase as described above
double Phase(double source_position) {
    double stationary_point_min = StationaryPointMinimum(source_position);
    double numerator_term = stationary_point_min - source_position;
    double phase = (numerator_term * numerator_term)/2 -\
                   log(stationary_point_min);
    return phase;
}

// Function computes the exponential component of the amplification factor
// for the wave optics case at the dimensionless frequency and source position
// with the exponentiation occuring with the specified precision. The resultant
// value is assigned to the given acb_t object.
//
// Input:
//      acb_t exponential_component : object for containing the complex value of
//                                    the exponential component
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//      slong precision : arithmetic precision to use in arbitrary precision
//                        calculations
void ExponentialComponent(acb_t exponential_component,
                          double dimensionless_frequency,
                          double source_position,
                          slong precision) {
    acb_t exponent;
    acb_init(exponent);

    // Calculate the real part of the exponential component. Sufficiently
    // simple as to not require arbitrary precision
    double exp_real = (M_PI * dimensionless_frequency)/4;

    // Calculate the imaginary part of the exponential component. Sufficiently
    // simple as to not require arbitrary precision
    double exp_imag_log_part = log(dimensionless_frequency/2);
    double exp_imag_phase_part = 2 * Phase(source_position);
    double exp_imag = (dimensionless_frequency/2) *\
                      (exp_imag_log_part - exp_imag_phase_part);

    // Calculate the exponential component with designated precision
    acb_set_d_d(exponent, exp_real, exp_imag);
    acb_exp(exponential_component, exponent, precision);

    // Memory clear up of defined acb_t
    acb_clear(exponent);
}

// Function computes the gamma function component of the amplification factor
// for the wave optics case at the specified dimensionless frequency with the
// specified precision. The resultant value is assigned to the given acb_t
// object.
//
// Input:
//      acb_t gamma_component : object for containing the complex value of the
//                              gamma function component
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      slong precision : arithmetic precision to use in arbitrary precision
//                        calculations
void GammaComponent(acb_t gamma_component,
                    double dimensionless_frequency,
                    slong precision) {
    acb_t gamma_argument;
    acb_init(gamma_argument);

    // Calculate the real and imaginary arguments to the gamma function
    // component
    double gamma_arg_real = 1;
    double gamma_arg_imag = -(dimensionless_frequency/2);

    acb_set_d_d(gamma_argument, gamma_arg_real, gamma_arg_imag);

    // Compute the gamma function with given arguments with the specified
    // precision
    acb_gamma(gamma_component, gamma_argument, precision);

    // Memory clear up of defined argument acb_t
    acb_clear(gamma_argument);
}

// Function computes the confluent hypergeometric component of the
// amplification factor for the wave optics case at the specified dimensionless
// frequency and source position with the specified precision. The resultant
// value is assigned to the given acb_t object
//
// Input:
//      acb_t hyper_component : object for containing the complex value of the
//                              confluent hypergeometric function component
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//      slong precision : arithmetic precision to use in arbitrary precision
//                        calculations
void HyperComponent(acb_t hyper_component,
                    double dimensionless_frequency,
                    double source_position,
                    slong precision) {
    // Calculate any of the non trivial real and imaginary parts of the
    // arguments to the confluent hypergeometric function
    double hyper_arg_a_imag = dimensionless_frequency/2;

    double source_position_sq = source_position * source_position;
    double hyper_arg_z_imag = hyper_arg_a_imag * source_position_sq;

    // Initialisation of the arguments to the confluent hypergeomtric function
    acb_t hyper_arg_a;
    acb_t hyper_arg_b;
    acb_t hyper_arg_z;

    acb_init(hyper_arg_a);
    acb_init(hyper_arg_b);
    acb_init(hyper_arg_z);

    acb_set_d_d(hyper_arg_a, 0, hyper_arg_a_imag);
    acb_one(hyper_arg_b);
    acb_set_d_d(hyper_arg_z, 0, hyper_arg_z_imag);

    // Compute the confluent hypergoemtric function value
    acb_hypgeom_1f1(hyper_component,
                    hyper_arg_a,
                    hyper_arg_b,
                    hyper_arg_z,
                    0,
                    precision);

    // Memory clear up of each defined acb_t
    acb_clear(hyper_arg_a);
    acb_clear(hyper_arg_b);
    acb_clear(hyper_arg_z);
}

// Function computes the amplification factor for an isolated axially symmetric
// point mass lens using full wave optics for given values of dimensionless
// frequency and source position with the given arithmetic precision. This
// value is given to the provided amplification factor object.
//
// Input:
//      acb_t amplification_factor : object for containing complex
//                                   amplification factor value
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//      slong precision : arithemtic precision to use in arbitrary precision
//                        calculations
void AmplificationFactorWave(acb_t amplification_factor,
                             double dimensionless_frequency,
                             double source_position,
                             slong precision) {
    // Calculate the exponential component of the amplification factor
    acb_t exponential_component;
    acb_init(exponential_component);
    ExponentialComponent(exponential_component,
                         dimensionless_frequency,
                         source_position,
                         precision);

    // Calculate the gamma function component of the amplification factor
    acb_t gamma_component;
    acb_init(gamma_component);
    GammaComponent(gamma_component,
                   dimensionless_frequency,
                   precision);

    // Calculate the confluent hypergeometric function component of the
    // amplification factor
    acb_t hyper_component;
    acb_init(hyper_component);
    HyperComponent(hyper_component,
                   dimensionless_frequency,
                   source_position,
                   precision);

    // Construct the final value by multiplying the components
    acb_t exponential_gamma;
    acb_init(exponential_gamma);

    acb_mul(exponential_gamma,
            exponential_component,
            gamma_component,
            precision);
    acb_mul(amplification_factor,
            exponential_gamma,
            hyper_component,
            precision);

    // Memory clear up of declared acb_t objects
    acb_clear(exponential_component);
    acb_clear(gamma_component);
    acb_clear(hyper_component);
    acb_clear(exponential_gamma);
}
