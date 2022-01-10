// Program definition files for calculating the amplification factor for the
// point mass lens model in both the wave and geometric optic styles
//
// Mick Wright 2021

#include "src/pointlens.h"

// Function computes the value of xm - the source position divided by a length
// normalisation constant for the phase constant corresponding to the minimum
// time delay i.e. that of the image that travels the shortest path to the
// observer. Function is given by xm = (y + sqrt(y^2+4))/2
double StationaryPointMinimum(double source_position) {
    double source_position_sq = source_position * source_position;
    double xm = (source_position + sqrt(source_position_sq+4))/2;
    return xm;
}

// Function computes the value of phi - the phase constant used to obtain the
// minimum time delay induced by the lensing. Function is given by
// phi = (xm(y)-y)^2 / 2 - log(xm(y))
double MinTimeDelayPhaseConstant(double source_position) {
    double xm_y = StationaryPointMinimum(source_position);
    double x_minus_y = xm_y - source_position;
    double phi = (x_minus_y*x_minus_y)/2 - log(xm_y);
    return phi;
}

// Function computes the amplification factor for an axially symmetric point
// mass lens using full wave optics for given values of dimensionless frequency
// and source position with arithmetic precision given by precision
void AmplificationFactorCalculation(acb_t amplification_factor,
                                    double dimensionless_frequency,
                                    double source_position,
                                    slong precision) {
    // Calculate the real part of the exponential component of the
    // amplification factor
    double exp_real = (M_PI * dimensionless_frequency)/4;

    // Calculate the imaginary part of the exponential component of the
    // amplification factor
    double exp_log_part = log(dimensionless_frequency/2);
    double exp_phi_part = 2 * MinTimeDelayPhaseConstant(source_position);
    double exp_imag = (dimensionless_frequency/2) * (exp_log_part-exp_phi_part);

    // From the calculations construct the exponential component
    acb_t exponent;
    acb_t exponential_component;

    acb_init(exponent);
    acb_init(exponential_component);

    acb_set_d_d(exponent, exp_real, exp_imag);
    acb_exp(exponential_component, exponent, precision);

    // Construct the real and imaginary arguments of the gamma component
    // and then calculate the value of the gamma component
    double gamma_arg_real = 1;
    double gamma_arg_imag = -(dimensionless_frequency/2);

    acb_t gamma_argument;
    acb_t gamma_component;

    acb_init(gamma_argument);
    acb_init(gamma_component);

    acb_set_d_d(gamma_argument, gamma_arg_real, gamma_arg_imag);

    acb_gamma(gamma_component, gamma_argument, precision);

    // For those of non-trivial value calculate the real and imaginary parts
    // of the arguments for the confluent hypergeometric function component
    // of the amplification factor
    double hyper_arg_a_imag = dimensionless_frequency/2;

    double source_position_sq = source_position * source_position;
    double hyper_arg_z_imag = hyper_arg_a_imag * source_position_sq;

    // Construct the arguments to the confluent hypergeometric function
    acb_t hyper_arg_a;
    acb_t hyper_arg_b;
    acb_t hyper_arg_z;

    acb_init(hyper_arg_a);
    acb_init(hyper_arg_b);
    acb_init(hyper_arg_z);

    acb_set_d_d(hyper_arg_a, 0, hyper_arg_a_imag);
    acb_one(hyper_arg_b);
    acb_set_d_d(hyper_arg_z, 0, hyper_arg_z_imag);

    // Calculate the confluent hypergeometric function value
    acb_t hyper_component;
    acb_init(hyper_component);

    acb_hypgeom_1f1(hyper_component,
                    hyper_arg_a,
                    hyper_arg_b,
                    hyper_arg_z,
                    0,
                    precision);

    // Construct the final value of the amplification factor by first
    // calculating the value of the exponential * gamma components and then
    // finally multiplying by the hyper component value
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

    // Clean up of declared acb_ts for memory management
    acb_clear(exponent);
    acb_clear(exponential_component);
    acb_clear(gamma_argument);
    acb_clear(gamma_component);
    acb_clear(hyper_arg_a);
    acb_clear(hyper_arg_b);
    acb_clear(hyper_arg_z);
    acb_clear(hyper_component);
    acb_clear(exponential_gamma);
}

// Function computes the mangification for the geometric optics approximation
// with the plus and minus images given by the state of mode. This is given by
// 1/2 +- (y^2 + 2)/(2y * sqrt(y^2 + 4)). Mode determines plus or minus
double Magnification(double source_position, int mode) {
    double magnification = 1./2.;
    double second_term_numerator = source_position*source_position + 2;
    double second_term_denominator = 2 * source_position * sqrt(
        source_position * source_position + 4);
    double second_term = second_term_numerator/second_term_denominator;

    if (mode == 1) {
        magnification += second_term;
    } else {
        magnification -= second_term;
    }

    return magnification;
}

// Function computes the time delay for the geometric optics approximation.
// This is given by:
// y * (sqrt(y^2 + 4)/2) + ln((sqrt(y^2+4)+y)/(sqrt(y^2+4)-y))
double TimeDelay(double source_position) {
    double first_term = source_position;
    double first_term_numerator = sqrt(source_position*source_position + 4);
    first_term *= first_term_numerator/2.;

    double second_term_numerator = first_term_numerator + source_position;
    double second_term_denominator = first_term_numerator - source_position;
    double second_term = log(second_term_numerator/second_term_denominator);

    double time_delay = first_term + second_term;

    return time_delay;
}

// Function computes the amplification factor for an axially symmetric point
// mass lens using the geometric optics approximation for given values of
// dimensionless frequency and source position
std::complex<double> AmplificationFactorGeometric(
    double dimensionless_frequency, double source_position) {
    // Using the complex_literals i operator
    using std::literals::complex_literals::operator""i;

    // First compute the magnifications
    double magnification_plus = Magnification(source_position, 1);
    double magnification_minus = Magnification(source_position, 0);

    // Compute the Time Delay
    double time_delay = TimeDelay(source_position);

    // Calculate the exponential term
    double exponent_imag  = dimensionless_frequency * time_delay;
    std::complex<double> exp_term = 1i * exponent_imag;
    exp_term = exp(exp_term);

    // Construct the final amplification factor
    double first_term = sqrt(abs(magnification_plus));
    std::complex<double> second_term =
        sqrt(abs(magnification_minus)) * exp_term;
    std::complex<double> geometric_factor = first_term - 1i * second_term;

    return geometric_factor;
}
