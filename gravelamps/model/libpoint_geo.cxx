// Function definitions for calculating the amplification factor for the
// isolated point mass lens model in the geometric optics regime
//
// Mick Wright 2021

#include "src/point.h"

// Function computes the magnification of the either the plus or minus image
// in the geometric optics approximation at the specified source position
// as determined by the state of mode.
// This magnification is given by:
// mu_{+-} = 1/2 +- (y^2 + 2)/(2y * sqrt(y^2 + 4))
//
// Input:
//      double source_position : dimensionless displacement from optical axis
//      int mode : determines whether to calculate plus or minus image based on
//                 whether 1 or 0 respectively
//
// Output:
//      double magnification : magnification of specified image as described
//                             above
double Magnification(double source_position, int mode) {
    double magnification = 1./2.;

    double source_position_sq = source_position * source_position;
    double source_position_numerator = source_position_sq + 2;
    double source_position_denominator = 2 * source_position *
                                         sqrt(source_position_sq + 4);
    double source_position_term = source_position_numerator/
                                  source_position_denominator;

    // Determine whether to add or subtract second term based on value of mode
    if (mode == 1) {
        magnification += source_position_term;
    } else {
        magnification -= source_position_term;
    }

    return magnification;
}

// Function computes the time delay for the geometric optics approximation at
// the specified source position. This is given by:
// y * (sqrt(y^2 + 4)/2) + ln((sqrt(y^2 + 4) + y)/sqrt(y^2 + 4) - y)
//
// Input:
//      double source_position : dimensionless displacement from optical axis
//
// Output:
//      double time_delay : dimensionless time delay induced by the lensing
double TimeDelay(double source_position) {
    double source_position_sq = source_position * source_position;
    double source_position_sqrt_term = sqrt(source_position_sq + 4.);

    double multiplication_term = source_position * source_position_sqrt_term/2;

    double log_term_numerator = source_position_sqrt_term + source_position;
    double log_term_denominator = source_position_sqrt_term - source_position;
    double log_term = log(log_term_numerator/log_term_denominator);

    double time_delay = multiplication_term + log_term;

    return time_delay;
}

// Function computes the amplification factor for an isolated axially symmetric
// point mass lens using the geometric optics approximation for given values of
// dimensionless frequency and source position.
//
// Input:
//      double dimensionless_frequency : dimensionless form of the freuqency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//
// Output:
//      std::complex<double> amplification_factor : complex value of the
//                                                  amplification factor
std::complex<double> AmplificationFactorGeometric(
    double dimensionless_frequency, double source_position) {
    // i Operator
    using std::literals::complex_literals::operator""i;

    // Compute the magnifications of both images
    double magnification_plus = Magnification(source_position, 1);
    double magnification_minus = Magnification(source_position, 0);

    // Compute magnification term
    double magnification_term = sqrt(abs(magnification_plus));

    // Compute the time delay
    double time_delay = TimeDelay(source_position);

    // Calculate the exponential term
    std::complex<double> exponential_term = 1i * dimensionless_frequency
                                            * time_delay;
    exponential_term = 1i * sqrt(abs(magnification_minus))
                       * exp(exponential_term);

    // Construct the final amplification factor
    std::complex<double> amplification_factor = magnification_term
                                                - exponential_term;

    return amplification_factor;
}
