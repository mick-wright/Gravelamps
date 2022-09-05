// Function definitions for calculating the amplification factor for the
// isolated singular isothermal sphere (SIS) lens model in the geometric optics
// approximation regime.
//
// Mick Wright 2021

#include "src/sis.h"

// Function computes the magnification of either the plus or minus image
// in the geometric optics approximation at the specified source position as
// determined by the state of mode. This magnification is given by
// mu_{+-} = +- 1 + 1/source_position
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
    double inverse_source_position = 1./source_position;

    double magnification;
    if (mode == 1) {
        magnification = 1 + inverse_source_position;
    } else {
        magnification = inverse_source_position - 1;
    }

    return magnification;
}

// Function computes the time delay for the geometric optics approximation at
// the specified source position. This is given by:
// td = 2y
//
// Input:
//      double source_position : dimensionless displacement from optical axis
//
// Output:
//      double time_delay : time delay induced by the lensing
double TimeDelay(double source_position) {
    double time_delay = 2 * source_position;
    return time_delay;
}

// Function computes the amplification factor for an isolated axially symmetric
// SIS lens using the geometric optics approximation for given values of
// dimensionless frequency and source position.
//
// Input:
//      double dimensionless_frequency : dimensionless form of the frequency
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

    // Calculate the mangification term
    double magnification_term = sqrt(abs(Magnification(source_position, 1)));

    if (source_position >= 1) {
        std::complex<double> amplification_factor = magnification_term;
        return amplification_factor;
    }

    // Compute the exponential term
    double magnification_minus = sqrt(abs(Magnification(source_position, 0)));
    double time_delay = TimeDelay(source_position);

    std::complex<double> exponential_term = 1i * dimensionless_frequency\
                                            * time_delay;
    exponential_term = 1i * magnification_minus * exp(exponential_term);

    // Construct the final amplification factor
    std::complex<double> amplification_factor = magnification_term\
                                                - exponential_term;

    return amplification_factor;
}
