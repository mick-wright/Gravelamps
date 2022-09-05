// Utility function definitions for assisting in the calculation of the
// amplification factor for the isolated singular isothermal sphere (SIS)
// lens model
//
// Mick Wright 2022

#include "src/sis.h"

// Wrapper converts the std::complex value of the amplification factor
// calculated by AmplificationFactorGeometric into a pair of doubles for
// compatibility with ctypes.
//
// Input:
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//
// Ouptut:
//      double* amplification_array : array length 2 containing the real and
//                                    imaginary parts of the amplification
//                                    factor value
double* PyAmplificationFactorGeometric(
    double dimensionless_frequency, double source_position) {

    // Calculate complex value using AmplificationFactorGeometric
    std::complex<double> amplification_factor =\
        AmplificationFactorGeometric(dimensionless_frequency, source_position);

    double* amplification_array = new double[2];
    amplification_array[0] = std::real(amplification_factor);
    amplification_array[1] = std::imag(amplification_factor);

    return amplification_array;
}
