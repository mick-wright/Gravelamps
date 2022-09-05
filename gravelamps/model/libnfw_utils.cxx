// Utility function definitions for assisting in the calculation of the
// amplification factor for the isolated Navarro, Frenk, White (NFW) lens
// model.
//
// Mick Wright 2022

#include "src/nfw.h"

// Wrapper function returning the real compononent of the phase required for
// the minimum time delay induced by the lensing to be zero. This is for
// compatibility with ctypes for python.
//
// Input:
//      double source_position : dimensionless displacement from optical axis
//      double scaling_constant : characteristic scale of the profile
//
// Output:
//      double real_phase : real component of the phase as calculated by the
//                          function Phase
double PyPhase(double source_position, double scaling_constant) {
    std::complex<double> phase = Phase(source_position, scaling_constant);
    double real_phase = std::real(phase);

    return real_phase;
}

// Wrapper function converts the vector calculated by ImagePositions to
// an array of doubles with the number of images inserted into the start
// of the array.
//
// Input:
//      double source_position : dimensionless displacement from optical axis
//      double scaling_constant : characteristic scale of the profile
//
// Output:
//      double* image_array : array containing the image positions as
//                            calculated by the ImagePositions function
//                            preceded by the number of images
double* PyImagePositions(double source_position, double scaling_constant) {
    std::vector<double> image_positions = ImagePositions(source_position,
                                                         scaling_constant);

    int number_of_images = image_positions.size() + 1;
    image_positions.insert(image_positions.begin(), number_of_images);

    double* image_array = new double[image_positions.size()];
    std::copy(image_positions.begin(), image_positions.end(), image_array);

    return image_array;
}

// Wrapper converts the std::complex value of the amplification factor
// calculated by AmplificationFactorGeometric into a pair of doubles for
// compatibility with ctypes.
//
// Input:
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//      double scaling_constant : characteristic scale of the profile
//
// Optional Input:
//      double* image positions : optional array of image positions
//      double phase : value of the phase required for the minimal time delay
//                     induced by the lensing to be zero
//      int number_of_images : length of the image_positions array
//
// Output:
//      double* amplification_array : array containing the real and imaginary
//                                    parts of the amplification factor value
double* PyAmplificationFactorGeometric(double dimensionless_frequency,
                                       double source_position,
                                       double scaling_constant,
                                       double * image_array = NULL,
                                       double phase = 0.0,
                                       int number_of_images = 0) {
    // Calculate complex value using AmplificationFactorGeometric
    std::complex<double> amplification_factor =\
        AmplificationFactorGeometric(dimensionless_frequency,
                                     source_position,
                                     scaling_constant,
                                     image_array,
                                     phase,
                                     number_of_images);

    double* amplification_array = new double[2];
    amplification_array[0] = std::real(amplification_factor);
    amplification_array[1] = std::imag(amplification_factor);

    return amplification_array;
}
