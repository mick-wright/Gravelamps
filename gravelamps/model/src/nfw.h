// Header file for the functions relating to calculating the amplification
// factor for the isolated Navarro, Frenk, White (NFW) lensing profile in both
// the wave and geometric optics regimes as well as uitility functions
//
// Mick Wright 2021

#ifndef GRAVELAMPS_MODEL_SRC_NFW_H_
#define GRAVELAMPS_MODEL_SRC_NFW_H_

#include <acb.h>
#include <acb_hypgeom.h>
#include <acb_calc.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <utility>
#include <complex>
#include <string>
#include <functional>
#include <sstream>

#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>

#include "./utils.h"

// Wave Optics Functions : libnfw_wave

// Function computes the value of the lensing potential for specified values
// of the image position and scaling constant. This value is zero in the case
// where the image position is equal to one.
std::complex<double> LensingPotential(double image_position,
                                      double scaling_constant);

// Function computes the time delay for given image and source positions
// at the given scale modified by a specified phase.
std::complex<double> TimeDelay(double image_position,
                               double source_position,
                               double scaling_constant,
                               double phase);

// Function computes the phase required such that the minimum time delay
// induced by the lensing is zero for the given scaling of the profile at the
// specified source position.
std::complex<double> Phase(double source_position, double scaling_constant);

// Function computes the value of the positive exponential term for the
// intermediate K function used in the calculation of the amplification factor
// for the specified dimensionless frequency, source position, and phase with
// the exponentiation being done at the specified precision. The resultant
// value is assigned to the given acb_t object.
void KFunctionPositiveExponential(acb_t positive_exponential_term,
                                  double dimensionless_frequency,
                                  double source_position,
                                  double phase,
                                  slong precision);

// Function computes the value of the negative exponential term for the
// intermediate K function used in the claculation of the amplification factor
// for the specified dimensionless frequency, image position, and scaling
// constant with the exponentiation being done at the specified precision.
// The resultant value is assigned to the given acb_t object.
void KFunctionNegativeExponential(acb_t negative_exponential_term,
                                  double dimensionless_frequency,
                                  double image_position,
                                  double scaling_constant,
                                  slong precision);

// Function computes the value of the bessel function term for the intermediate
// K function used in the calculation of the amplification factor for the
// specified dimensionless frequency, source position, and image position with
// the function claculated to the specified precision. The resultant value is
// assigned to the given acb_t object.
void KFunctionBessel(acb_t bessel_term,
                     double dimensionless_frequency,
                     double source_position,
                     double image_position,
                     slong precision);

// Function computes the value of the intermediate K function used in the
// calculation of the amplification factor for the wave optics case for the
// specified dimensionless frequency, source position, image position,
// scaling constant, and phase wiht the functions calculated wiht specified
// precision. The value of the function is assigned to the given acb_t object.
void KFunctionCalculation(acb_t k_function_value,
                          double dimensionless_frequency,
                          double source_position,
                          double image_position,
                          double scaling_constant,
                          double phase,
                          slong precision);

// Function computes the value of the integrand for the integral in the
// calculation of the amplification factor for the specified parameters.
// These parameters are passed as a single object to the function and then
// extracted for consistency with arb's integration methodology. The order
// parameter required by said methodology is not implemented, however,
// the functions are calculated with precision as specified. The resultant
// value is assigned to the pointer provided to the function.
int AmplificationWaveIntegrand(acb_ptr integrand,
                               const acb_t integration_parameter,
                               void * parameter_set,
                               slong order,
                               slong precision);

// Function computs the value of the exponential term of the first correction
// to the amplification factor for the specified value of dimensionless
// frequency and maximum image position at the specified precision. The
// resultant value is assigned to the given acb_t object.
void FirstCorrectionExponential(acb_t exponential_term,
                                double dimensionless_frequency,
                                double max_image_position,
                                slong precision);

// Function computes the value of the first correction term to the
// amplification factor for the specified values of the dimensionless
// frequency, source position, scaling constant, and phase at the specifed
// precision. This correction term relies on the maximal image position value
// being considered in the integration. The resultant value is assigned to the
// given acb_t object.
void AmplificationWaveFirstCorrection(acb_t first_correction_term,
                                      double dimensionless_frequency,
                                      double source_position,
                                      double max_image_position,
                                      double scaling_constant,
                                      double phase,
                                      slong precision);

// Function computes the value of the derivative term of the second correction
// term to the amplification factor at the specified dimensionless frequency,
// source position, scaling constant, image position, and phase values. The
// derivative is calculated using a central finite differences method. The
// resultant value is assigned to the given acb_t object.
void SecondCorrectionDerivative(acb_t derivative_term,
                                double dimensionless_frequency,
                                double source_position,
                                double max_image_position,
                                double scaling_constant,
                                double phase,
                                slong precision);

// Function computes the value of the exponential term of the second correction
// term to the amplification factor at the specified values of dimensionless
// frequency and image position with the desired precision in the
// exponentiation. The resultant value is assigned to the given acb_t object.
void SecondCorrectionExponential(acb_t exponential_term,
                                 double dimensionless_frequency,
                                 double max_image_position,
                                 slong precision);

// Function computes the value of the second correction term to the
// amplification factor at the specified dimensionless frequency, source
// position, scaling constant, and phase values. The correction term also
// relies on the maximal image position being considered as part of the
// integration in the main amplification factor calculation. The resultant
// value is assigned to the given acb_t object.
void AmplificationWaveSecondCorrection(acb_t second_correction_term,
                                       double dimensionless_frequency,
                                       double source_position,
                                       double max_image_position,
                                       double scaling_constant,
                                       double phase,
                                       slong precision);

// Function computes the amplification factor for an isolated axially symmetric
// Navarro, Frenk, White (NFW) lens for given values of dimensionless frequency,
// source position, and scaling constant with functions calculated with the
// specified precision. The infinite integral is approximated as a finite
// integral up to the specified value. The resultant value is assigned to the
// given amplification factor object.
void AmplificationFactorWave(acb_t amplification_factor,
                             double dimensionless_frequency,
                             double source_position,
                             double scaling_constant,
                             double max_image_position,
                             slong precision);

// Geometric Optics Function : libnfw_geo

// Function computes the value of the lens equation for a given set of image
// position, source position, and scaling constant of the profile.
double LensEquation(double image_position,
                    double source_position,
                    double scaling_constant);

// Function computes the image positions for a given source position by means
// of root finding the lens equation in the negative, zero, and positive cases.
std::vector<double> ImagePositions(double source_position,
                                   double scaling_constant);

// Function computes the NFW mass associated with a given image position and
// scaling constant.
std::complex<double> ImageMass(double image_position,
                               double scaling_constant);

// Function computes the surface density associated with an image at the
// specified position for the given scaling.
std::complex<double> SurfaceDensity(double image_position,
                                    double scaling_constant);

// Function computes the magnification of an image at the specified position
// for the given scaling
std::complex<double> Magnification(double image_position,
                                   double scaling_constant);

// Function computes the value of the morse factor for an image based on the
// dimensionless position and scaling constant. This value can take one of
// three distinct forms based on the sign of the surface density.
double MorseFactor(double image_position, double scaling_constant);

// Function computes the contribution of a single image to the amplification
// factor in the geometric optics approximation with the specified dimnensionless
// frequency, source position, scaling constant, and phase.
std::complex<double> ImageContribution(double image_position,
                                       double dimensionless_frequency,
                                       double source_position,
                                       double scaling_constant,
                                       double phase);

// Function computes the amplification factor for an isolated axially symmetric
// Navarro, Frenk, White (NFW) at a given scaling using the geometric optics
// appxoimation for given values of dimensionless frequency and source
// position. Optionally one may give the image positions and phase if these
// have been precomputed to optimize for speed.
std::complex<double> AmplificationFactorGeometric(
    double dimensionless_frequency,
    double source_position,
    double scaling_constant,
    double * image_array,
    double phase,
    int number_of_images);

// Utility Functions : libnfw_utils

// C Functions that are used by the python side of Gravelamps. These are
// written so as to be compatible with the python builtin ctypes
extern "C" {
    // Wrapper function returning the real component of the phase required for
    // the minimum time delay induced by the lensing to be zero. This is for
    // compatibility with ctypes for python.
    double PyPhase(double source_position, double scaling_constant);

    // Wrapper function converts the vector calculated by ImagePositions to
    // an array of doubles with the number of images inserted into the start
    // of the array.
    double* PyImagePositions(double source_position, double scaling_constant);

    // Wrapper converts the std::complex value of the amplification factor
    // calculated by AmplificationFactorGeometric into a pair of doubles for
    // compatibility with ctypes
    double* PyAmplificationFactorGeometric(double dimensionless_frequency,
                                           double source_position,
                                           double scaling_constant,
                                           double* image_array,
                                           double phase,
                                           int number_of_images);
}

// Interpolator Generator Functions : libnfw_generator

// Function takes in values of dimensionless frequency and source position and
// calculates the amplification factor for the model with given scaling using
// either wave or geometric optics depending upon whether the dimensionless
// frequency is above the specified switch. It returns in either case a pair
// of doubles contianing the real and imaginary components.
double* AmplificationFactor(double dimensionless_frequency,
                            double source_position,
                            double scaling_constant,
                            double max_image_position,
                            slong precision,
                            slong geo_switch,
                            double* image_positions,
                            double phase,
                            int number_of_images);


// Function takes in a pair of files containing a grid of values of
// dimensionless frequency and source position and generates a pair of files
// containing the corresponding grid of the real and imaginary parts of the
// amplification factor. If these files exist already it will read them in
// and use them to avoid repeating calculations, allowing for the process to be
// interrupted.
extern "C" {
    int GenerateLensData(char* dimensionless_frequency_file,
                         char* source_position_file,
                         char* amplification_factor_real_file,
                         char* amplification_factor_imag_file,
                         double scaling_constant,
                         double max_image_position,
                         slong precision,
                         slong geo_switch);
}

#endif  // GRAVELAMPS_MODEL_SRC_NFW_H_
