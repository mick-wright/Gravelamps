// Header file for the functions relating to calculating the amplification
// factor for the isolated point mass lensing profile in both the wave and
// geometric optics regimes as well as utility functions.
//
// Mick Wright 2021

#ifndef GRAVELAMPS_MODEL_SRC_POINT_H_
#define GRAVELAMPS_MODEL_SRC_POINT_H_

#include <acb.h>
#include <acb_hypgeom.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <utility>
#include <string>
#include <sstream>
#include <complex>
#include <iomanip>

#include "./utils.h"

// Wave Optics Functions : libpoint_wave

// Function computes the value of the displacement as a function of the
// dimensionless source position used in the calculation of the phase required
// for zero minimal time delay.
double StationaryPointMinimum(double source_position);

// Function computes the value of the phase required for the minimum time delay
// induced by the lensing to be zero.
double Phase(double source_position);

// Function computes the exponential component of the amplification factor for
// the wave optics case at the dimensionless frequency and source position
// with the exponentiation occuring with the specified precision. The resultant
// value is assigned to the given acb_t object.
void ExponentialComponent(acb_t exponential_component,
                          double dimensionless_frequency,
                          double source_position,
                          slong precision);

// Function computes the gamma function component of the amplification factor
// for the wave optics case at the specified dimensionless frequency with the
// specified precision. The resultant value is assigned to the given acb_t
// object.
void GammaComponent(acb_t gamma_component,
                    double dimensionless_frequency,
                    slong precision);

// Function computes the confluent hypergeometric component of the
// amplification factor for the wave optics case at the specified dimensionless
// frequency and source position with the specified precision. The resultant
// value is assigned to the given acb_t object.
void HyperComponent(acb_t hyper_component,
                    double dimensionless_frequency,
                    double source_position,
                    slong precision);

// Function computes the amplification factor for an isolated axially symmetric
// point mass lens using full wave optics for given values of dimensionless
// frequency and source position with the given arithmetic precision. This
// value is given to the provided acb_t object.
void AmplificationFactorWave(acb_t amplification_factor,
                             double dimensionless_frequency,
                             double source_position,
                             slong precision);

// Geometric Optics Functions : libpoint_geo

// Function computes the magnification of either the plus or minus image in the
// geometric optics approximation at the specified source position as
// determined by the state of mode.
double Magnification(double source_position, int mode);

// Function computes the time delay for the geometric optics approximation at
// the specified source position.
double TimeDelay(double source_position);

// Function computes the amplification factor for an axially symmetric point
// mass lens using the geometric optics approximation for given values of
// dimensionless frequency and source position
std::complex<double> AmplificationFactorGeometric(
    double dimensionless_frequency, double source_position);


// Utility Functions : libpoint_utils

// C Functions that are used by the python side of Gravelamps. These are
// written so as to be compatible with the python builtin ctypes.
extern "C" {
    // Wrapper function converting the std::complex value of the amplification
    // factor calculated by AmplificationFactorGeometric into a pair of doubles
    // for compatibility with ctypes.
    double* PyAmplificationFactorGeometric(
        double dimensionless_frequency, double source_position);
}

// Generator Functions : libpoint_generator

// Function takes in values of dimensionless frequency and source position and
// calculates the amplification factor for the model using either wave or
// geometric optics depending upon wehther the dimensionless frequency is above
// the specified switch. It returns in either case a pair of doubles containing
// the real and imaginary components.
double* AmplificationFactor(double dimensionless_frequency,
                            double source_position,
                            slong precision,
                            slong geo_switch);

// Function takes in a pair of files containing a grid of values of
// dimensionless frequency and source position and generates a pair of files
// contianing the correspoinding grid of the real and imaginary parts of the
// amplification factor. If these files exist already it will read them in and
// use them to avoid repeating calculations, allowing for th eprocess to be
// interrupted.

extern "C" {
    int GenerateLensData(char* dimensionless_frequency_file,
                         char* source_position_file,
                         char* amplification_factor_real_file,
                         char* amplification_factor_imag_file,
                         int precision_int,
                         int geo_switch_int);
}

#endif  // GRAVELAMPS_MODEL_SRC_POINT_H_
