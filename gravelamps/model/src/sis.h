// Header file for the functions relating to calculating the amplification
// factor for the isolated singular isothermal sphere (SIS) lensing profile in
// both the wave and geometric optics regimes as well as utility functions.
//
// Mick Wright 2021

#ifndef GRAVELAMPS_MODEL_SRC_SIS_H_
#define GRAVELAMPS_MODEL_SRC_SIS_H_

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
#include <string>
#include <complex>
#include <sstream>
#include <iomanip>

#include "./utils.h"

// Wave Optics Functions : libsis_wave

// Function computes the value of the phase required for the minimum time delay
// induced by the lensing to be zero.
double Phase(double source_position);

// Function computes the prefactor component for the amplification factor for
// the wave optics case at the dimensionless frequency and source position
// specified with the exponentiation occuring with the specified precision.
// The resultant value is assigned to the given acb_t object.
void PrefactorCalculation(acb_t prefactor,
                          double dimensionless_frequency,
                          double source_position,
                          slong precision);

// Function computes the gamma function component of the nth summation term for
// the given value of n with the desired precision. The resultant value is
// assigned to the given acb_t object.
void GammaComponent(acb_t gamma_term,
                    int n,
                    slong precision);

// Function computes the power component of the nth summation term for the
// given value of dimensionless frequency and n with the desired precision. The
// resultant value is assigned to the given acb_t object.
void PowerComponent(acb_t power_term,
                    double dimensionless_frequency,
                    int n,
                    slong precision);

// Function computes the confluent hypergeometric function component of the nth
// summation term for the given dimensionless frequency, source position, and n
// with the desired precision. The resultant value is assigned to the given
// acb_t object.
void HyperComponent(acb_t hyper_term,
                    double dimensionless_frequency,
                    double source_position,
                    int n,
                    slong precision);

// Function computes the nth summation term for the amplification factor for
// the wave optics case at the dimensionless frequency and source position
// specified with the desired precision. The resultant value is assigned to the
// given acb_t object.
void SummationNTerm(acb_t summation_term,
                    double dimensionless_frequency,
                    double source_position,
                    slong precision,
                    int n);

// Function computes the amplification factor for an isolated axially symmetric
// SIS lens using full wave optics for given values of dimensionless frequency
// and source position with the given arithmetic precision. The infinite sum
// in the analytical expression is approximated up to a specified term. The
// resultant value is given to the provided amplification factor object.
void AmplificationFactorWave(acb_t amplification_factor,
                             double dimensionless_frequency,
                             double source_position,
                             slong summation_upper_limit,
                             slong precision);

// Geometric Optics Functions : libsis_geo

// Function computes the magnification of either the plus or minus image in the
// geometric optics approximation at the sepcified source position as
// determined by the state of mode.
double Magnification(double source_position, int mode);

// Function computes the time delay for the geometric optics approximation at
// the specified source position
double TimeDelay(double source_position);

// Function computes the amplification factor for an isolated axially symmetric
// SIS lens using the geometric optics approximation for given values of
// dimensionless frequency and source position
std::complex<double> AmplificationFactorGeometric(
    double dimensionless_frequency, double source_position);

// Utility Functions : libsis_utils

// C Functions that are used by the python side of Gravelamps. These are
// written so as to be compatible with the python builtin ctypes.
extern "C" {
    // Wrapper function converting the std::complex value of the amplification
    // factor calculated by AmplificationFactorGeometric into a pair of doubles
    // for compatibility with ctypes.
    double* PyAmplificationFactorGeometric(
        double dimensionless_frequency, double source_position);
}

// Generator Functions : libsis_generator

// Function takes in values of dimensionless frequency and source position and
// calculates the amplification factor for the model using either wave or
// geometric optics depending upon whether the dimensionless frequency is above
// the specified switch. It returns in either case a pair of doubles containing
// the real and imaginary components
double* AmplificationFactor(double dimensionless_frequency,
                            double source_position,
                            slong summation_upper_limit,
                            slong precision,
                            slong geo_switch);

// Function takes in a pair of files containing a grid of values of
// dimensionless frequency and source position and generates a pair of files
// containing the corresponding grid of the real and imaginary parts of the
// amplification factor. If these files exist already it will read them in and
// use them to avoid repeating calculations, allowing for the process to be
// inerrupted.
extern "C" {
    int GenerateLensData(char* dimensionless_frequency_file,
                         char* source_position_file,
                         char* amplification_factor_real_file,
                         char* amplification_factor_imag_file,
                         int summation_upper_limit_int,
                         int precision_int,
                         int geo_switch_int);
}

#endif  // GRAVELAMPS_MODEL_SRC_SIS_H_
