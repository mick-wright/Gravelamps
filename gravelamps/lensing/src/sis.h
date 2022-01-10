// Header File for sis.cc - a set of functions to calculate the amplification
// factor - split into real and imaginary matrices for a given set of source
// positions and dimensionless frequencies
//
// Mick Wright 2021

#ifndef GRAVELAMPS_LENSING_SIS_H_
#define GRAVELAMPS_LENSING_SIS_H_

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

#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_calc.h"

// Function computesthe amplification factor for an axially symmetric singular
// isothermal sphere (SIS) lens using full wave optics for given values of
// dimensionless frequency and source position. It does this using a summation
// method with the infinite sum approximated up to a given threshold and with
// a given arithmetic precision
void AmplificationFactorCalculation(acb_t amplification_factor,
                                    double dimensionless_frequency,
                                    double source_position,
                                    slong summation_upper_limit,
                                    slong precision);

// Function computes the amplification factor for an axially symmetric singular
// isothermal sphere (SIS) style lens using the geometric optics approximation
// for given values of dimensionless frequency and source position
std::complex<double> AmplificationFactorGeometric(
    double dimensionless_frequency, double source_position);

extern "C" {

// Wrapper function converting the amplification factor result calculated by
// AmplificationFactorGeometric to a pair of floats for compatibility with
// ctypes for python
double* AFGRealOnly(
    double dimensionless_frequency, double source_position);

// Function destroys the object placed within, used for deallocating the memory
// used by the above function
void destroyObj(double* object);
}

#endif  // GRAVELAMPS_LENSING_SIS_H_
