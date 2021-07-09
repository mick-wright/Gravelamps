// Header File for pointlens.cc - a set of functions to calculate the
// amplification factor matrices from a given set of impact parameters and
// dimensionless frequencies
//
// Mick Wright 2021

#ifndef GRAVELAMPS_LENSING_POINTLENS_H_
#define GRAVELAMPS_LENSING_POINTLENS_H_

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <utility>
#include <string>
#include <complex>

#include "acb.h"
#include "acb_hypgeom.h"

// Function computes the value of xm - the impact parameter divided by a length
// normalisation constant for the phase constant corresponding to the minimum
// time delay i.e. that of the image that travels the shortest path to the
// observer
double StationaryPointMinimum(double impact_parameter);

// Function computes the value of phi - the phase constant used to obtain the
// minimum time delay induced by the lensing
double MinTimeDelayPhaseConstant(double impact_parameter);

// Function computes the amplification factor for an axially symmetric point
// mass lens with full wave optics for given values of dimensionless frequency
// and impact parameter with arithmetic precision given by precision
void AmplificationFactorCalculation(acb_t amplification_factor,
                                    double dimensionless_frequency,
                                    double impact_parameter,
                                    slong precision);

// Function computes the magnification for the geometric optics approximation
// with the plus and minus images given by the state of mode
double Magnification(double impact_parameter, int mode);

// Function computes the time delay for the geometric optics approximation
double TimeDelay(double impact_parameter);

// Function computes the amplification factor for an axially symmetric point
// mass lens using the geometric optics approximation for given values of
// dimensionless frequency and impact parameter
std::complex<double> AmplificationFactorGeometric(
    double dimensionless_frequency, double impact_parameter);

// Function constructs two matrices containing the real and imaginary parts of
// the value of the amplification factor function based upon two vectors
// containing values of dimensionless frequency and impact parameter and
// returns these inside of a pair object
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
    AmplificationFactorMatrices(std::vector<double> dimensionless_frequency,
                                std::vector<double> impact_parameter,
                                slong precision);

#endif  // GRAVELAMPS_LENSING_POINTLENS_H_
