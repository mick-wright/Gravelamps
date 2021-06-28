// Header File for sis.cc - a set of functions to calculate the amplification
// factor - split into real and imaginary matrices for a given set of impact
// parameters and dimensionless frequencies
//
// Mick Wright 2021

#ifndef GWLENSING_LENSING_SIS_H_
#define GWLENSING_LENSING_SIS_H_

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <utility>
#include <string>

#include "acb.h"
#include "acb_hypgeom.h"

// Function computes the value of phi - the phase constant used to obtain the
// minimum time delay induced by the lensing
double MinTimeDelayPhaseConstant(double impact_parameter);

// Function computesthe amplification factor for an axially symmetric singular
// isothermal sphere (SIS) lens for given values of dimensionless frequency and
// impact parameter with arithmetic precision given by precision. The infinte
// sum is calculated to an upper value given by sum_threshold
void AmplificationFactorCalculation(acb_t amplification_factor,
                                    double dimensionless_frequency,
                                    double impact_parameter,
                                    int sum_threshold,
                                    slong precision);

// Function constructs two matrices containing the real and imaginary parts of
// the value of the amplification factor function based upon two vectors
// containing values of dimensionless frequency and impact parameter and
// returns these inside of a pair object
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
    AmplificationFactorMatrices(std::vector<double> dimensionless_frequency,
                                std::vector<double> impact_parameter,
                                int sum_threshold,
                                slong precision);

#endif  // GWLENSING_LENSING_SIS_H_
