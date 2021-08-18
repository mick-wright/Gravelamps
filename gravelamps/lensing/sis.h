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

#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_calc.h"

// Function computes the value of the intermediate function k(w,y,z) for the
// amplification factor calculation.
void IntermediateFunctionCalculation(acb_t intermediate_function_value,
                                     acb_t dimensionless_frequency,
                                     acb_t source_position,
                                     const acb_t integration_parameter,
                                     slong precision);

// Function computes teh value of the integrand being integrated in the
// amplification factor calculation
int SisIntegrand(acb_ptr integrand,
                 const acb_t integration_parameter,
                 void * parameter_set,
                 slong order,
                 slong precision);

// Function computes the value of the first correction term for the
// amplification factor
void FirstCorrectionTerm(acb_t first_correction_term,
                         acb_t dimensionless_frequency,
                         acb_t source_position,
                         acb_t integration_upper_limit,
                         slong precision);

// Function computes the value of the second correction term for the
// amplification factor calculation
void SecondCorrectionTerm(acb_t second_correction_term,
                          acb_t dimensionless_frequency,
                          acb_t source_position,
                          acb_t integration_upper_limit,
                          slong precision);

// Function computesthe amplification factor for an axially symmetric singular
// isothermal sphere (SIS) lens using full wave optics for given values of
// dimensionless frequency and source position with arithmetic precision given
// by precision. The infinte integral is calculated up to a finite limit given
// by integration_upper_limit.
void AmplificationFactorCalculation(acb_t amplification_factor,
                                    double dimensionless_frequency,
                                    double source_position,
                                    double integration_upper_limit,
                                    slong precision);

// Function computes the amplification factor for an axially symmetric singular
// isothermal sphere (SIS) style lens using the geometric optics approximation
// for given values of dimensionless frequency and source position
std::complex<double> AmplificationFactorGeometric(
    double dimensionless_frequency, double source_position);

// Function constructs two matrices containing the real and imaginary parts of
// the value of the amplification factor function based upon two vectors
// containing values of dimensionless frequency and source position and
// returns these inside of a pair object
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
    AmplificationFactorMatrices(std::vector<double> dimensionless_frequency,
                                std::vector<double> source_position,
                                double integration_upper_limit,
                                slong precision,
                                slong approx_switch);

#endif  // GRAVELAMPS_LENSING_SIS_H_
