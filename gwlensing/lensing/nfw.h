// Header File for nfwlens.cc - a set of functions to calculate the
// amplification factor - split into real and imaginary matrices for a given
// set of impact parameters and dimensionless frequencies
//
// Mick Wright 2021

#ifndef GWLENSING_LENSING_NFW_H_
#define GWLENSING_LENSING_NFW_H_

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <utility>
#include <complex>
#include <string>

#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_calc.h"

// Function computes the value of psi - the value of the lensing potential
void LensingPotential(acb_t lensing_potential,
                      const acb_t scaled_surface_density,
                      acb_t scaling_constant,
                      slong precision);

// Function computes the intermediate function k(w,y,x,ks) for the
// amplification factor calculation. The function k is given by
// -iw*(exp[iw(y^2/2)]*J0(wy*sqrt(2x))*exp(-iw*psi(sqrt(2x),ks)))
void IntermediateFunctionCalculation(acb_t intermediate_function_value,
                                     acb_t dimensionless_frequency,
                                     acb_t impact_parameter,
                                     const acb_t integration_parameter,
                                     acb_t scaling_constant,
                                     slong precision);

// Function computes the value of the integrand being integrated in the
// amplification factor calculation. The order parameter is unused but is
// required by the integration process
int NfwIntegrand(acb_ptr integrand,
                 const acb_t integration_parameter,
                 void * parameter_set,
                 slong order,
                 slong precision);

// Function computes the value of the first correction term for the
// amplification factor.
void FirstCorrectionTerm(acb_t first_correction_term,
                         acb_t dimensionless_frequency,
                         acb_t impact_parameter,
                         acb_t integration_upper_limit,
                         acb_t scaling_constant,
                         slong precision);

// Function computes the value of the second correction term for the
// amplification factor
void SecondCorrectionTerm(acb_t second_correction_term,
                          acb_t dimensionless_frequency,
                          acb_t impact_parameter,
                          acb_t integration_upper_limit,
                          acb_t scaling_constant,
                          slong precision);

// Function computes the amplification factor for an axially symmetric Navarro,
// Frenk, and White (NFW) lens for given values of dimensionless frequency and
// impact parameter with arithmetic precision given by precision. The infinite
// integral is approximated by calculating the finite integral with upper limit
// given by integration_upper_limit
void AmplificationFactorCalculation(acb_t amplification_factor,
                                    double dimensionless_frequency,
                                    double impact_parameter,
                                    double scaling_constant,
                                    double integration_upper_limit,
                                    slong precision);

// Function constructs two matrices containing the real and imaginary parts of
// the value of the amplification factor function based upon two vectors
// containing values of dimensionless frequency and impact parameter and
// returns these inside of a pair object
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
    AmplificationFactorMatrices(std::vector<double> dimensionless_frequency,
                                std::vector<double> impact_parameter,
                                double scaling_constant,
                                double integration_upper_limit,
                                slong precision);

#endif  // GWLENSING_LENSING_NFW_H_
