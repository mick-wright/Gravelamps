// Header File for nfwlens.cc - a set of functions to calculate the
// amplification factor - split into real and imaginary matrices for a given
// set of source positions and dimensionless frequencies
//
// Mick Wright 2021

#ifndef GRAVELAMPS_LENSING_NFW_H_
#define GRAVELAMPS_LENSING_NFW_H_

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

#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_calc.h"

// Function computes the value of psi - the value of the lensing potential
void LensingPotential(acb_t lensing_potential,
                      const acb_t scaled_surface_density,
                      acb_t scaling_constant,
                      slong precision);

std::complex<double> LensingPotential(double scaled_surface_density,
                                      double scaling_constant);

// Function computes the first two terms of the time delay function for
// finding the minimum which yields the phase guage needed
std::complex<double> TimeDelayPartial(double scaled_surface_density,
                                      double source_position,
                                      double scaling_constant);

// Function computes the phase needed for a minimum time delay of zero
std::complex<double> MinTimeDelayPhase(double source_position,
                                       double scaling_constant);

// Function computes the lens equation value
double LensEquation(double scaled_surface_density,
                    double source_position,
                    double scaling_constant);

// Function computes the image position for a given source position by means
// of root finding the lens equation, looking for the negative root,
// the positive root, and a root near zero
std::vector<double> ImagePositions(double source_position,
                                   double scaling_constant);

// Function computes the time delay for a given image position and source
// position and minimum phase delay. This is the complete version of the
// above partial function
std::complex<double> TimeDelay(double image_position,
                               double source_position,
                               double scaling_constant,
                               double phase_minimum);

// Function computes the NFW mass for a given image position and scaling
// constant
std::complex<double> ImagePositionMass(double image_position,
                                       double scaling_constant);

// Function computes the surface density for a given image position and
// scaling constant
std::complex<double> SurfaceDensity(double image_position,
                                    double scaling_constant);

// Function computes the determinant of the matrix used to calculate the
// magnification
std::complex<double> MatrixDeterminant(double image_position,
                         double scaling_constant);

// Function computes the magnification for a given image position and scaling
// constant
std::complex<double> Magnification(double image_position, double scaling_constant);

// Function computes the intermediate function k(w,y,x,ks) for the
// amplification factor calculation. The function k is given by
// -iw*(exp[iw(y^2/2)]*J0(wy*sqrt(2x))*exp(-iw*psi(sqrt(2x),ks)))
void IntermediateFunctionCalculation(acb_t intermediate_function_value,
                                     acb_t dimensionless_frequency,
                                     acb_t source_position,
                                     const acb_t integration_parameter,
                                     acb_t scaling_constant,
                                     double minimum_phase,
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
                         acb_t source_position,
                         acb_t integration_upper_limit,
                         acb_t scaling_constant,
                         double minimum_phase,
                         slong precision);

// Function computes the value of the second correction term for the
// amplification factor
void SecondCorrectionTerm(acb_t second_correction_term,
                          acb_t dimensionless_frequency,
                          acb_t source_position,
                          acb_t integration_upper_limit,
                          acb_t scaling_constant,
                          double minimum_phase,
                          slong precision);

// Function computes the amplification factor for an axially symmetric Navarro,
// Frenk, and White (NFW) lens using full wave optics for given values of
// dimensionless frequency and source position with arithmetic precision given
// by precision. The infinite integral is approximated by calculating the
// finite integral with upper limit given by integration_upper_limit
void AmplificationFactorCalculation(acb_t amplification_factor,
                                    double dimensionless_frequency,
                                    double source_position,
                                    double scaling_constant,
                                    double integration_upper_limit,
                                    slong precision);

// Function computes the amplification factor for an axially symmetric Navarro,
// Frenk, White (NFW) lens using the geometric optics approximation for given
// values of dimensionless frequency and source position
std::complex<double> AmplificationFactorGeometric(
    double dimensionless_frequency,
    double source_position,
    double scaling_constant);

extern "C" {

// Wrapper function converting the result of MinTimeDelayPhase to a real value
// for compatibility with ctypes for python
double MinTimeDelayPhaseReal(double source_position,
                             double scaling_constant);

// Wrapper fiunction converting the vector caclulated by ImagePositions to an
// array of doubles for compatibility with ctypes for python
double* ImagePositionArray(double source_position,
                            double scaling_constant);

// Wrapper function converting the amplification factor result calculated by
// AmplificationFactorGeometric to a pair of floats for compatibility with
// ctypes for python
double* AFGRealOnly(double dimensionless_frequency,
                    double source_position,
                    double scaling_constant);

// Function destroys the object placed within, used for deallocating the memory
// used by the above function
void destroyObj(double* object);


// Simplified version of AmplificationFactorGeometric that takes in the image
// positions and min time delay phase to speed computation
double* SimpleAmpFac(double dimensionless_frequency,
                      double source_position,
                      double scaling_constant,
                      double * image_positions,
                      double min_time_delay_phase,
                      int number_of_images);
}

#endif  // GRAVELAMPS_LENSING_NFW_H_
