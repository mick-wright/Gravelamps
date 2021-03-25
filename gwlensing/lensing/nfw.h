#ifndef NFW_H
#define NFW_H

#include <cmath>
#include "acb.h" 
#include "acb_hypgeom.h" 
#include "acb_calc.h" 
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm> 
#include <complex>
#include <boost/math/differentiation/finite_difference.hpp> 

//Function to calculate the lensing potential 
void acb_psi(acb_t psi, const acb_t x, acb_t ks, slong prec); 

//Intermediate function in calculating the amplification factor 
void acb_k_calc(acb_t k, acb_t w, acb_t y, const acb_t x, acb_t ks, slong prec); 

//std::complex version of k - used for differentation
std::complex<double> k_complex(double wVal, double yVal, double xVal, double ksVal); 

//Calculator for the Integrand in calculating the amplification factor
int f_nfwIntegrand(acb_ptr res, const acb_t x, void * param, slong order, slong prec); 

//Calculator for Correction Term One in calculating the amplification factor
void acb_f2t1(acb_t termOne, acb_t w, acb_t y, acb_t upperLimitX, acb_t ks, slong prec); 

//Calculator for Term Two in calculating the amplification factor
void acb_f2t2(acb_t termTwo, acb_t w, acb_t y, acb_t upperLimitX, acb_t ks, slong prec); 

//Function to Calculate the amplification factor
void ampFacCalc(acb_t ampFac, double wVal, double yVal, double ksVal, double upperLimitVal, slong prec); 

//Function to Construct the Amplification Factor Files
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> ampFacMatrices(std::vector<double> w, std::vector<double> y, double ksVal, double upperLimitVal, slong prec);

#endif 
