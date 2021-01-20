//Header File for pointlens.cpp - a set of functions to calculate amplification factor matrices from a given set of impact parameters and dimensionless frequencies

//Mick Wright
//DoC: 20/01/2021 

#ifndef POINTLENS_H
#define POINTLENS_H 

#include <cmath> 
#include <acb.h>
#include <acb_hypgeom.h>
#include <iostream> 
#include <fstream> 
#include <vector>

//Function computes the value of xm - the impact parameter divided by a length normalisation constant for the phase constant corresponding to the minimum time delay i.e. that of the image that travels the shortest path to the observer 
double xm(double y); 

//Function computes the value of phi - the phase constant used to obtain the minium time delay induced by the lensing
double phi(double y); 

//Function computes the amplification factor for an axially symmetric point mass lens for given values of dimensionless frequency and impact parameter 
void ampFacCalc(acb_t ampFac, double w, double y); 

//Function constructs two matrices containing the real and imaginary parts of the value of the amplification function based upon two vectors containing values of dimensionless frequency and impact parameter and returns these inside of a pair object
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> ampFacMatrices(std::vector<double> w, std::vector<double> y); 

#endif
