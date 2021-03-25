#ifndef SIS_H
#define SIS_H

#include <cmath>
#include "acb.h" 
#include "acb_hypgeom.h" 
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm> 

double phi(double y);

void ampFacCalc(acb_t ampFac, double w, double y); 

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> ampFacMatrices(std::vector<double> w, std::vector<double> y);

#endif 
