#include "pointlens.h"

double xm(double y){

double xm = (y + sqrt(y*y + 4))/2;
return xm; 

}

double phi(double y){

double phiNum = xm(y) - y; 
double phi = (phiNum*phiNum)/2 - log(xm(y));  
return phi; 

}

void ampFacCalc(acb_t ampFac, double w, double y){

double expReal = (M_PI*w)/4; 
double expImag = (w/2)*(log(w/2)-(2*phi(y))); 

acb_t exponent;
acb_t expComp; 
acb_init(exponent);
acb_init(expComp); 
acb_set_d_d(exponent, expReal, expImag); 
acb_exp(expComp, exponent, 50); 

double gammaReal = 1;
double gammaImag = -(w/2); 

acb_t gammaArg; 
acb_t gammaComp; 
acb_init(gammaArg);
acb_init(gammaComp); 
acb_set_d_d(gammaArg, gammaReal, gammaImag); 
acb_gamma(gammaComp, gammaArg, 50); 

double hyperArgAImag = w/2; 
double hyperArgZImag = hyperArgAImag*(y*y); 

acb_t hyperArgA;
acb_t hyperArgB;
acb_t hyperArgZ; 
acb_t hyperComp; 
acb_init(hyperArgA); 
acb_init(hyperArgB); 
acb_init(hyperArgZ); 
acb_init(hyperComp); 
acb_set_d_d(hyperArgA, 0, hyperArgAImag); 
acb_one(hyperArgB); 
acb_set_d_d(hyperArgZ, 0, hyperArgZImag); 
acb_hypgeom_1f1(hyperComp, hyperArgA, hyperArgB, hyperArgZ, 0, 100);

acb_t expGamma;
acb_init(expGamma);
acb_mul(expGamma, expComp, gammaComp, 50); 
acb_mul(ampFac, expGamma, hyperComp, 50); 
}

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> ampFacMatrices(std::vector<double> w, std::vector<double> y){

long ySize = y.size();
long wSize = w.size();

std::vector<std::vector<double>> ampFacReal(ySize, std::vector<double> (wSize)); 
std::vector<std::vector<double>> ampFacImag(ySize, std::vector<double> (wSize));

for (int i=0; i<ySize; i++){
for (int j=0; j<wSize; j++){ 

acb_t ampFac;
acb_init(ampFac); 
ampFacCalc(ampFac, w[j], y[i]); 

ampFacReal[i][j] = arf_get_d(arb_midref(acb_realref(ampFac)), ARF_RND_NEAR); 
ampFacImag[i][j] = arf_get_d(arb_midref(acb_imagref(ampFac)), ARF_RND_NEAR); 

}
}

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> ampFacMats(ampFacReal, ampFacImag); 

return ampFacMats;
}

