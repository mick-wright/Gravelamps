#include <cmath> 
#include <acb.h> 
#include <acb_hypgeom.h> 

double xm(double y){

double xm = (y + sqrt(y*y + 4))/2;
return xm; 

}

double phi(double y){

double phiNum = xm(y) - y; 
double phi = (phiNum*phiNum)/2 - log(xm(y));  
return phi; 

}

void ampFac(double w, double y){

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
acb_hypgeom_1f1(hyperComp, hyperArgA, hyperArgB, hyperArgZ, 0, 10000);

acb_t expGamma;
acb_t ampFac; 
acb_init(expGamma);
acb_init(ampFac); 
acb_mul(expGamma, expComp, gammaComp, 50); 
acb_mul(ampFac, expGamma, hyperComp, 50); 

acb_printn(expComp, 50, 0); flint_printf("\n"); 
acb_printn(gammaComp, 50, 0); flint_printf("\n"); 
acb_printn(hyperComp, 50, 0); flint_printf("\n"); 

acb_printn(ampFac, 50, 0); flint_printf("\n"); 
}

int main(){
double w = 132;
double y = 0.4;

ampFac(w,y); 
}

