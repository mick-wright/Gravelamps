#include <cmath> 
#include <acb.h> 
#include <acb_hypgeom.h> 
#include <iostream>
#include <fstream> 
#include <vector> 
#include <iterator>
#include <chrono>

double xm(double y){

double xm = (y + sqrt(y*y + 4))/2;
return xm; 

}

double phi(double y){

double phiNum = xm(y) - y; 
double phi = (phiNum*phiNum)/2 - log(xm(y));  
return phi; 

}

void ampFac(FILE* outfile, double w, double y){

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
acb_t ampFac; 
acb_init(expGamma);
acb_init(ampFac); 
acb_mul(expGamma, expComp, gammaComp, 50); 
acb_mul(ampFac, expGamma, hyperComp, 50); 

acb_fprintn(outfile, ampFac, 20, 0); flint_fprintf(outfile, "\t"); 
}

int main(){

std::ifstream wfile("w.dat");
std::istream_iterator<double> wstart(wfile), wend; 
std::vector<double> wVals(wstart, wend); 

std::ifstream yfile("y.dat"); 
std::istream_iterator<double> ystart(yfile), yend; 
std::vector<double> yVals(ystart, yend); 

std::cout << wVals.size() << std::endl; 
std::cout << yVals.size() << std::endl;

FILE * outfile; 
outfile = fopen("c++_f.dat", "w"); 

auto t1 = std::chrono::high_resolution_clock::now(); 
for (auto i:yVals){
for (auto j:wVals){
ampFac(outfile,i,j); 
}
flint_fprintf(outfile, "\n"); 
}
auto t2 = std::chrono::high_resolution_clock::now(); 

std::cout << "This took " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " milliseconds" << std::endl;

}
