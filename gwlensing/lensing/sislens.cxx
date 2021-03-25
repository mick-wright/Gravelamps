#include "sis.h" 

double phi(double y){
return y+1/2.; 
}

void ampFacCalc(acb_t ampFac, double w, double y, int sumThreshold){

double prefactorExpImag = (w/2)*(y*y + 2*(phi(y))); 
acb_t prefactorExp; 
acb_init(prefactorExp); 
acb_set_d_d(prefactorExp, 0, prefactorExpImag); 
acb_exp(prefactorExp, prefactorExp, 500); 

acb_t sumTerm; 
acb_init(sumTerm); 
acb_zero(sumTerm); 

for (int n=0; n<=sumThreshold; n++){ 

//Gamma Term
double gammaCompReal = 1 + n/2.;
acb_t gammaComp; 
acb_init(gammaComp); 
acb_set_d_d(gammaComp, gammaCompReal, 0);
acb_gamma(gammaComp, gammaComp, 500); 
arb_t factorial; 
arb_init(factorial); 
arb_fac_ui(factorial, n, 500); 
acb_t cFactorial; 
acb_init(cFactorial); 
acb_set_arb(cFactorial, factorial); 
acb_div(gammaComp, gammaComp, cFactorial, 500); 

//Power Term
double powerTermPrefactorReal = 2*w; 
double powerTermPower = n/2.; 
arb_t pTermPower; 
arb_init(pTermPower); 
arb_set_d(pTermPower, powerTermPower); 
acb_t powerTerm;
acb_init(powerTerm); 
acb_set_d_d(powerTerm, 0, -1*powerTermPrefactorReal); 
acb_pow_arb(powerTerm, powerTerm, pTermPower, 500); 

//Hyper Term
double hyperArgAReal = 1 + n/2.; 
double hyperArgBReal = 1; 
double hyperArgZImag = -(1./2.) * w * y * y;

acb_t hyperArgA; 
acb_t hyperArgB; 
acb_t hyperArgZ; 
acb_init(hyperArgA);
acb_init(hyperArgB);
acb_init(hyperArgZ); 
acb_set_d_d(hyperArgA, hyperArgAReal, 0); 
acb_set_d_d(hyperArgB, hyperArgBReal, 0); 
acb_set_d_d(hyperArgZ, 0, hyperArgZImag); 
acb_t hyperTerm;
acb_init(hyperTerm); 
acb_hypgeom_1f1(hyperTerm, hyperArgA, hyperArgB, hyperArgZ, 0, 500); 

acb_t sumTemp; 
acb_init(sumTemp); 
acb_mul(sumTemp, gammaComp, powerTerm, 500); 
acb_mul(sumTemp, sumTemp, hyperTerm, 500); 

acb_add(sumTerm, sumTerm, sumTemp, 500);

acb_clear(gammaComp);
arb_clear(factorial);
acb_clear(cFactorial);
arb_clear(pTermPower);
acb_clear(powerTerm);
acb_clear(hyperArgA);
acb_clear(hyperArgB);
acb_clear(hyperArgZ);
acb_clear(hyperTerm);
acb_clear(sumTemp); 
}

acb_mul(ampFac, prefactorExp, sumTerm, 500); 

acb_clear(prefactorExp);
acb_clear(sumTerm); 
}

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> ampFacMatrices(std::vector<double> w, std::vector<double> y){

long ySize = y.size();
long wSize = w.size();

std::vector<std::vector<double>> ampFacReal(ySize, std::vector<double> (wSize));
std::vector<std::vector<double>> ampFacImag(ySize, std::vector<double> (wSize));

#pragma omp parallel for collapse(2) schedule(dynamic) 
for (int i=0; i<ySize; i++){
for (int j=0; j<wSize; j++){
acb_t ampFac;
acb_init(ampFac);
ampFacCalc(ampFac, w[j], y[i], 1000);
ampFacReal[i][j] = arf_get_d(arb_midref(acb_realref(ampFac)), ARF_RND_NEAR);
ampFacImag[i][j] = arf_get_d(arb_midref(acb_imagref(ampFac)), ARF_RND_NEAR);
}
}

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> ampFacMats(ampFacReal, ampFacImag);

return ampFacMats;
}

int main(int argc, char* argv[]){

std::string wFileName = argv[1];
std::string yFileName = argv[2];

std::ifstream wFile (wFileName);
std::istream_iterator<double> wStart(wFile), wEnd;
std::vector<double> wVals(wStart, wEnd);

std::ifstream yFile (yFileName);
std::istream_iterator<double> yStart(yFile), yEnd;
std::vector<double> yVals(yStart, yEnd);

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> fVals = ampFacMatrices(wVals,yVals);

std::ofstream fRealFile ("fReal.dat", std::ofstream::out);
std::ofstream fImagFile ("fImag.dat", std::ofstream::out);

int rowLen = fVals.first.size();
int colLen = fVals.first[0].size();

for(int i=0; i<rowLen; i++){
for(int j=0; j<colLen; j++){
fRealFile << fVals.first[i][j] << "\t";
fImagFile << fVals.second[i][j] << "\t";
}
fRealFile << "\n";
fImagFile << "\n";
}

return 0; 
}
