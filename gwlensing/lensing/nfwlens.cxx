#include "nfw.h"

void acb_psi(acb_t psi, const acb_t x, acb_t ks, slong prec){ 
//Function to calculate psi(x, ks) - the lensing potential 

//Get the Value of X so as to determine what value to use 
double xVal = arf_get_d(arb_midref(acb_realref(x)), ARF_RND_NEAR);

//Initialise Component Parts of Psi
acb_t psiPrefactor, psiFirstTerm, psiSecondTerm;
acb_init(psiPrefactor);
acb_init(psiFirstTerm);
acb_init(psiSecondTerm); 

//Initialise Modifiers for Equation 
acb_t one, two; 
acb_init(two);
acb_init(one); 
acb_set_d(two, 2); 
acb_one(one); 

//Check for if x=1, if so return 0; 
if (xVal == 1){
acb_zero(psi);

acb_clear(psiPrefactor);
acb_clear(psiFirstTerm);
acb_clear(psiSecondTerm);
acb_clear(one);
acb_clear(two); 
return; 
}

//Calculate the prefactor (ks/2) 
acb_div(psiPrefactor, ks, two, prec); 

//Calculate the first term, log(x/2)^2 
acb_div(psiFirstTerm, x, two, prec); 
acb_log(psiFirstTerm, psiFirstTerm, prec); 
acb_sqr(psiFirstTerm, psiFirstTerm, prec); 

//Check whether x is greater or less than one; 
if (xVal<1){
//If less than one, psiSecondTerm = -ArcTanh(Sqrt(1-x^2))^2 

acb_one(psiSecondTerm); 
acb_submul(psiSecondTerm, x, x, prec);
acb_sqrt(psiSecondTerm, psiSecondTerm, prec); 
acb_atanh(psiSecondTerm, psiSecondTerm, prec); 
acb_sqr(psiSecondTerm, psiSecondTerm, prec);
acb_neg(psiSecondTerm, psiSecondTerm); 
}
else{
//If greater than one, psiSecondTerm = ArcTan(Sqrt(x^2-1))^2 

acb_sqr(psiSecondTerm, x, prec);
acb_sub(psiSecondTerm, psiSecondTerm, one, prec); 
acb_sqrt(psiSecondTerm, psiSecondTerm, prec);
acb_atan(psiSecondTerm, psiSecondTerm, prec); 
acb_sqr(psiSecondTerm, psiSecondTerm, prec); 
}

//Construct psi = psiPrefactor*(psiFirstTerm+psiSecondTerm) 
acb_add(psi, psiFirstTerm, psiSecondTerm, prec);
acb_mul(psi, psiPrefactor, psi, prec);

acb_clear(psiPrefactor);
acb_clear(psiFirstTerm);
acb_clear(psiSecondTerm);
acb_clear(one);
acb_clear(two); 
}

void acb_k_calc(acb_t k, acb_t w, acb_t y, const acb_t x, acb_t ks, slong prec){
//Calculates the intermediate function k(w,y,x,ks) for the amplification factor cacalculation. Function k is given by -iw*(exp[iw(y^2/2)]*J0(wy*sqrt(2x))*exp(-iw*psi(sqrt(2x), ks))) 

//Initialise Components to Construct K
acb_t kPrefactor, kFirstExpTerm, kBesselTerm, kSecondExpTerm; 
acb_init(kPrefactor);
acb_init(kFirstExpTerm);
acb_init(kBesselTerm);
acb_init(kSecondExpTerm); 

//Initialise an acb for two and zero for use
acb_t two, zero; 
acb_init(two);
acb_init(zero); 
acb_set_d(two, 2);
acb_zero(zero); 

//Prefactor is given by -iw
acb_mul_onei(kPrefactor, w);
acb_neg(kPrefactor, kPrefactor); 

//First Exp Term is exp(iw*(y^2/2)) 
acb_sqr(kFirstExpTerm, y, prec); 
acb_div(kFirstExpTerm, kFirstExpTerm, two, prec); 
acb_mul(kFirstExpTerm, w, kFirstExpTerm, prec);
acb_mul_onei(kFirstExpTerm, kFirstExpTerm); 
acb_exp(kFirstExpTerm, kFirstExpTerm, prec); 

//Bessel Term is J0(wy*sqrt(2x))
acb_mul(kBesselTerm, two, x, prec); 
acb_sqrt(kBesselTerm, kBesselTerm, prec); 
acb_mul(kBesselTerm, kBesselTerm, y, prec);
acb_mul(kBesselTerm, kBesselTerm, w, prec); 
acb_hypgeom_bessel_j(kBesselTerm, zero, kBesselTerm, prec); 

//Second Exp Term is exp(-iw*psi(sqrt(2x),ks)) 
acb_mul(kSecondExpTerm, x, two, prec); 
acb_sqrt(kSecondExpTerm, kSecondExpTerm, prec); 
acb_psi(kSecondExpTerm, kSecondExpTerm, ks, prec); 
acb_mul(kSecondExpTerm, kSecondExpTerm, w, prec); 
acb_mul_onei(kSecondExpTerm, kSecondExpTerm); 
acb_neg(kSecondExpTerm, kSecondExpTerm); 
acb_exp(kSecondExpTerm, kSecondExpTerm, prec); 

//Construct by Multiplying Each Term
acb_mul(k, kPrefactor, kFirstExpTerm, prec);
acb_mul(k, k, kBesselTerm, prec);
acb_mul(k, k, kSecondExpTerm, prec);

//Clear Up
acb_clear(kPrefactor);
acb_clear(kBesselTerm); 
acb_clear(kFirstExpTerm);
acb_clear(kSecondExpTerm);
acb_clear(two);
acb_clear(zero); 
}

int f_nfwIntegrand(acb_ptr res, const acb_t x, void * param, slong order, slong prec){
//Intermediate Function to be Integrated in order to calculate the NFW amplification factor given by k(w,y,x,ks)*exp(iwx) 

//Extract w, y, and ks from param
std::vector<double> wyks = ((std::vector<double> *) param)[0]; 
double wVal = wyks[0];
double yVal = wyks[1]; 
double ksVal = wyks[2]; 

//Set up acb_t for w y and ks
acb_t w, y, ks; 
acb_init(w);
acb_init(y);
acb_init(ks); 
acb_set_d(w, wVal);
acb_set_d(y, yVal);
acb_set_d(ks, ksVal); 

//Set up First Term - k(w,y,x,ks); 
acb_t intFirstTerm;
acb_init(intFirstTerm);
acb_k_calc(intFirstTerm, w, y, x, ks, prec); 

//Set up Second Term - exp(iwx); 
acb_t intSecondTerm;
acb_init(intSecondTerm); 
acb_mul(intSecondTerm, w, x, prec); 
acb_mul_onei(intSecondTerm, intSecondTerm); 
acb_exp(intSecondTerm, intSecondTerm, prec); 

//Send result to res
acb_mul(res, intFirstTerm, intSecondTerm, prec);

//Clean up and return
acb_clear(w);
acb_clear(y);
acb_clear(ks);
acb_clear(intFirstTerm);
acb_clear(intSecondTerm); 

return 0; 
}

void acb_f2t1(acb_t termOne, acb_t w, acb_t y, acb_t upperLimitX, acb_t ks, slong prec){
//Calculates First Correction Term in the Amplification Factor Calculation given by -((k(w,y,xUpperLimit,ks)*exp(iw*xUpperLimit))/iw) 

//Initialise Terms
acb_t kTerm, expTerm, denomTerm;
acb_init(kTerm);
acb_init(expTerm);

//Get the k Term
acb_k_calc(kTerm, w, y, upperLimitX, ks, prec);

//Get the exponential Term - exp(iw*UpperLimitX)
acb_mul(expTerm, w, upperLimitX, prec);
acb_mul_onei(expTerm, expTerm); 
acb_exp(expTerm, expTerm, prec); 

//Calculate the Denominator - iw 
acb_init(denomTerm); 
acb_mul_onei(denomTerm, w); 

//Construct Term One
acb_mul(termOne, kTerm, expTerm, prec); 
acb_div(termOne, termOne, denomTerm, prec);
acb_neg(termOne, termOne); 
}

void acb_f2t2(acb_t termTwo, acb_t w, acb_t y, acb_t upperLimitX, acb_t ks, slong prec){
//Calculates the Second Term in the Amplification Factor Calculation given by d(k(w,y,x,ks)*exp(iwz))/dx)/(iw)**2 

//Initialise f(x+h) f(x-h) df/dx
acb_t fXPlusH, fXMinusH, dfdx; 
acb_init(fXPlusH);
acb_init(fXMinusH); 
acb_init(dfdx); 

//Initialise h - step size value for differentiation 
acb_t h; 
acb_init(h);
acb_set_d(h,0.00001); 

//Initialise (x+h) (x-h) (2h) 
acb_t XPlusH, XMinusH, twoH; 
acb_init(XPlusH);
acb_init(XMinusH);
acb_init(twoH); 

//Get values for the above
acb_add(twoH, h, h, prec); 
acb_add(XPlusH, upperLimitX, h, prec); 
acb_sub(XMinusH, upperLimitX, h, prec);

//Get f(x+h), f(x-h); 
acb_k_calc(fXPlusH, w, y, XPlusH, ks, prec);
acb_k_calc(fXMinusH, w, y, XMinusH, ks, prec); 

//Calculate f(x+h)-f(x-h)/2h 
acb_sub(dfdx, fXPlusH, fXMinusH, prec); 
acb_div(dfdx, dfdx, twoH, prec); 

//Calculate Exp Term - exp(iwx) 
acb_t expTerm;
acb_init(expTerm);
acb_mul(expTerm, w, upperLimitX, prec); 
acb_mul_onei(expTerm, expTerm); 
acb_exp(expTerm, expTerm, prec); 

//Calculate Denominator - (iw)**2
acb_t denomTerm;
acb_init(denomTerm); 
acb_mul_onei(denomTerm, w);
acb_sqr(denomTerm, denomTerm, prec); 

//Construct Term Two 
acb_mul(termTwo, dfdx, expTerm, prec);
acb_div(termTwo, termTwo, denomTerm, prec); 
}

void ampFacCalc(acb_t ampFac, double wVal, double yVal, double ksVal, double upperLimitVal, slong prec){
//Function to calculate the NFW amplification factor 

//Construct Vector of Lensing Parameters to Pass to Integrand Function
std::vector<double> params {wVal, yVal, ksVal}; 

//Construct acbs for the Lensing Parameters
acb_t w, y, ks;
acb_init(w);
acb_init(y);
acb_init(ks); 

acb_set_d(w,wVal);
acb_set_d(y,yVal);
acb_set_d(ks,ksVal); 

//Set goal and tolerance for the integration
slong goal = prec; 

mag_t tol; 
mag_set_ui_2exp_si(tol, 1, -1*prec); 

//Options
acb_calc_integrate_opt_t options;
acb_calc_integrate_opt_init(options);
options -> use_heap = 1;
options -> depth_limit = 128*prec;
options -> eval_limit = prec*prec*prec;

//Create objects for the upper and lower limits of the integration 
acb_t lowerLimit, upperLimit; 
acb_init(lowerLimit);
acb_init(upperLimit);
acb_set_d(lowerLimit, 0.000001); 
acb_set_d(upperLimit, upperLimitVal); 

//Calculate integration term
acb_t intTerm;
acb_init(intTerm); 
acb_calc_integrate(intTerm, f_nfwIntegrand, &params, lowerLimit, upperLimit, goal, tol, options, prec); 

//Calculate the first term 
acb_t termOne; 
acb_init(termOne);
acb_f2t1(termOne, w, y, upperLimit, ks, prec); 

//Calculate the second term
acb_t termTwo;
acb_init(termTwo); 
acb_f2t2(termTwo, w, y, upperLimit, ks, prec); 

//Construct the Amplification Factor as the sum of the three terms 
acb_add(ampFac, intTerm, termOne, prec);
acb_add(ampFac, ampFac, termTwo, prec); 
}

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> ampFacMatrices(std::vector<double> w, std::vector<double> y, double ksVal, double upperLimitVal, slong prec){

long ySize = y.size();
long wSize = w.size();

std::vector<std::vector<double>> ampFacReal(ySize, std::vector<double> (wSize));
std::vector<std::vector<double>> ampFacImag(ySize, std::vector<double> (wSize));

#pragma omp parallel for collapse(2) schedule(dynamic)
for (int i=0; i<ySize; i++){
for (int j=0; j<wSize; j++){
acb_t ampFac;
acb_init(ampFac);
ampFacCalc(ampFac, w[j], y[i], ksVal, upperLimitVal, prec);
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

double ksVal = atof(argv[3]); 
double upperLimitVal = atof(argv[4]); 
slong prec = atoi(argv[5]); 

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> fVals = ampFacMatrices(wVals,yVals,ksVal,upperLimitVal,prec);

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
