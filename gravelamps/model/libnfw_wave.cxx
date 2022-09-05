// Function definitions for calculating the amplification factor the isolated
// Navarro, Frenk, and White (NFW) lens model in the wave optics regime
//
// Mick Wright 2021

#include "src/nfw.h"

// Function computes the value of the lensing potential for specified values
// of the image position and scaling constant. This value is zero in the case
// where the image position is equal to one.
//
// Input:
//      double image_position : dimensionless position of the lensed image
//      double scaling_constant : characteristic scale of the profile
//
// Output:
//      std::complex<double> lensing_potential : complex value of the lensing
//                                               potential for the given
//                                               image position and scale
std::complex<double> LensingPotential(double image_position,
                                      double scaling_constant) {
    // Return zero if the image position is 1
    if (image_position == 1) {
        return 0;
    }

    // Compute the scaling constant term
    double scaling_constant_term = scaling_constant/2;

    // Compute the log term
    std::complex<double> image_complex = image_position;
    std::complex<double> log_term = log(image_complex/2.);
    log_term = log_term * log_term;

    // Compute the trigonemtric term based upon whether the image position is
    // greater or less than one
    std::complex<double> trig_term;
    if (image_position < 1) {
        trig_term = 1. - image_complex * image_complex;
        trig_term = -1. * atanh(sqrt(trig_term)) * atanh(sqrt(trig_term));
    } else {
        trig_term = image_complex * image_complex - 1.;
        trig_term = atan(sqrt(trig_term)) * atan(sqrt(trig_term));
    }

    // Construct the final result
    std::complex<double> lensing_potential = scaling_constant_term
                                             * (log_term + trig_term);

    return lensing_potential;
}

// Function computes the time delay for given image and source positions
// at the given scale modified by a specified phase.
//
// Input:
//      double image_position : dimensionless position of the lensed image
//      double source_position : dimensionless displacement from optical axis
//      double scaling_constant : characteristic scale of the profile
//      double phase : phase to modify such that the minimum time delay is zero
//
// Output:
//      std::complex<double> time_delay : calculated time delay between images
std::complex<double> TimeDelay(double image_position,
                               double source_position,
                               double scaling_constant,
                               double phase) {
    // Compute the positional term
    double positional_term = abs(image_position - source_position);
    positional_term = (positional_term * positional_term)/2.;

    // Compute the lensing potential
    std::complex<double> lensing_potential =
        LensingPotential(abs(image_position), scaling_constant);

    // Construct the final time delay
    std::complex<double> time_delay = 
        positional_term - lensing_potential + phase;

    return time_delay;
}

// Function computes the phase required such that the minimum time delay
// induced by the lensing is zero for the given scaling of the profile at the
// specified source position.
//
// Input:
//      double source_position : dimensionless displacement from optical axis
//      double scaling_constant : characteristic scale of the profile
//
// Output:
//      std::complex<double> phase : phase as described above
std::complex<double> Phase(double source_position, double scaling_constant) {
    using boost::math::tools::brent_find_minima;

    // Setting maximal precision
    int bits = std::numeric_limits<double>::digits;

    // Creating a lambda version of the TimeDelay function to narrow it to
    // a single dimensional function in image position
    auto LambdaTimeDelay =
        [source_position, scaling_constant](double image_position){
            return std::real(TimeDelay(image_position,
                                       source_position,
                                       scaling_constant,
                                       0));
        };

    // Feed this into the minima finding algorithm
    std::pair<double, double> minimal_pair =
        brent_find_minima(LambdaTimeDelay, 0.01, 10.0, bits);

    // Use that first value to compute the complex value to return
    std::complex<double> phase = TimeDelay(minimal_pair.first,
                                           source_position,
                                           scaling_constant,
                                           0);

    return -1.*phase;
}

// Function computes the value of the positive exponential term for the
// intermediate K function used in the calculation of the amplification factor
// for the specified dimensionless frequency, source position, and phase  with
// the exponentiation being done at the specified precision. The resultant value
// is assigned to the given acb_t object.
//
// Input:
//      acb_t positive_exponential_term : object for containing the complex
//                                        value of the positiuve exponential
//                                        term of the function
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from the optical
//                               axis
//      double phase : phase such that the minimum time delay induced by the
//                     lensing is zero
//      slong precision : arithmetic precision to use in arbitrary precision
//                        calculations
void KFunctionPositiveExponential(acb_t positive_exponential_term,
                                  double dimensionless_frequency,
                                  double source_position,
                                  double phase,
                                  slong precision) {
    double source_position_sq = source_position * source_position;
    double positive_exponential_inner_bracket =
        source_position_sq/2. + phase;
    double positive_exponential_term_imag =
        dimensionless_frequency * positive_exponential_inner_bracket;

    acb_set_d_d(positive_exponential_term, 0,
                positive_exponential_term_imag);
    acb_exp(positive_exponential_term, positive_exponential_term, precision);
}

// Function computes the value of the negative exponential term for the
// intermediate K function used in the calculation of the amplification factor
// for the specified dimensionless frequency, image position, and scaling
// constant with the exponentiation being done at the specified precision.
// The resultant value is assigned to the given acb_t object.
//
// Input:
//      acb_t negative_exponential_term : object for containing the complex
//                                        value of the negative exponential
//                                        term of the function
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double image_position : dimensionless position of image considered
//      double scaling_constant : characteristic scale of the profile
//      slong precision : arithmetic precision to use in arbitrary precision
//                        calculation
void KFunctionNegativeExponential(acb_t negative_exponential_term,
                                  double dimensionless_frequency,
                                  double image_position,
                                  double scaling_constant,
                                  slong precision) {
    using std::literals::complex_literals::operator""i;

    double lensing_potential_arg = sqrt(2 * image_position);
    std::complex<double> lensing_potential =
        LensingPotential(lensing_potential_arg, scaling_constant);

    std::complex<double> exponential_arg = -1i * dimensionless_frequency
                                           * lensing_potential;

    acb_set_d_d(negative_exponential_term,
                std::real(exponential_arg),
                std::imag(exponential_arg));
    acb_exp(negative_exponential_term, negative_exponential_term, precision);
}

// Function computes the value of the bessel function term for the intermediate
// K function used in the calculation of the amplification factor for the
// specified dimensionless frequency, source position, and image position with
// the function calculated to the specified precision. The resultant value is
// assigned to the given acb_t object.
//
// Input:
//      acb_t bessel_term : object for containing the complex value of the
//                          bessel function term of the function
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//      double image_position : dimensionless posiition of image considered
//      slong precision : arithmetic precision to use in arbitrary precision
//                        calculation
void KFunctionBessel(acb_t bessel_term,
                     double dimensionless_frequency,
                     double source_position,
                     double image_position,
                     slong precision) {
    double bessel_term_real = dimensionless_frequency * source_position
                              * sqrt(2 * image_position);

    // Zero needed for bessel function level
    acb_t zero;
    acb_init(zero);
    acb_zero(zero);

    acb_set_d(bessel_term, bessel_term_real);
    acb_hypgeom_bessel_j(bessel_term, zero, bessel_term, precision);

    // Clearing declared acb to release memory
    acb_clear(zero);
}

// Function computes the value of the intermediate K function used in the
// calculation of the amplification factor for the wave optics case for
// specified dimensionless frequency, source position, image position
// scaling constant, and phase with the functions calculated with specified
// precision. The value of the function is assigned to the
// given acb_t object.
//
// Input:
//      acb_t intermediate_function_value : object for containing the complex
//                                          value of the K function
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//      double image_position : dimensionless position of image considered
//      double scaling_constant : characteristic scale of the profile
//      double phase : phase such that the minimum time delay induced by the
//                     lensing is zero
//      slong precision : arithmetic precision to use in arbitrary precision
//                        calculation
void KFunctionCalculation(acb_t k_function_value,
                          double dimensionless_frequency,
                          double source_position,
                          double image_position,
                          double scaling_constant,
                          double phase,
                          slong precision) {
    // Construct the dimensionless frequency term
    double dimensionless_frequency_term_imag = -1 * dimensionless_frequency;

    acb_t dimensionless_frequency_term;
    acb_init(dimensionless_frequency_term);
    acb_set_d_d(dimensionless_frequency_term, 0,
                dimensionless_frequency_term_imag);

    // Construct the positive exponential term
    acb_t positive_exponential_term;
    acb_init(positive_exponential_term);
    KFunctionPositiveExponential(positive_exponential_term,
                                 dimensionless_frequency,
                                 source_position,
                                 phase,
                                 precision);

    // Construct the bessel function term
    acb_t bessel_term;
    acb_init(bessel_term);
    KFunctionBessel(bessel_term,
                    dimensionless_frequency,
                    source_position,
                    image_position,
                    precision);

    // Construct the negative exponential term
    acb_t negative_exponential_term;
    acb_init(negative_exponential_term);
    KFunctionNegativeExponential(negative_exponential_term,
                                 dimensionless_frequency,
                                 image_position,
                                 scaling_constant,
                                 precision);

    // Construct the final value of the intermediate function
    acb_mul(k_function_value,
            dimensionless_frequency_term,
            positive_exponential_term,
            precision);
    acb_mul(k_function_value,
            k_function_value,
            bessel_term,
            precision);
    acb_mul(k_function_value,
            k_function_value,
            negative_exponential_term,
            precision);

    // Clear declared acbs to free memory
    acb_clear(dimensionless_frequency_term);
    acb_clear(positive_exponential_term);
    acb_clear(bessel_term);
    acb_clear(negative_exponential_term);
}

// Function computes the value of the integrand for the integral in the
// calculation of the amplification factor for the specified parameters.
// These parameters are passed as a single object to the function and then
// extracted for consistency with arb's integration methodology. The order
// parameter required by said methodology is not implemented, however
// the functions are calculated with precision as specified. The resultant
// value is assigned to the pointer provided to the function.
//
// Input:
//      acb_ptr integrand : pointer to object containing the complex integrand
//                          value
//      const acb_t integration_parameter : parameter over which the
//                                          integration is occurring. This is
//                                          over image positions.
//      void * parameter_set : object containing the set of parameters to be
//                             used to determine the integrand function value.
//                             These parameters are the dimensionless frequency
//                             , source position, scaling constant, and phase
//                             value
//      slong order : Unimplemented parameter to specify whether to calculate
//                    the integral as the first n coefficients of a Taylor
//                    expansion
//      slong precision : arithemtic precision to use in arbitrary precision
//                        calculations
int AmplificationWaveIntegrand(acb_ptr integrand,
                               const acb_t integration_parameter,
                               void * parameter_set,
                               slong order,
                               slong precision) {
    // Extracting the parameters from the parameter_set object
    std::vector<double> parameter_vector =
        ((std::vector<double> *) parameter_set)[0];

    double dimensionless_frequency = parameter_vector[0];
    double source_position = parameter_vector[1];
    double scaling_constant = parameter_vector[2];
    double phase = parameter_vector[3];

    // Extracting the value of the integration parameter for use in functions
    double image_position =
        arf_get_d(arb_midref(acb_realref(integration_parameter)), ARF_RND_NEAR);

    // Caclulate the K function term
    acb_t k_function_term;
    acb_init(k_function_term);
    KFunctionCalculation(k_function_term,
                         dimensionless_frequency,
                         source_position,
                         image_position,
                         scaling_constant,
                         phase,
                         precision);

    // Calculate the exponential term
    acb_t exponential_term;
    acb_init(exponential_term);
    acb_set_d_d(exponential_term, 0, dimensionless_frequency);
    acb_mul(exponential_term, exponential_term, integration_parameter, precision);
    acb_exp(exponential_term, exponential_term, precision);

    // Calculate the final value
    acb_mul(integrand, k_function_term, exponential_term, precision);

    // Clear declared acb objects to clear memory
    acb_clear(k_function_term);
    acb_clear(exponential_term);

    return 0;
}

// Function computes the value of the exponential term of the first correction
// to the amplification factor for the specified value of dimensionlesss
// frequency and maximum image position at the specified precision. The
// resultant value is assigned to the given acb_t object.
//
// Input:
//      acb_t exponential_term : object for containing the complex value of the
//                               exponential term
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double max_image_position : the maximal value of image position used as
//                                  the upper limit of the integration being
//                                  performed in the calculation of the
//                                  amplification factor
//      slong precision : arithmetic precision to be used in arbitrary precision
//                        calculations
void FirstCorrectionExponential(acb_t exponential_term,
                                double dimensionless_frequency,
                                double max_image_position,
                                slong precision) {
    double exponential_term_imag = dimensionless_frequency * max_image_position;

    acb_set_d_d(exponential_term, 0, exponential_term_imag);
    acb_exp(exponential_term, exponential_term, precision);
}

// Function computes the value of the first correction term to the
// amplification factor for the specified values of the dimensionless
// frequency, source position, scaling constant and phase at the specified
// precision. This correction term relies on the maximal image position value
// being considered in the integration. The resultant value is assigned to the
// given acb_t object.
//
// Input:
//      acb_t first_correction_term : object for containing the complex value
//                                    of the first correction term to the
//                                    amplification factor
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//      double max_image_position : the maximal value of image position used as
//                                  the upper limit of the integration being
//                                  performed in the calculation of the
//                                  amplification factor
//      double scaling_constant : charcateristic scale of the profile
//      double phase : phase required for the minimum time delay induced by the
//                     lensing to be zero
//      slong precision : arithmetic precision to be used in arbitrary precision
//                        calculations
void AmplificationWaveFirstCorrection(acb_t first_correction_term,
                                      double dimensionless_frequency,
                                      double source_position,
                                      double max_image_position,
                                      double scaling_constant,
                                      double phase,
                                      slong precision) {
    // Compute the K-function value
    acb_t k_function_term;
    acb_init(k_function_term);
    KFunctionCalculation(k_function_term,
                         dimensionless_frequency,
                         source_position,
                         max_image_position,
                         scaling_constant,
                         phase,
                         precision);

    // Compute the exponential term
    acb_t exponential_term;
    acb_init(exponential_term);
    FirstCorrectionExponential(exponential_term,
                               dimensionless_frequency,
                               max_image_position,
                               precision);

    // Create the denominator term
    acb_t dimensionless_frequency_term;
    acb_init(dimensionless_frequency_term);
    acb_set_d_d(dimensionless_frequency_term, 0, dimensionless_frequency);

    // Construct the correction term
    acb_mul(first_correction_term,
            k_function_term,
            exponential_term,
            precision);
    acb_div(first_correction_term,
            first_correction_term,
            dimensionless_frequency_term,
            precision);
    acb_neg(first_correction_term, first_correction_term);

    // Clear the declared acb_t objects to free memory
    acb_clear(k_function_term);
    acb_clear(exponential_term);
    acb_clear(dimensionless_frequency_term);
}

// Function computes the value of the derivative term of the second correction
// term to the amplification factor at the specified dimensionles frequency,
// source position, scaling constant, image position, and phase values. The
// derivative is calculated using a central finite differences method. The
// resultant value is assigned to the given acb_t object.
//
// Input:
//      acb_t derivative_term : object for containing the complex value of the
//                              derivative term of the second correction term
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//      double max_image_position : the maximal image position value being used
//                                  as the upper limit of integration in the
//                                  calculation of the amplification factor
//      double scaling_constant : characteristic scale of the profile
//      double phase : phase required for the minimum time delay induced by the
//                     lensing to be zero
//      slong precision : arithmetic precision to use in arbitrary precision
//                        calculations
void SecondCorrectionDerivative(acb_t derivative_term,
                                double dimensionless_frequency,
                                double source_position,
                                double max_image_position,
                                double scaling_constant,
                                double phase,
                                slong precision) {
    // Setting the stepsize
    const int radix = std::numeric_limits<double>::radix;
    const int mantissa = std::numeric_limits<double>::digits;
    double stepsize = std::pow(1.0/radix, mantissa/2);

    // Calculating f(x+h) & f(x-h)
    double x_plus = max_image_position + stepsize;
    double x_minus = max_image_position - stepsize;

    acb_t x_plus_term;
    acb_init(x_plus_term);
    KFunctionCalculation(x_plus_term,
                         dimensionless_frequency,
                         source_position,
                         x_plus,
                         scaling_constant,
                         phase,
                         precision);

    acb_t x_minus_term;
    acb_init(x_minus_term);
    KFunctionCalculation(x_minus_term,
                         dimensionless_frequency,
                         source_position,
                         x_minus,
                         scaling_constant,
                         phase,
                         precision);

    // Setting the denominator of the finite differences
    double denominator_value = 2 * stepsize;

    arb_t denominator;
    arb_init(denominator);
    arb_set_d(denominator, denominator_value);

    // Calculate the final value
    acb_sub(derivative_term, x_plus_term, x_minus_term, precision);
    acb_div_arb(derivative_term, derivative_term, denominator, precision);

    // Clear declared acb_t and arb_t objects to free memory
    acb_clear(x_plus_term);
    acb_clear(x_minus_term);
    arb_clear(denominator);
}

// Function computes the value of the exponential term of the second correction
// term to the amplification factor at the specified values of dimensionless
// frequency and image position with the desired precision in the
// exponentiation. The resultant value is assigned to the given acb_t object.
//
// Input:
//      acb_t exponential_term : object for containing the complex value of the
//                               second correction term
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double max_image_position : the maximal image position value being used
//                                  as the upper limit of integration in the
//                                  calculation of the amplification factor
//      slong precision : arithmetic precision to use in arbitrary precision
//                        calculations
void SecondCorrectionExponential(acb_t exponential_term,
                                 double dimensionless_frequency,
                                 double max_image_position,
                                 slong precision) {
    double exponential_term_imag = dimensionless_frequency * max_image_position;

    acb_set_d_d(exponential_term, 0, exponential_term_imag);
    acb_exp(exponential_term, exponential_term, precision);
}

// Function computes the value of the second correction term to the
// amplification factor at the specified dimensionless frequency, source
// position, scaling constant, and phase values. The correction term also
// relies on the maximal image position being considered as part of the
// integration in the main amplification factor calculation. The resultant
// value is assigned to the given acb_t object.
//
// Input:
//      acb_t second_correction_term : object for containing the complex value
//                                     of the second correction term
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//      double max_image_position : the maximal image position value being used
//                                  as the upper limit of integration in the
//                                  calculation of the amplification factor
//      double scaling_constant : characteristic scale of the profile
//      double phase : phase required for the minimum time delay induced by the
//                     lensing to be zero
//      slong precision : arithemetic precision to use in arbitrary precision
//                        calculations
void AmplificationWaveSecondCorrection(acb_t second_correction_term,
                                       double dimensionless_frequency,
                                       double source_position,
                                       double max_image_position,
                                       double scaling_constant,
                                       double phase,
                                       slong precision) {
    // Calculate the derivative term
    acb_t derivative_term;
    acb_init(derivative_term);
    SecondCorrectionDerivative(derivative_term,
                               dimensionless_frequency,
                               source_position,
                               max_image_position,
                               scaling_constant,
                               phase,
                               precision);

    // Calculate exponential term
    acb_t exponential_term;
    acb_init(exponential_term);
    SecondCorrectionExponential(exponential_term,
                                dimensionless_frequency,
                                max_image_position,
                                precision);

    // Calculate the dimensionless frequency term
    double dimensionless_frequency_isq = -1 * dimensionless_frequency
                                         * dimensionless_frequency;

    acb_t dimensionless_frequency_term;
    acb_init(dimensionless_frequency_term);
    acb_set_d(dimensionless_frequency_term, dimensionless_frequency);

    // Construct the final correction term
    acb_mul(second_correction_term,
            derivative_term,
            exponential_term,
            precision);
    acb_div(second_correction_term,
            second_correction_term,
            dimensionless_frequency_term,
            precision);

    // Clear declared acb_t objects to free memory
    acb_clear(derivative_term);
    acb_clear(exponential_term);
    acb_clear(dimensionless_frequency_term);
}

// Function computes the amplification factor for an isolated axially symmetric
// Navarro, Frenk, White (NFW) lens for given values of dimensionless frequency,
// source position, and scaling constant with functions calculated with the
// specified precision. The infinite integral is approximated as a finite
// integral up to the specified value. The resultant value is assigned to the
// given amplification factor object.
//
// Input:
//      acb_t amplification_factor : object for containing the complex value of
//                                   the amplification factor
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//      double scaling_constant : characteristic scale of the profile
//      double max_image_position : the maximal image position value used as
//                                  the upper limit of the integration being
//                                  performed to calculate the amplification
//                                  factor
//      slong precision : arithmetic precision to use in arbitrary precision
//                        calculations
void AmplificationFactorWave(acb_t amplification_factor,
                             double dimensionless_frequency,
                             double source_position,
                             double scaling_constant,
                             double max_image_position,
                             slong precision) {
    // Calculate the phase for the source position and scaling constant
    double phase = std::real(Phase(source_position, scaling_constant));

    // Create a vector containing the lensing parameters to be passed to the
    // integral
    std::vector<double> parameter_set {dimensionless_frequency,
                                       source_position,
                                       scaling_constant,
                                       phase};

    // Setting options for integration
    slong goal = precision;
    mag_t tol;
    mag_set_ui_2exp_si(tol, 1, -1.*precision);

    acb_calc_integrate_opt_t options;
    acb_calc_integrate_opt_init(options);

    options -> use_heap = 1;
    options -> depth_limit = 128*precision;
    options -> eval_limit = precision*precision*precision;

    // Setting lower limit of integration - due to an error at precisely zero
    // this is set just above
    double lower_limit_val = 0.000001;

    acb_t lower_limit;
    acb_init(lower_limit);
    acb_set_d(lower_limit, lower_limit_val);

    acb_t upper_limit;
    acb_init(upper_limit);
    acb_set_d(upper_limit, max_image_position);

    // Calculate the integration term
    acb_t integration_term;
    acb_init(integration_term);
    acb_calc_integrate(integration_term,
                       AmplificationWaveIntegrand,
                       &parameter_set,
                       lower_limit,
                       upper_limit,
                       goal,
                       tol,
                       options,
                       goal);

    // Caclulate the first correction term
    acb_t first_correction_term;
    acb_init(first_correction_term);
    AmplificationWaveFirstCorrection(first_correction_term,
                                     dimensionless_frequency,
                                     source_position,
                                     max_image_position,
                                     scaling_constant,
                                     phase,
                                     precision);

    // Calculate the second correction term
    acb_t second_correction_term;
    acb_init(second_correction_term);
    AmplificationWaveSecondCorrection(second_correction_term,
                                      dimensionless_frequency,
                                      source_position,
                                      max_image_position,
                                      scaling_constant,
                                      phase,
                                      precision);

    // Construct the final value of the amplification factor
    acb_add(amplification_factor,
            integration_term,
            first_correction_term,
            precision);
    acb_add(amplification_factor,
            amplification_factor,
            second_correction_term,
            precision);

    // Clear declared acb_t objects to clear memory
    acb_clear(lower_limit);
    acb_clear(upper_limit);
    acb_clear(integration_term);
    acb_clear(first_correction_term);
    acb_clear(second_correction_term);
}

