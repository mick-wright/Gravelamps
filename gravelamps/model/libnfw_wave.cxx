// Function definitions for calculating the amplification factor the isolated
// Navarro, Frenk, and White (NFW) lens model in the wave optics regime
//
// Mick Wright 2021

#include "src/nfw.h"

// Function computes the value of psi - the value of the lensing potential
void LensingPotential(acb_t lensing_potential,
                      const acb_t scaled_surface_density,
                      acb_t scaling_constant,
                      slong precision) {
    // For comparison purposes, we need a double version of the scaled
    // surface density so as to decide on the which case we are in
    double scaled_surf_den_val = arf_get_d(
        arb_midref(acb_realref(scaled_surface_density)), ARF_RND_NEAR);

    // In the case of the scaled surface density being 1, the value of psi
    // is zero, so we can quickly return this
    if (scaled_surf_den_val == 1) {
        acb_zero(lensing_potential);
        return;
    }

    // Otherwise begin constructing the lensing potential function. This is
    // done in three seperate term. A prefactor which is simply half of the
    // scaling constant. The first main term, which is given by the square of
    // the log of half of the scaled surface density - and the third term which
    // is trigonometric in nature and dependant upon whether the scaled surface
    // density is greater or less than one.
    acb_t lensing_potential_prefactor;
    acb_t lensing_potential_log_term;
    acb_t lensing_potential_trig_term;
    acb_t two;

    acb_init(lensing_potential_prefactor);
    acb_init(lensing_potential_log_term);
    acb_init(lensing_potential_trig_term);
    acb_init(two);

    acb_set_d(two, 2);

    // The prefactor is given by the scaling_constant/2
    acb_div(lensing_potential_prefactor, scaling_constant, two, precision);

    // The logarithmic term is given by log(scaled_surface_density/2)^2
    acb_div(
        lensing_potential_log_term, scaled_surface_density, two, precision);
    acb_log(lensing_potential_log_term, lensing_potential_log_term, precision);
    acb_sqr(lensing_potential_log_term, lensing_potential_log_term, precision);

    // Perform the check as to whether the scaled surface density is greater or
    // less than one
    if (scaled_surf_den_val < 1) {
        // If less than one the value of the trigonometric term is given by
        // -ArcTanh(sqrt(1 - scaled_surface_density^2))^2
        acb_one(lensing_potential_trig_term);
        acb_submul(lensing_potential_trig_term,
                   scaled_surface_density,
                   scaled_surface_density,
                   precision);
        acb_sqrt(lensing_potential_trig_term,
                 lensing_potential_trig_term,
                 precision);
        acb_atanh(lensing_potential_trig_term,
                  lensing_potential_trig_term,
                  precision);
        acb_sqr(lensing_potential_trig_term,
                lensing_potential_trig_term,
                precision);
        acb_neg(lensing_potential_trig_term, lensing_potential_trig_term);
    } else {
        // If greater than one the value of the trigonometric term is given by
        // ArcTan(sqrt(scaled_surface_density^2 - 1))^2
        acb_t one;
        acb_init(one);
        acb_one(one);

        acb_sqr(
            lensing_potential_trig_term, scaled_surface_density, precision);
        acb_sub(lensing_potential_trig_term,
                lensing_potential_trig_term,
                one,
                precision);
        acb_sqrt(lensing_potential_trig_term,
                 lensing_potential_trig_term,
                 precision);
        acb_atan(lensing_potential_trig_term,
                 lensing_potential_trig_term,
                 precision);
        acb_sqr(lensing_potential_trig_term,
                lensing_potential_trig_term,
                precision);

        acb_clear(one);
    }

    // With each part of the lensing potential calculated, construct the final
    // value:
    // lensing_potential = prefactor*(log_term + trig_term)
    acb_add(lensing_potential,
            lensing_potential_log_term,
            lensing_potential_trig_term,
            precision);
    acb_mul(lensing_potential,
            lensing_potential_prefactor,
            lensing_potential,
            precision);

    // For memory management, destroy the constructed acb values
    acb_clear(lensing_potential_prefactor);
    acb_clear(lensing_potential_log_term);
    acb_clear(lensing_potential_trig_term);
    acb_clear(two);
}

// Function computes the value of the intermediate function k(w,y,z,ks) for the
// amplification factor calculation. The function k is given by
// -iw*(exp[iw(y^2/2 + phimin)]*J0(wy*sqrt(2x))*exp(-iw*psi(sqrt(2x),ks)))
void IntermediateFunctionCalculation(acb_t intermediate_function_value,
                                     acb_t dimensionless_frequency,
                                     acb_t source_position,
                                     const acb_t integration_parameter,
                                     acb_t scaling_constant,
                                     double minimum_phase,
                                     slong precision) {
    // Initialise the components of the function - as can be seen there are
    // four - a prefactor (-iw), a first expontential term (exp(iw*(y^2/2))), a
    // bessel term (J0(wy*sqrt(2x))), and a second expontential term
    // (exp(-iw * psi(sqrt(2x), ks)))
    acb_t prefactor;
    acb_t first_exponential_term;
    acb_t bessel_term;
    acb_t second_exponential_term;

    acb_init(prefactor);
    acb_init(first_exponential_term);
    acb_init(bessel_term);
    acb_init(second_exponential_term);

    // In addition, we need to initalise a zero and a two value for use in the
    // calculations
    acb_t two;
    acb_t zero;

    acb_init(zero);
    acb_init(two);

    acb_set_d(two, 2.);
    acb_zero(zero);

    // Finally we need to calculate phimin - the phase required for a minimum
    // time delay of zero;
    acb_t minimum_phase_acb;
    acb_init(minimum_phase_acb);
    acb_set_d(minimum_phase_acb, minimum_phase);

    // Construct the prefactor term
    acb_mul_onei(prefactor, dimensionless_frequency);
    acb_neg(prefactor, prefactor);

    // Construct the first exponential term
    acb_sqr(first_exponential_term, source_position, precision);
    acb_div(first_exponential_term, first_exponential_term, two, precision);
    acb_add(first_exponential_term,
            first_exponential_term,
            minimum_phase_acb,
            precision);
    acb_mul(first_exponential_term,
            dimensionless_frequency,
            first_exponential_term,
            precision);
    acb_mul_onei(first_exponential_term, first_exponential_term);
    acb_exp(first_exponential_term, first_exponential_term, precision);

    // Construct the bessel term
    acb_mul(bessel_term, two, integration_parameter, precision);
    acb_sqrt(bessel_term, bessel_term, precision);
    acb_mul(bessel_term, bessel_term, source_position, precision);
    acb_mul(bessel_term, bessel_term, dimensionless_frequency, precision);
    acb_hypgeom_bessel_j(bessel_term, zero, bessel_term, precision);

    // Construct the second exponential term
    acb_mul(second_exponential_term, integration_parameter, two, precision);
    acb_sqrt(second_exponential_term, second_exponential_term, precision);
    LensingPotential(second_exponential_term,
                     second_exponential_term,
                     scaling_constant,
                     precision);
    acb_mul(second_exponential_term,
            second_exponential_term,
            dimensionless_frequency,
            precision);
    acb_mul_onei(second_exponential_term, second_exponential_term);
    acb_neg(second_exponential_term, second_exponential_term);
    acb_exp(second_exponential_term, second_exponential_term, precision);

    // Construct the value of the intermediate function by multiplying the
    // terms together
    acb_mul(intermediate_function_value,
            prefactor,
            first_exponential_term,
            precision);
    acb_mul(intermediate_function_value,
            intermediate_function_value,
            bessel_term,
            precision);
    acb_mul(intermediate_function_value,
            intermediate_function_value,
            second_exponential_term,
            precision);

    // Memory Management - clear up the declared acbs
    acb_clear(prefactor);
    acb_clear(first_exponential_term);
    acb_clear(second_exponential_term);
    acb_clear(bessel_term);
    acb_clear(zero);
    acb_clear(two);
    acb_clear(minimum_phase_acb);
}

// Function computes the value of the integrand being integrated in the
// amplification factor calculation. The order parameter is unused but is
// required by the integration process.
int NfwIntegrand(acb_ptr integrand,
                 const acb_t integration_parameter,
                 void * parameter_set,
                 slong order,
                 slong precision) {
    // The parameter_set contains a vector which itself contains the
    // dimensionless frequency, source position, and scaling constant. These
    // need to be extracted and then placed into acb types for the rest of the
    // calculation
    std::vector<double> parameter_vector = (
        (std::vector<double> *) parameter_set)[0];
    double dimensionless_frequency_value = parameter_vector[0];
    double source_position_value = parameter_vector[1];
    double scaling_constant_value = parameter_vector[2];
    double minimum_phase = parameter_vector[3];

    acb_t dimensionless_frequency;
    acb_t source_position;
    acb_t scaling_constant;

    acb_init(dimensionless_frequency);
    acb_init(source_position);
    acb_init(scaling_constant);

    acb_set_d(dimensionless_frequency, dimensionless_frequency_value);
    acb_set_d(source_position, source_position_value);
    acb_set_d(scaling_constant, scaling_constant_value);

    // The integrand is a combination of two terms - the first is the
    // intermediate function k(w,y,x,ks) and the second is exp(iwx).
    // We must first then construct each of these terms and then
    // multiply them
    acb_t intermediate_function_term;
    acb_t exponential_term;

    acb_init(intermediate_function_term);
    acb_init(exponential_term);

    // Calculation of the k function term
    IntermediateFunctionCalculation(intermediate_function_term,
                                    dimensionless_frequency,
                                    source_position,
                                    integration_parameter,
                                    scaling_constant,
                                    minimum_phase,
                                    precision);

    // Calculation of the exponential term
    acb_mul(exponential_term,
            dimensionless_frequency,
            integration_parameter,
            precision);
    acb_mul_onei(exponential_term, exponential_term);
    acb_exp(exponential_term, exponential_term, precision);

    // Multiplication of the terms to the resultant acb at integrand
    acb_mul(
        integrand, intermediate_function_term, exponential_term, precision);

    // Memory Management - clear up the declared acbs inside the fucntion
    acb_clear(dimensionless_frequency);
    acb_clear(source_position);
    acb_clear(scaling_constant);
    acb_clear(intermediate_function_term);
    acb_clear(exponential_term);

    // Now return a zero to indicate function's successful completion
    return 0;
}

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


// Function computes the value of the first correction term for the
// amplification factor. The first correction term is given by
// -((k(w,y,x_upper_limit,ks) * exp(iw*x_upper_limit))/iw)
void FirstCorrectionTerm(acb_t first_correction_term,
                         acb_t dimensionless_frequency,
                         acb_t source_position,
                         acb_t integration_upper_limit,
                         acb_t scaling_constant,
                         double minimum_phase,
                         slong precision) {
    // The function is calculated by splitting the calculation into three parts
    // an intermediate function term, an exponential term, and the denominator
    // term. These must be initialised as acbs
    acb_t intermediate_function_term;
    acb_t exponential_term;
    acb_t denominator_term;

    acb_init(intermediate_function_term);
    acb_init(exponential_term);
    acb_init(denominator_term);

    // Calculate the intermediate function term
    IntermediateFunctionCalculation(intermediate_function_term,
                                    dimensionless_frequency,
                                    source_position,
                                    integration_upper_limit,
                                    scaling_constant,
                                    minimum_phase,
                                    precision);

    // Calculate the exponential term
    acb_mul(exponential_term,
            dimensionless_frequency,
            integration_upper_limit,
            precision);
    acb_mul_onei(exponential_term, exponential_term);
    acb_exp(exponential_term, exponential_term, precision);

    // Calculate the denominator term
    acb_mul_onei(denominator_term, dimensionless_frequency);

    // Construct the correction term value
    acb_mul(first_correction_term,
            intermediate_function_term,
            exponential_term,
            precision);
    acb_div(first_correction_term,
            first_correction_term,
            denominator_term,
            precision);
    acb_neg(first_correction_term, first_correction_term);

    // Memory Management - clear the declared acbs inside the function
    acb_clear(intermediate_function_term);
    acb_clear(exponential_term);
    acb_clear(denominator_term);
}


// Function computes the value of the second correction term for the
// amplification factor. This term is given by
// d(k(w,y,x,ks)*exp(iwz))/dx/(iw)^2
void SecondCorrectionTerm(acb_t second_correction_term,
                          acb_t dimensionless_frequency,
                          acb_t source_position,
                          acb_t integration_upper_limit,
                          acb_t scaling_constant,
                          double minimum_phase,
                          slong precision) {
    // The function will be computed by constructing three terms, the
    // derivative term, the exponential term, and the denominator term.
    // To first calculate the derivative term, a central finite differences
    // method is being applied with a set step-size of 0.00001. This means that
    // df/dx = f(x+h)-f(x-h)/2h.
    // Firstly we must initialise the necessary values of h, 2h, x+h, x-h,
    // f(x+h), f(x-h), and df/dx
    acb_t stepsize;
    acb_t upper_limit_plus;
    acb_t upper_limit_minus;
    acb_t double_stepsize;
    acb_t function_value_plus;
    acb_t function_value_minus;
    acb_t derivative_term;

    acb_init(stepsize);
    acb_init(upper_limit_plus);
    acb_init(upper_limit_minus);
    acb_init(double_stepsize);
    acb_init(function_value_plus);
    acb_init(function_value_minus);
    acb_init(derivative_term);

    // Calculate the values of the easily constructed values
    acb_set_d(stepsize, 0.00001);
    acb_add(double_stepsize, stepsize, stepsize, precision);
    acb_add(upper_limit_plus, integration_upper_limit, stepsize, precision);
    acb_sub(upper_limit_minus, integration_upper_limit, stepsize, precision);

    // Calculate the two function values
    IntermediateFunctionCalculation(function_value_plus,
                                    dimensionless_frequency,
                                    source_position,
                                    upper_limit_plus,
                                    scaling_constant,
                                    minimum_phase,
                                    precision);
    IntermediateFunctionCalculation(function_value_minus,
                                    dimensionless_frequency,
                                    source_position,
                                    upper_limit_minus,
                                    scaling_constant,
                                    minimum_phase,
                                    precision);

    // Now calculate the derivative term by using the central differences
    // method
    acb_sub(
        derivative_term, function_value_plus, function_value_minus, precision);
    acb_div(derivative_term, derivative_term, double_stepsize, precision);

    // Calculate the exponential term
    acb_t exponential_term;
    acb_init(exponential_term);

    acb_mul(exponential_term,
            dimensionless_frequency,
            integration_upper_limit,
            precision);
    acb_mul_onei(exponential_term, exponential_term);
    acb_exp(exponential_term, exponential_term, precision);

    // Calculate the denominator
    acb_t denominator_term;
    acb_init(denominator_term);

    acb_mul_onei(denominator_term, dimensionless_frequency);
    acb_sqr(denominator_term, denominator_term, precision);

    // Construct the full correction term
    acb_mul(second_correction_term,
            derivative_term,
            exponential_term,
            precision);
    acb_div(second_correction_term,
            second_correction_term,
            denominator_term,
            precision);

    // Memory Management - clear the declared acbs
    acb_clear(stepsize);
    acb_clear(upper_limit_plus);
    acb_clear(upper_limit_minus);
    acb_clear(double_stepsize);
    acb_clear(function_value_plus);
    acb_clear(function_value_minus);
    acb_clear(derivative_term);
    acb_clear(exponential_term);
    acb_clear(denominator_term);
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
                       NfwIntegrand,
                       &parameter_set,
                       lower_limit,
                       upper_limit,
                       goal,
                       tol,
                       options,
                       goal);

    acb_t w;
    acb_t y;
    acb_t ks;

    acb_init(w);
    acb_init(y);
    acb_init(ks);

    acb_set_d(w, dimensionless_frequency);
    acb_set_d(y, source_position);
    acb_set_d(ks, scaling_constant);

    // Caclulate the first correction term
    acb_t first_correction_term;
    acb_init(first_correction_term);
    FirstCorrectionTerm(first_correction_term,
                        w,
                        y,
                        upper_limit,
                        ks,
                        phase,
                        precision);

    // Calculate the second correction term
    acb_t second_correction_term;
    acb_init(second_correction_term);
    SecondCorrectionTerm(second_correction_term,
                         w,
                         y,
                         upper_limit,
                         ks,
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
    acb_clear(w);
    acb_clear(y);
    acb_clear(ks);
}



