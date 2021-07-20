// Program definition files for calculating the amplification factor for the
// Navarro, Frenk, and White (NFW) lens and the main loop for performing the
// calculation based upon these
//
// Mick Wright 2021

#include "nfw.h"

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

std::complex<double> LensingPotential(double scaled_surface_density,
                                      double scaling_constant) {
    // If the scaled surface density is 1 - return 0
    if (scaled_surface_density == 1) {
        return 0;
    }

    // create a complex number containing the scaled surface density
    std::complex<double> scaled_complex = scaled_surface_density;

    // Otherwise, compute the prefactor and log term
    double prefactor = scaling_constant/2.;
    std::complex<double> log_term = log(scaled_complex/2.);
    log_term = log_term * log_term;

    std::complex<double> trig_term;
    // Compute the trig term based upon whether the scaled surface density is
    // greater or less than one
    if (scaled_surface_density < 1) {
        trig_term = 1. - scaled_surface_density * scaled_surface_density;
        trig_term = -1. * atanh(sqrt(trig_term)) * atanh(sqrt(trig_term));
    } else {
        trig_term = scaled_surface_density * scaled_surface_density - 1.;
        trig_term = atan(sqrt(trig_term)) * atan(sqrt(trig_term));
    }

    // Construct the final result
    std::complex<double> lensing_potential = prefactor * (log_term + trig_term);
    return lensing_potential;
}

// Function computes the first two terms of the time delay function for finding
// the minimum which yields the phase guage needed
std::complex<double> TimeDelayPartial(double scaled_surface_density,
                                      double impact_parameter,
                                      double scaling_constant) {
    // Compute the power term, (x-y)**2/2
    double power_term = scaled_surface_density - impact_parameter;
    power_term = (power_term * power_term)/2.;

    // Calculate the psi term
    std::complex<double> psi_term = LensingPotential(scaled_surface_density,
                                                     scaling_constant);

    // Compute the final result
    std::complex<double> partial_time_delay = power_term - psi_term;

    return partial_time_delay;
}

// Function computes the minimum time delay phase
std::complex<double> MinTimeDelayPhase(double impact_parameter, double scaling_constant) {
    // Starting at zero, look for the first minimum in a step size of 0.0001
    double test_value = 0;
    double step_size = 0.0001;
    std::complex<double> time_delay = 0;
    std::complex<double> previous_time_delay = 1000000;

    // Loop through, calculating the partial time delay and checking against
    // previous value to see if we've passed a minimum. If so, return the
    // negative. If not, continue on
    while (true) {
        time_delay = TimeDelayPartial(
            test_value, impact_parameter, scaling_constant);

        if (std::real(time_delay) > std::real(previous_time_delay)) {
             return -1.*previous_time_delay;
        }

        previous_time_delay = time_delay;
        test_value += step_size;
    }
}

// Function computes the lens equation
double LensEquation(double scaled_surface_density,
                    double impact_parameter,
                    double scaling_constant) {
    // Compute the real part of the derivative of the value of psi(x,ks). This
    // is computed using a central finite differences method with step size
    // 0.00001
    double step_size = 0.00001;
    double scaled_minus = scaled_surface_density - step_size;
    double scaled_plus = scaled_surface_density + step_size;
    std::complex<double> psi_minus = LensingPotential(
        scaled_minus, scaling_constant);
    std::complex<double> psi_plus = LensingPotential(
        scaled_plus, scaling_constant);
    double differentiated_psi = std::real(
        (psi_plus - psi_minus)/(2 * step_size));

    // Compute the value of the lensing equation
    double lens_eq_value = scaled_surface_density - differentiated_psi;
    lens_eq_value -= impact_parameter;

    return lens_eq_value;
}

// Function computes the image positions for a given impact parameter
// by means of root finding the lens equation in the three regimes
std::vector<double> ImagePositions(double impact_parameter,
                                   double scaling_constant) {
    using boost::math::tools::bisect;
    using boost::math::tools::eps_tolerance;

    // Setting the desired accuracy as three less digits than the maximal
    // precision of the double type. This is to give good accuracy whilst
    // preventing thrashing around of the iterations at the end
    int digits = std::numeric_limits<double>::digits;
    int get_digits = digits - 3;
    eps_tolerance<double> tolerance(get_digits);

    // Create lambda expression version of the LensEquation to make it
    // one parameter only
    auto LensEq =
        [impact_parameter, scaling_constant](double scaled_surface_density) {
            return LensEquation(scaled_surface_density,
                                impact_parameter,
                                scaling_constant);
        };

    // Find the three possible roots - these are a negative root, one around
    // zero, and one positive root. If y is greater than the critical value
    // there will only be the positive root and the others will fail. If so
    // set them to the same value
    std::pair<double, double> root_pos_pair = bisect(
        LensEq, 0.1, 12., tolerance);
    double root_pos = root_pos_pair.first + (
        root_pos_pair.second - root_pos_pair.first)/2.;

    double root_neg = root_pos;
    double root_zero = root_pos;

    try {
        std::pair<double, double> root_neg_pair = bisect(
            LensEq, -12., -0.1, tolerance);
        root_neg = root_neg_pair.first + (
           root_neg_pair.second - root_neg_pair.first)/2.;
    } catch(...) {
    }

    try {
        std::pair<double, double> root_zero_pair = bisect(
            LensEq, -0.1, 0.1, tolerance);
        root_zero = root_zero_pair.first + (
            root_zero_pair.second - root_zero_pair.first)/2.;
    } catch(...) {
    }

    // Create the vector containing the roots - assuming that none of the roots
    // have been returned as the same value
    std::vector<double> image_positions = {root_pos};

    if (root_neg != root_pos) {
        image_positions.push_back(root_neg);
    }

    if (root_zero != root_pos && root_zero != root_neg) {
        image_positions.push_back(root_zero);
    }

    return image_positions;
}

// Function computes the time delay for a given image position and impact
// parameter and minimum phase delay. This is the complete version of the above
// partial function
std::complex<double> TimeDelay(double image_position,
                               double impact_parameter,
                               double scaling_constant,
                               double phase_minimum) {
    // Compute the value of the lensing potential for the given image position
    std::complex<double> lensing_potential = LensingPotential(
        abs(image_position), scaling_constant);

    // Compute the first term 0.5 * (x-y)**2
    double power_term = abs(image_position - impact_parameter);
    power_term = (power_term * power_term)/2.;

    // Construct the final time delay
    std::complex<double> time_delay =
        power_term - lensing_potential + phase_minimum;

    return time_delay;
}

// Function computes the NFW mass for a given image position and scaling
// constant
std::complex<double> ImagePositionMass(double image_position,
                                       double scaling_constant) {
    // The Mass is given as the combination of a trig term and a log term.
    // First calculate the trig term which is dependent upon the vakye if the
    // image position
    std::complex<double> trig_term;
    std::complex<double> image_complex = image_position;
    std::complex<double> image_position_term;

    // The three cases are x = 1, x < 1, or x > 1. In the case of the first
    // the solution is the most simple, the trig term is also 1. In the case
    // where x < 1, the trig term is arctanh(sqrt(1-x^2))/sqrt(1-x^2). In the
    // case of x > 1 the trig term is arctan(sqrt(x^2-1))/sqrt(x^2-1).
    if (image_position == 1) {
        trig_term = 1.;
    } else if (image_position < 1) {
        image_position_term = sqrt(1. - image_complex * image_complex);
        trig_term = atanh(image_position_term)/image_position_term;
    } else {
        image_position_term = sqrt(image_complex * image_complex - 1.);
        trig_term = atan(image_position_term)/image_position_term;
    }

    // The logarithimic term is log(x/2)
    std::complex<double> log_term = log(image_complex/2.);

    // Compute the final result - ks * (log_term + trig_term)
    std::complex<double> image_position_mass =
        scaling_constant * (log_term + trig_term);

    return image_position_mass;
}

// Function computes the surface density for a given image position and
// scaling constant
std::complex<double> SurfaceDensity(double image_position,
                                    double scaling_constant) {
    // Compute the prefactor - ks/2.
    double prefactor = scaling_constant/2.;

    // Calculate the interior term. This is dependent upon the image position
    // value. If x = 1, then the interior term is just 1/3. If x > 1, then the
    // value is given by -2/((x^2 - 1)^3/2) * atan(sqrt((x-1)/(1-x))) *
    // 1/(x^2 - 1). If x < 1, then the value is given by 2/((1-x^2)^3/2) *
    // arctanh(sqrt((1-x)/(x-1))) * 1/(1 - x^2)
    std::complex<double> interior_term;
    std::complex<double> image_complex = image_position;
    std::complex<double> interior_prefactor;
    std::complex<double> trig_term;
    std::complex<double> final_term;

    if (image_position == 1) {
        interior_term = 1./3.;
    } else if (image_position > 1) {
        interior_prefactor = image_complex*image_complex - 1.;
        interior_prefactor = pow(interior_prefactor, 1.5);
        interior_prefactor = -2./interior_prefactor;

        trig_term = (image_complex - 1.)/(1. + image_complex);
        trig_term = atan(sqrt(trig_term));

        final_term = 1./(image_complex*image_complex - 1.);

        interior_term = interior_prefactor * trig_term + final_term;
    } else {
       interior_prefactor = 1. - image_complex*image_complex;
       interior_prefactor = pow(interior_prefactor, 1.5);
       interior_prefactor = 2./interior_prefactor;

       trig_term = (1. - image_complex)/(image_complex + 1.);
       trig_term = atanh(sqrt(trig_term));

       final_term = 1./(1. - image_complex*image_complex);

       interior_term = interior_prefactor * trig_term - final_term;
    }

    // Construct the surface density by multiplying the interior term by the
    // prefactor
    std::complex<double> surface_density = prefactor * interior_term;

    return surface_density;
}

// Function computes the determinant of the matrix used to calculate the
// magnification
std::complex<double> MatrixDeterminant(double image_position,
                                       double scaling_constant) {
    // Compute the first term 1 - ImagePositionMass(abs(x), ks)/x^2
    std::complex<double> first_term = ImagePositionMass(
        abs(image_position), scaling_constant)/(image_position*image_position);
    first_term = 1. - first_term;

    // Compute the first bracket term - ImagePositionMass(abs(x,ks)/abs(x)^2
    std::complex<double> bracket_term_one = ImagePositionMass(
        abs(image_position), scaling_constant)/(
            abs(image_position) * abs(image_position));

    // Compute the second bracket term - 2 * SurfaceDensity(abs(x), ks)
    std::complex<double> bracket_term_two = 2. * SurfaceDensity(
        abs(image_position), scaling_constant);

    // Compute the full bracket term (1. + term_one - term_two)
    std::complex<double> bracket_term =
        1. + bracket_term_one - bracket_term_two;

    // Compute the Determinant as the first term * bracket_term
    std::complex<double> determinant = first_term * bracket_term;

    return determinant;
}

// Function computes the magnification for a given image position and scaling
// constant
std::complex<double> Magnification(double image_position,
                                   double scaling_constant) {
    // The magnification is simply 1/MatrixDeterminant(x,ks)
    return 1./MatrixDeterminant(image_position, scaling_constant);
}

// Function computes the value of the intermediate function k(w,y,z,ks) for the
// amplification factor calculation. The function k is given by
// -iw*(exp[iw(y^2/2 + phimin)]*J0(wy*sqrt(2x))*exp(-iw*psi(sqrt(2x),ks)))
void IntermediateFunctionCalculation(acb_t intermediate_function_value,
                                     acb_t dimensionless_frequency,
                                     acb_t impact_parameter,
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
    acb_t zero;
    acb_t two;

    acb_init(zero);
    acb_init(two);

    acb_zero(zero);
    acb_set_d(two, 2);

    // Finally we need to calculate phimin - the phase required for a minimum
    // time delay of zero;
    acb_t minimum_phase_acb;
    acb_init(minimum_phase_acb);
    acb_set_d(minimum_phase_acb, minimum_phase);

    // Construct the prefactor term
    acb_mul_onei(prefactor, dimensionless_frequency);
    acb_neg(prefactor, prefactor);

    // Construct the first exponential term
    acb_sqr(first_exponential_term, impact_parameter, precision);
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

    // Construct the Bessel Term
    acb_mul(bessel_term, two, integration_parameter, precision);
    acb_sqrt(bessel_term, bessel_term, precision);
    acb_mul(bessel_term, bessel_term, impact_parameter, precision);
    acb_mul(bessel_term, bessel_term, dimensionless_frequency, precision);
    acb_hypgeom_bessel_j(bessel_term, zero, bessel_term, precision);

    // Construct the second exponential term
    acb_mul(second_exponential_term, two, integration_parameter, precision);
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
    acb_clear(bessel_term);
    acb_clear(second_exponential_term);
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
    // dimensionless frequency, impact parameter, and scaling constant. These
    // need to be extracted and then placed into acb types for the rest of the
    // calculation
    std::vector<double> parameter_vector = (
        (std::vector<double> *) parameter_set)[0];
    double dimensionless_frequency_value = parameter_vector[0];
    double impact_parameter_value = parameter_vector[1];
    double scaling_constant_value = parameter_vector[2];
    double minimum_phase = parameter_vector[3];

    acb_t dimensionless_frequency;
    acb_t impact_parameter;
    acb_t scaling_constant;

    acb_init(dimensionless_frequency);
    acb_init(impact_parameter);
    acb_init(scaling_constant);

    acb_set_d(dimensionless_frequency, dimensionless_frequency_value);
    acb_set_d(impact_parameter, impact_parameter_value);
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
                                    impact_parameter,
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
    acb_clear(impact_parameter);
    acb_clear(scaling_constant);
    acb_clear(intermediate_function_term);
    acb_clear(exponential_term);

    // Now return a zero to indicate function's successful completion
    return 0;
}

// Function computes the value of the first correction term for the
// amplification factor. The first correction term is given by
// -((k(w,y,x_upper_limit,ks) * exp(iw*x_upper_limit))/iw)
void FirstCorrectionTerm(acb_t first_correction_term,
                         acb_t dimensionless_frequency,
                         acb_t impact_parameter,
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
                                    impact_parameter,
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
                          acb_t impact_parameter,
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
                                    impact_parameter,
                                    upper_limit_plus,
                                    scaling_constant,
                                    minimum_phase,
                                    precision);
    IntermediateFunctionCalculation(function_value_minus,
                                    dimensionless_frequency,
                                    impact_parameter,
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

// Function computes the amplification factor for an axially symmetric Navarro,
// Frenk, and White (NFW) lens for given values of dimensionless frequency and
// impact parameter with arithmetic precision given by precision. The inifinte
// integral is approximated by calculating the finite integral with upper limit
// given by integration_upper_limit
void AmplificationFactorCalculation(acb_t amplification_factor,
                                    double dimensionless_frequency,
                                    double impact_parameter,
                                    double scaling_constant,
                                    double integration_upper_limit,
                                    slong precision) {
    // Calculate the minimum time delay phase
    double minimum_phase = std::real(MinTimeDelayPhase(impact_parameter,
                                                       scaling_constant));

    // The integrand function requires that the lensing parameters be passed to
    // it in the form of a single vector
    std::vector<double> parameter_set {dimensionless_frequency,
                                       impact_parameter,
                                       scaling_constant,
                                       minimum_phase};

    // Construct acbs for the lensing parameters to be used in the more complex
    // calculations
    acb_t dimensionless_frequency_acb;
    acb_t impact_parameter_acb;
    acb_t scaling_constant_acb;

    acb_init(dimensionless_frequency_acb);
    acb_init(impact_parameter_acb);
    acb_init(scaling_constant_acb);

    acb_set_d(dimensionless_frequency_acb, dimensionless_frequency);
    acb_set_d(impact_parameter_acb, impact_parameter);
    acb_set_d(scaling_constant_acb, scaling_constant);

    // Set the goal and the tolerance for the integration, these are based upon
    // the precision specified
    slong goal = precision;

    mag_t tolerance;
    mag_set_ui_2exp_si(tolerance, 1, -1*precision);

    // Integration options - these are set such that the precision changes the
    // depth limit as 128 times the precision, and the evaluation limit as the
    // cube of the precision
    acb_calc_integrate_opt_t integration_options;
    acb_calc_integrate_opt_init(integration_options);

    integration_options -> use_heap = 1;
    integration_options -> depth_limit = 128 * precision;
    integration_options -> eval_limit = precision*precision*precision;

    // Create acbs for the lower and upper limits of the integration. The
    // lower limit of the integration should be zero, however, the calculation
    // is unable to process this, and so to avoid issues, the lower limit is
    // set to be 0.000001.
    acb_t lower_limit;
    acb_t upper_limit;

    acb_init(lower_limit);
    acb_init(upper_limit);

    acb_set_d(lower_limit, 0.000001);
    acb_set_d(upper_limit, integration_upper_limit);

    // Create and calculate the integration term using the integrand function
    acb_t integration_term;
    acb_init(integration_term);
    acb_calc_integrate(integration_term,
                       NfwIntegrand,
                       &parameter_set,
                       lower_limit,
                       upper_limit,
                       goal,
                       tolerance,
                       integration_options,
                       precision);

    acb_printn(integration_term, 50, 0); flint_printf("\n");

    // Calculate the first correction term
    acb_t first_correction_term;
    acb_init(first_correction_term);
    FirstCorrectionTerm(first_correction_term,
                        dimensionless_frequency_acb,
                        impact_parameter_acb,
                        upper_limit,
                        scaling_constant_acb,
                        minimum_phase,
                        precision);

    // Calculate the second correction term
    acb_t second_correction_term;
    acb_init(second_correction_term);
    SecondCorrectionTerm(second_correction_term,
                         dimensionless_frequency_acb,
                         impact_parameter_acb,
                         upper_limit,
                         scaling_constant_acb,
                         minimum_phase,
                         precision);

    // From the integration term, and the two correction terms, construct the
    // final value of the amplification factor
    acb_add(amplification_factor,
            integration_term,
            first_correction_term,
            precision);
    acb_add(amplification_factor,
            amplification_factor,
            second_correction_term,
            precision);

    // Memory Management - clear the acbs declared in the function
    //
    acb_clear(dimensionless_frequency_acb);
    acb_clear(impact_parameter_acb);
    acb_clear(scaling_constant_acb);
    acb_clear(lower_limit);
    acb_clear(upper_limit);
    acb_clear(integration_term);
    acb_clear(first_correction_term);
    acb_clear(second_correction_term);
}

// Function computes the amplification factor for an axially symmetric Navarro,
// Frenk, White (NFW) style lens using the geometric optics approximation
// for given values of dimensionless frequency and imapct parameter
std::complex<double> AmplificationFactorGeometric(
    double dimensionless_frequency,
    double impact_parameter,
    double scaling_constant) {
    // Using the i operator
    using std::literals::complex_literals::operator""i;

    // First compute the Image Positions
    std::vector<double> image_positions = ImagePositions(impact_parameter,
                                                         scaling_constant);

    int number_of_images = image_positions.size();

    // Next compute the phase needed for zero minimum time delay
    double minimum_phase = std::real(MinTimeDelayPhase(impact_parameter,
                                                       scaling_constant));

    std::cout << minimum_phase << std::endl;

    // The amplification factor is constructed as the sum of
    // sqrt(u(xj, ks)) * exp(iw*Tj - I*pi*nj(xj, ks) for each image
    // position
    std::complex<double> geometric_factor;

    for (int i=0; i < number_of_images; i++) {
        // Compute the time delay
        std::complex<double> time_delay = TimeDelay(image_positions[i],
                                                    impact_parameter,
                                                    scaling_constant,
                                                    std::real(minimum_phase));;

        // Compute the Magnification
        std::complex<double> magnification = Magnification(image_positions[i],
                                                           scaling_constant);
	magnification = abs(magnification);

        // Calculate the morse phase for the image position. This is based on
        // the surface density for the position.
        std::complex<double> double_surface_density;
        double_surface_density = 2. * SurfaceDensity(abs(image_positions[i]),
                                                     scaling_constant);
        double testing_condition = std::real(double_surface_density);

        double morse_factor;
        if (testing_condition == 0) {
            morse_factor = 0.5;
        } else if (testing_condition > 0) {
            morse_factor = 0;
        } else {
            morse_factor = 1;
        }

        // With each of the parameters, we can now construct the amplification
        // factor contribution
        std::complex<double> factor_contribution;
        factor_contribution = (1i * dimensionless_frequency * time_delay)
                              - (1i * M_PI * morse_factor);
        factor_contribution = exp(factor_contribution);

        factor_contribution = factor_contribution * sqrt(magnification);

        geometric_factor += factor_contribution;
    }

    return geometric_factor;
}

// Function constructs two matrices containing the real and imaginary parts of
// the value of the amplification factor function based upon two vectors
// containing values of dimensionless frequency and impact parameter and
// returns these inside of a pair object
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
    AmplificationFactorMatrices(std::vector<double> dimensionless_frequency,
                                std::vector<double> impact_parameter,
                                double scaling_constant,
                                double integration_upper_limit,
                                slong precision,
                                slong approx_switch) {
    // Calculate the sizes of the dimensionless frequency and impact parameter
    // vectors, and then construct two correctly sized matrices
    int impact_parameter_size = impact_parameter.size();
    int dimensionless_frequency_size = dimensionless_frequency.size();

    std::vector<std::vector<double>> amp_fac_real(
        impact_parameter_size,
        std::vector<double>(dimensionless_frequency_size));
    std::vector<std::vector<double>> amp_fac_imag(
        impact_parameter_size,
        std::vector<double>(dimensionless_frequency_size));

    // Create the loop that goes through and calculates the amplification
    // factor values over which they will be calculated using full wave
    // optics calculations. The schedule here is dynamic due to the fact
    // that the calculation takes differing amounts of time at differing
    // points
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i=0; i < impact_parameter_size; i++) {
        for (int j=0; j < approx_switch; j++) {
            acb_t amplification_factor;
            acb_init(amplification_factor);
            AmplificationFactorCalculation(amplification_factor,
                                           dimensionless_frequency[j],
                                           impact_parameter[i],
                                           scaling_constant,
                                           integration_upper_limit,
                                           precision);
            amp_fac_real[i][j] = arf_get_d(
                arb_midref(acb_realref(amplification_factor)), ARF_RND_NEAR);
            amp_fac_imag[i][j] = arf_get_d(
                arb_midref(acb_imagref(amplification_factor)), ARF_RND_NEAR);
        }
    }

    // Create the loop that goes through and calculates the amplification
    // factor values over which they will be calculated using the geometric
    // optics approximation
    if (approx_switch < dimensionless_frequency_size) {
        #pragma omp parallel for collapse(2) schedule(dynamic)
        for (int i=0; i < impact_parameter_size; i++) {
            for (int j=approx_switch; j < dimensionless_frequency_size; j++) {
                std::complex<double> geometric_factor;
                geometric_factor =
                    AmplificationFactorGeometric(dimensionless_frequency[j],
                                                 impact_parameter[i],
                                                 scaling_constant);

                amp_fac_real[i][j] = std::real(geometric_factor);
                amp_fac_imag[i][j] = std::imag(geometric_factor);
            }
        }
    }

    // Put the two matrices inside of a pair object and return them
    std::pair<std::vector<std::vector<double>>,
              std::vector<std::vector<double>>> amplification_factor_matrices(
                  amp_fac_real, amp_fac_imag);

    return amplification_factor_matrices;
}

// Main Function - this takes in eight arguemnts:
//     dimensionless_frequency_file - input file containing dimensionless
//                                    frequency values
//     impact_parameter_file - input file containing impact parameter values
//     amplification_factor_real_file - output file containing the real parts of
//                                      the amplification factor values
//     amplification_factor_imag_file - output file containing the imaginary
//                                      parts of the amplification factor values
//     scaling_constant - value of the scaling constant (ks) for the profile
//     integration_upper_limit - value to calculate the integral up to
//     precision - integer value used as the arithmetic precision
//     approx_switch - integer value which is the location in the dimensionless
//                     frequency array where the geometric optics approximation
//                     should take over
//
// Function takes in a dimensionless frequency and impact parameter file each
// containing a vector of numbers which it will then generate a pair of matrices
// of amplification factor values from these, outputting the real and imaginary
// components of this to two files. The integration from zero to infinity is
// approximated by using a user specified finite upper limit and an
// infinitesimal but non-zero lower limit. The arithmetic is done to a
// specified precision.
int main(int argc, char* argv[]) {
    // Read in all of the filenames
    std::string dimensionless_frequency_file = argv[1];
    std::string impact_parameter_file = argv[2];
    std::string amplification_factor_real_file = argv[3];
    std::string amplification_factor_imag_file = argv[4];

    // Read in the scaling constant, integration upper limit and precision
    double scaling_constant = atof(argv[5]);
    double integration_upper_limit = atof(argv[6]);
    slong precision = atoi(argv[7]);

    // Read in the position for the geometric optics switch
    slong approx_switch = atoi(argv[8]);

    // For each of the dimensionless frequency and impact parameter files
    // open the file for reading, generate an iterator object which will get
    // all of the enclosed values and put them inside of a vector
    std::ifstream dim_freq_fstream(dimensionless_frequency_file);
    std::istream_iterator<double> dim_freq_start(dim_freq_fstream),
                                  dim_freq_end;
    std::vector<double> dimensionless_frequency(dim_freq_start, dim_freq_end);

    std::ifstream imp_par_fstream(impact_parameter_file);
    std::istream_iterator<double> imp_par_start(imp_par_fstream), imp_par_end;
    std::vector<double> impact_parameter(imp_par_start, imp_par_end);

    // With the vectors generated, now perform the main loop of calculating the
    // amplification factor values and return these as a pair of matrices
    std::pair<std::vector<std::vector<double>>,
              std::vector<std::vector<double>>> amp_fac_matrices;
    amp_fac_matrices = AmplificationFactorMatrices(dimensionless_frequency,
                                                   impact_parameter,
                                                   scaling_constant,
                                                   integration_upper_limit,
                                                   precision,
                                                   approx_switch);

    // Open the amplification factor files for writing
    std::ofstream amp_fac_real_fstream(amplification_factor_real_file);
    std::ofstream amp_fac_imag_fstream(amplification_factor_imag_file);

    // Loop over the matrices and write the values of the amplification factor
    // to the files
    int number_rows = amp_fac_matrices.first.size();
    int number_columns = amp_fac_matrices.first[0].size();

    for (int i=0; i < number_rows; i++) {
        for (int j=0; j < number_columns; j++) {
            amp_fac_real_fstream << amp_fac_matrices.first[i][j] << "\t";
            amp_fac_imag_fstream << amp_fac_matrices.second[i][j] << "\t";
        }
        amp_fac_real_fstream << "\n";
        amp_fac_imag_fstream << "\n";
    }
}
