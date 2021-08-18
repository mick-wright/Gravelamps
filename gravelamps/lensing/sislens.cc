// Program definition files for calculating the amplification factor for the
// Singular Isothermal Sphere (SIS) lens and the main loop for performing the
// calculation based upon these
//
// Mick Wright

#include "sis.h"

// Function computes the value of the intermediate function k(w,y,z) for the
// amplification factor calculation. The function k is given by
// -iw*exp(iw(y^2/2 + phi(y)))*J0(wy*sqrt(z))*exp(-iw*sqrt(2z))
void IntermediateFunctionCalculation(acb_t intermediate_function_value,
                                     acb_t dimensionless_frequency,
                                     acb_t source_position,
                                     const acb_t integration_parameter,
                                     slong precision) {
    // Initialise the components of the function - as can be seen there are
    // four - a prefactor, the first exponential term, the bessel term, and
    // finally, the second exponential term
    acb_t prefactor;
    acb_t first_exponential_term;
    acb_t bessel_term;
    acb_t second_exponential_term;

    acb_init(prefactor);
    acb_init(first_exponential_term);
    acb_init(bessel_term);
    acb_init(second_exponential_term);

    // In addition, we need a zero and a two value for use in the calculations
    acb_t two;
    acb_init(two);
    acb_set_d(two, 2);

    acb_t zero;
    acb_init(zero);
    acb_zero(zero);

    // Construct the prefactor term
    acb_mul_onei(prefactor, dimensionless_frequency);
    acb_neg(prefactor, prefactor);

    // Calculate the value of phi - y + 1/2
    acb_t phi;
    acb_init(phi);
    acb_one(phi);
    acb_div(phi, phi, two, precision);
    acb_add(phi, phi, source_position, precision);

    // Construct the first exponential term
    acb_sqr(first_exponential_term, source_position, precision);
    acb_div(first_exponential_term, first_exponential_term, two, precision);
    acb_add(first_exponential_term, first_exponential_term, phi, precision);
    acb_mul(first_exponential_term,
            dimensionless_frequency,
            first_exponential_term,
            precision);
    acb_mul_onei(first_exponential_term, first_exponential_term);
    acb_exp(first_exponential_term, first_exponential_term, precision);

    // Construct the Bessel Term
    acb_mul(bessel_term, two, integration_parameter, precision);
    acb_sqrt(bessel_term, bessel_term, precision);
    acb_mul(bessel_term, bessel_term, source_position, precision);
    acb_mul(bessel_term, bessel_term, dimensionless_frequency, precision);
    acb_hypgeom_bessel_j(bessel_term, zero, bessel_term, precision);

    // Construct the second term
    acb_mul(second_exponential_term, two, integration_parameter, precision);
    acb_sqrt(second_exponential_term, second_exponential_term, precision);
    acb_mul(second_exponential_term,
            second_exponential_term,
            dimensionless_frequency,
            precision);
    acb_mul_onei(second_exponential_term, second_exponential_term);
    acb_neg(second_exponential_term, second_exponential_term);
    acb_exp(second_exponential_term, second_exponential_term, precision);

    // Construct the value of the intermediate function by multiplication
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
    acb_clear(phi);
    acb_clear(prefactor);
    acb_clear(first_exponential_term);
    acb_clear(bessel_term);
    acb_clear(second_exponential_term);
    acb_clear(zero);
    acb_clear(two);
}

// Function computes the value of the integrand being integrated in the
// amplification factor calculation. The order parameter is unused but is
// required by the integration process.
int SisIntegrand(acb_ptr integrand,
                 const acb_t integration_parameter,
                 void * parameter_set,
                 slong order,
                 slong precision) {
    // The parameter set contains a vector which itself contains the
    // dimensionless frequency and source position. These need to be extracted
    // and then placed into acb types for the rest of the calculation
    std::vector<double> parameter_vector = (
        (std::vector<double> *) parameter_set)[0];
    double dimensionless_frequency_value = parameter_vector[0];
    double source_position_value = parameter_vector[1];

    acb_t dimensionless_frequency;
    acb_t source_position;

    acb_init(dimensionless_frequency);
    acb_init(source_position);

    acb_set_d(dimensionless_frequency, dimensionless_frequency_value);
    acb_set_d(source_position, source_position_value);

    // The integrand is a combination of two terms - the first is the
    // intermediate function k(w,y,x) and the second is exp(iwx). We must first
    // then construct each of these terms and then multiply them
    acb_t intermediate_function_term;
    acb_t exponential_term;

    acb_init(intermediate_function_term);
    acb_init(exponential_term);

    // Calculation of the k function term
    IntermediateFunctionCalculation(intermediate_function_term,
                                    dimensionless_frequency,
                                    source_position,
                                    integration_parameter,
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

    // Memory Management - clear up the declared acbs inside the function
    acb_clear(dimensionless_frequency);
    acb_clear(source_position);
    acb_clear(intermediate_function_term);
    acb_clear(exponential_term);

    // Now return a zero to indicate function's successful completion
    return 0;
}

// Function computes the value of the first correction term for the
// amplification factor. The first correction term is given by
// -((k(w,y,x_upper_limit) * exp(iw*x_upper_limit))/iw)
void FirstCorrectionTerm(acb_t first_correction_term,
                         acb_t dimensionless_frequency,
                         acb_t source_position,
                         acb_t integration_upper_limit,
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
                                    precision);

    // Calculate the exponential term
    acb_mul(exponential_term,
            dimensionless_frequency,
            integration_upper_limit,
            precision);
    acb_mul_onei(exponential_term, exponential_term);
    acb_exp(exponential_term, exponential_term, precision);

    // Calculate the denomiantor term
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
// d(k(w,y,x)*exp(iwz)/dx/(iw)^2
void SecondCorrectionTerm(acb_t second_correction_term,
                          acb_t dimensionless_frequency,
                          acb_t source_position,
                          acb_t integration_upper_limit,
                          slong precision) {
    // The function will be computed by constructing three terms, the
    // derivative term, the exponential term, and the denominator term.
    // To first calculate the derivative term, a central finite differences
    // method will be applied with a set step-size of 0.00001. This means that
    // df/dx = f(x+h)-f(x-h)/2h.
    // Firstly we must initialise the necessary values of h, 2h, x+-h, f(x+-h),
    // and df/dx
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
                                    precision);
    IntermediateFunctionCalculation(function_value_minus,
                                    dimensionless_frequency,
                                    source_position,
                                    upper_limit_minus,
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

// Function computes the amplification factor for an axially symmetric Singular
// Isothermal Sphere (SIS) lens using full wave optics for given values of
// dimensionless frequency and source position with arithmetic precision given
// by precision. The infinite integral is approximated by calculating the
// finite integral with upper limit given by integration_upper_limit.
void AmplificationFactorCalculation(acb_t amplification_factor,
                                    double dimensionless_frequency,
                                    double source_position,
                                    double integration_upper_limit,
                                    slong precision) {
    // The integrand function requires that the lensing parameters be passed to
    // it in the form of a single vector
    std::vector<double> parameter_set {dimensionless_frequency,
                                       source_position};

    // Construct abs for the lensing parameters to be used in the more complex
    // calculations
    acb_t dimensionless_frequency_acb;
    acb_t source_position_acb;

    acb_init(dimensionless_frequency_acb);
    acb_init(source_position_acb);

    acb_set_d(dimensionless_frequency_acb, dimensionless_frequency);
    acb_set_d(source_position_acb, source_position);

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

    // Create acbs for the lower and upper limits of the integration. The lower
    // limit of the integration should be zero, however, the calculation is
    // unable to process this, and so to avoid issues, the lower limit is set
    // to be 0.000001
    acb_t lower_limit;
    acb_t upper_limit;

    acb_init(lower_limit);
    acb_init(upper_limit);

    acb_set_d(lower_limit, 0.000001);
    acb_set_d(upper_limit, integration_upper_limit);

    // Create and calculate the integration function using the integrand
    // function
    acb_t integration_term;
    acb_init(integration_term);
    acb_calc_integrate(integration_term,
                       SisIntegrand,
                       &parameter_set,
                       lower_limit,
                       upper_limit,
                       goal,
                       tolerance,
                       integration_options,
                       precision);

    // Calculate the first correction term
    acb_t first_correction_term;
    acb_init(first_correction_term);
    FirstCorrectionTerm(first_correction_term,
                        dimensionless_frequency_acb,
                        source_position_acb,
                        upper_limit,
                        precision);

    // Calculate the second correction term
    acb_t second_correction_term;
    acb_init(second_correction_term);
    SecondCorrectionTerm(second_correction_term,
                         dimensionless_frequency_acb,
                         source_position_acb,
                         upper_limit,
                         precision);

    // Fromn the integration term, and the two correciton terms, construct the
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
    acb_clear(dimensionless_frequency_acb);
    acb_clear(source_position_acb);
    acb_clear(lower_limit);
    acb_clear(upper_limit);
    acb_clear(integration_term);
    acb_clear(first_correction_term);
    acb_clear(second_correction_term);
}

// Function computes the amplification factor for an axially symmetric singular
// isothermal sphere (SIS) style lens using the geometric optics approximation
// for given values of dimensionless frequency and source position
std::complex<double> AmplificationFactorGeometric(
    double dimensionless_frequency, double source_position) {
    // Using the complex_literals i operator
    using std::literals::complex_literals::operator""i;

    // Compute the positive magnification and square root
    double mag_plus = sqrt(abs(1. + 1./source_position));

    // If the source position is greater than 1, return the root of the
    // positive magnification
    if (source_position >= 1) {
        std::complex<double> geometric_factor = mag_plus;
        return geometric_factor;
    }

    // Compute the negative magnification and square root
    double mag_minus = sqrt(abs(1. - 1./source_position));

    // Calculate the time delay
    double time_delay = 2 * source_position;

    // Calculate the exponential term
    std::complex<double> exponent = 1i * dimensionless_frequency * time_delay;
    std::complex<double> exp_term = exp(exponent);

    // Construct the final amplification factor
    std::complex<double> geometric_factor =
        mag_plus - 1i * mag_minus * exp_term;

    return geometric_factor;
}

// Function constructs two matrices containing the real and imaginary parts of
// the value of the amplification factor function based upon two vectors
// containing values of dimensionless frequency and source position and
// returns these inside of a pair object
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
    AmplificationFactorMatrices(std::vector<double> dimensionless_frequency,
                                std::vector<double> source_position,
                                double integration_upper_limit,
                                slong precision,
                                slong approx_switch) {
    // Calculate the sizes of the dimensionless frequency and imapct parameter
    // vectors, and then construct two correctly size matrices
    int source_position_size = source_position.size();
    int dimensionless_frequency_size = dimensionless_frequency.size();

    std::vector<std::vector<double>> amp_fac_real(
        source_position_size,
        std::vector<double>(dimensionless_frequency_size));
    std::vector<std::vector<double>> amp_fac_imag(
        source_position_size,
        std::vector<double>(dimensionless_frequency_size));

    // Create the loop that goes through and calculates the amplification
    // factor values over which they will be calculated using full wave
    // optics calculations. The scheudule here is dynanmic due to the fact
    // that the calculation takes differing amounts of time at differring
    // points

    acb_t amplification_factor;
    acb_init(amplification_factor);

    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i=0; i < source_position_size; i++) {
        for (int j=0; j < approx_switch; j++) {
            AmplificationFactorCalculation(amplification_factor,
                                           dimensionless_frequency[j],
                                           source_position[i],
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
        for (int i=0; i < source_position_size; i++) {
            for (int j=approx_switch; j < dimensionless_frequency_size; j++) {
                std::complex<double> geometric_factor;
                geometric_factor =
                    AmplificationFactorGeometric(dimensionless_frequency[j],
                                                 source_position[i]);

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

// Main Function - this takes in seven arguments:
//     dimensionless_frequency_file - input file containing dimensionless
//                                    frequency values
//     source_position_file - input file containing source position values
//     amplification_factor_real_file - output file containing the real parts of
//                                      the amplification factor values
//     amplification_factor_imag_file - output file containing the imaginary
//                                      parts of the amplification factor values
//     integration_upper_limit - value to calculate the integral up to
//     precision - integer value used as the arithmetic preciusion in the
//                 amplification factor calculations
//     approx_switch - integer value which is the location in the dimensionless
//                     frequency array where the geometric optics approximation
//                     should take over
//
// Function takes in a dimensionless frequency and source position file each
// containing a vector of numbers which it will then generate a pair of matrices
// of amplification factor values from these, outputting the real and imaginary
// components of this to two files. The infinite summation is approximated by
// using a finite threshold value and the arithmetic is done to a speciofed
// precision value
int main(int argc, char* argv[]) {
    // Read in all of the filenames
    std::string dimensionless_frequency_file = argv[1];
    std::string source_position_file = argv[2];
    std::string amplification_factor_real_file = argv[3];
    std::string amplification_factor_imag_file = argv[4];

    // Read in the summation threshold and arithmetic precision values
    double integration_upper_limit = atof(argv[5]);
    slong precision = atoi(argv[6]);

    // Read in the position for the geometric optics switch
    slong approx_switch = atoi(argv[7]);

    // For each of the dimensionless frequency and source position files
    // open the file for reading, generate an iterator object which will get
    // all of the enclosed values and put them inside of a vector
    std::ifstream dim_freq_fstream(dimensionless_frequency_file);
    std::istream_iterator<double> dim_freq_start(dim_freq_fstream),
                                  dim_freq_end;
    std::vector<double> dimensionless_frequency(dim_freq_start, dim_freq_end);

    std::ifstream sour_pos_fstream(source_position_file);
    std::istream_iterator<double> sour_pos_start(sour_pos_fstream), sour_pos_end;
    std::vector<double> source_position(sour_pos_start, sour_pos_end);

    // With the vectors generated, now perform the main loop of calculating
    // the amplification factor values and return these as a pair of matrices
    std::pair<std::vector<std::vector<double>>,
              std::vector<std::vector<double>>> amp_fac_matrices;
    amp_fac_matrices = AmplificationFactorMatrices(
        dimensionless_frequency, source_position,
        integration_upper_limit, precision, approx_switch);

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
