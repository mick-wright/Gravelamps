// Program definition files for calculating the amplification factor for the
// Singular Isothermal Sphere (SIS) lens and the main loop for performing the
// calculation based upon these
//
// Mick Wright

#include "sis.h"

// Function computes the value of phi - the phase constant used to obtain the
// minimum time delay induced by the lensing. Function is given by:
// phi = y + 1/2
double MinTimeDelayPhaseConstant(double impact_parameter) {
    return impact_parameter + 1/2.;
}

// Function computes the amplification factor for an axially symmetric singular
// isothermal sphere (SIS) lens for given values of dimensionless frequency and
// impact parameter with arithmetic precision given by precision. The infinite
// sum is calculated to an upper value given by sum_threshold
void AmplificationFactorCalculation(acb_t amplification_factor,
                                    double dimensionless_frequency,
                                    double impact_parameter,
                                    int sum_threshold,
                                    slong precision) {
    // Calculate the value of the exponential prefactor of the expression
    double prefactor_dim_freq_term = dimensionless_frequency/2.;
    double prefactor_imp_par_term = impact_parameter * impact_parameter;
    double prefactor_phi_term = 2 * MinTimeDelayPhaseConstant(
        impact_parameter);
    double prefactor_imp_phi = prefactor_imp_par_term + prefactor_phi_term;
    double prefactor_exponent_imag =
        prefactor_dim_freq_term * prefactor_imp_phi;

    // Construct the prefactor term
    acb_t prefactor_term;
    acb_init(prefactor_term);
    acb_set_d_d(prefactor_term, 0, prefactor_exponent_imag);
    acb_exp(prefactor_term, prefactor_term, precision);

    // Initialise the Summation term
    acb_t summation_term;
    acb_init(summation_term);
    acb_zero(summation_term);

    // Loop from 0 to the value of the sum_threshold calculating the value of
    // the ongoing summation
    for (int n=0; n <= sum_threshold; n++) {
        // Construct the Gamma function component
        double gamma_argument_real = 1 + n/2.;

        acb_t gamma_component;
        acb_init(gamma_component);
        acb_set_d_d(gamma_component, gamma_argument_real, 0);
        acb_gamma(gamma_component, gamma_component, precision);

        // Construct the n factorial component
        arb_t n_factorial_real;
        arb_init(n_factorial_real);
        arb_fac_ui(n_factorial_real, n, precision);

        acb_t n_factorial;
        acb_init(n_factorial);
        acb_set_arb(n_factorial, n_factorial_real);

        // Calculate the full gamma term - gamma component/n factorial
        acb_div(gamma_component, gamma_component, n_factorial, precision);

        // Construct the Power term
        double power_term_prefactor_real = 2 * dimensionless_frequency;
        double power_term_power_real = n/2.;

        arb_t power_term_power;
        arb_init(power_term_power);
        arb_set_d(power_term_power, power_term_power_real);

        acb_t power_term;
        acb_init(power_term);
        acb_set_d_d(power_term, 0, -1*power_term_prefactor_real);
        acb_pow_arb(power_term, power_term, power_term_power, precision);

        // Construct the confluent hypergeometric function term
        double hyper_arg_z_imag =
            -0.5 * dimensionless_frequency * impact_parameter*impact_parameter;

        acb_t hyper_arg_a;
        acb_t hyper_arg_b;
        acb_t hyper_arg_z;

        acb_init(hyper_arg_a);
        acb_init(hyper_arg_b);
        acb_init(hyper_arg_z);

        acb_set_d_d(hyper_arg_a, gamma_argument_real, 0);
        acb_one(hyper_arg_b);
        acb_set_d_d(hyper_arg_z, 0, hyper_arg_z_imag);

        acb_t hyper_term;
        acb_init(hyper_term);
        acb_hypgeom_1f1(hyper_term,
                        hyper_arg_a,
                        hyper_arg_b,
                        hyper_arg_z,
                        0,
                        precision);

        // Construct the value to add to the sum
        acb_t sum_temp;
        acb_init(sum_temp);
        acb_mul(sum_temp, gamma_component, power_term, precision);
        acb_mul(sum_temp, sum_temp, hyper_term, precision);

        // Add the sum_temp value to the total sum
        acb_add(summation_term, summation_term, sum_temp, precision);

        // Clear the acb_ts for memory managment
        acb_clear(gamma_component);
        arb_clear(n_factorial_real);
        acb_clear(n_factorial);
        arb_clear(power_term_power);
        acb_clear(power_term);
        acb_clear(hyper_arg_a);
        acb_clear(hyper_arg_b);
        acb_clear(hyper_arg_z);
        acb_clear(hyper_term);
        acb_clear(sum_temp);
    }

    // Construct the final term by multiplying the prefactor and summation term
    acb_mul(amplification_factor, prefactor_term, summation_term, precision);

    // Clear the prefactor and summation term acbs
    acb_clear(prefactor_term);
    acb_clear(summation_term);
}

// Function constructs two matrices containing the real and imaginary parts of
// the value of the amplification factor function based upon two vectors
// containing values of dimensionless frequency and impact parameter and
// returns these inside of a pair object
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
    AmplificationFactorMatrices(std::vector<double> dimensionless_frequency,
                                std::vector<double> impact_parameter,
                                int sum_threshold,
                                slong precision) {
    // Calculate the sizes of the dimensionless frequency and imapct parameter
    // vectors, and then construct two correctly size matrices
    int impact_parameter_size = impact_parameter.size();
    int dimensionless_frequency_size = dimensionless_frequency.size();

    std::vector<std::vector<double>> amp_fac_real(
        impact_parameter_size,
        std::vector<double>(dimensionless_frequency_size));
    std::vector<std::vector<double>> amp_fac_imag(
        impact_parameter_size,
        std::vector<double>(dimensionless_frequency_size));

    // Construct the main loop - looping through the dimensionless frequency
    // and impact parameter values to calculate the amplification factor value
    // for the given combination and storing the real and imaginary parts
    // inside of the matrices. Scheduler is dynamic because the amount of time
    // necessary to calculate the value is not static over the range of values

    acb_t amplification_factor;
    acb_init(amplification_factor);

    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i=0; i < impact_parameter_size; i++) {
        for (int j=0; j < dimensionless_frequency_size; j++) {
            AmplificationFactorCalculation(amplification_factor,
                                           dimensionless_frequency[j],
                                           impact_parameter[i],
                                           sum_threshold,
                                           precision);
            amp_fac_real[i][j] = arf_get_d(
                arb_midref(acb_realref(amplification_factor)), ARF_RND_NEAR);
            amp_fac_imag[i][j] = arf_get_d(
                arb_midref(acb_imagref(amplification_factor)), ARF_RND_NEAR);
        }
    }

    // Put the two matrices inside of a pair object and return them
    std::pair<std::vector<std::vector<double>>,
              std::vector<std::vector<double>>> amplification_factor_matrices(
                  amp_fac_real, amp_fac_imag);

    return amplification_factor_matrices;
}

// Main Function - this takes in six arguments:
//     dimensionless_frequency_file - input file containing dimensionless
//                                    frequency values
//     impact_parametter_file - input file containing impact parameter values
//     amplification_factor_real_file - output file containing the real parts of
//                                      the amplification factor values
//     amplification_factor_imag_file - output file containing the imaginary
//                                      parts of the amplification factor values
//     sum_threshold - integer value used as the threshold to stop summing at
//     precision - integer value used as the arithmetic preciusion in the
//                 amplification factor calculations
//
// Function takes in a dimensionless frequency and impact parameter file each
// containing a vector of numbers which it will then generate a pair of matrices
// of amplification factor values from these, outputting the real and imaginary
// components of this to two files. The infinite summation is approximated by
// using a finite threshold value and the arithmetic is done to a speciofed
// precision value
int main(int argc, char* argv[]) {
    // Read in all of the filenames
    std::string dimensionless_frequency_file = argv[1];
    std::string impact_parameter_file = argv[2];
    std::string amplification_factor_real_file = argv[3];
    std::string amplification_factor_imag_file = argv[4];

    // Read in the summation threshold and arithmetic precision values
    int sum_threshold = atoi(argv[5]);
    slong precision = atoi(argv[6]);

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

    // With the vectors generated, now perform the main loop of calculating
    // the amplification factor values and return these as a pair of matrices
    std::pair<std::vector<std::vector<double>>,
              std::vector<std::vector<double>>> amp_fac_matrices;
    amp_fac_matrices = AmplificationFactorMatrices(
        dimensionless_frequency, impact_parameter, sum_threshold, precision);

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
