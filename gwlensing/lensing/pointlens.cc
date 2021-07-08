// Program definition files for calculating the amplification factor for the
// point mass lens and the main loop for performing the calculation based upon
// these
//
// Mick Wright 2021

#include "pointlens.h"

// Function computes the value of xm - the impact parameter divided by a length
// normalisation constant for the phase constant corresponding to the minimum
// time delay i.e. that of the image that travels the shortest path to the
// observer. Function is given by xm = (y + sqrt(y^2+4))/2
double StationaryPointMinimum(double impact_parameter) {
    double impact_parameter_sq = impact_parameter * impact_parameter;
    double xm = (impact_parameter + sqrt(impact_parameter_sq+4))/2;
    return xm;
}

// Function computes the value of phi - the phase constant used to obtain the
// minimum time delay induced by the lensing. Function is given by
// phi = (xm(y)-y)^2 / 2 - log(xm(y))
double MinTimeDelayPhaseConstant(double impact_parameter) {
    double xm_y = StationaryPointMinimum(impact_parameter);
    double x_minus_y = xm_y - impact_parameter;
    double phi = (x_minus_y*x_minus_y)/2 - log(xm_y);
    return phi;
}

// Function computes the amplification factor for an axially symmetric point
// mass lens for given values of dimensionless frequency and impact parameter
// with arithemtic precision given by precision
void AmplificationFactorCalculation(acb_t amplification_factor,
                                    double dimensionless_frequency,
                                    double impact_parameter,
                                    slong precision) {
    // Calculate the real part of the exponential component of the
    // amplification factor
    double exp_real = (M_PI * dimensionless_frequency)/4;

    // Calculate the imaginary part of the exponential component of the
    // amplification factor
    double exp_log_part = log(dimensionless_frequency/2);
    double exp_phi_part = 2 * MinTimeDelayPhaseConstant(impact_parameter);
    double exp_imag = (dimensionless_frequency/2) * (exp_log_part-exp_phi_part);

    // From the calculations construct the exponential component
    acb_t exponent;
    acb_t exponential_component;

    acb_init(exponent);
    acb_init(exponential_component);

    acb_set_d_d(exponent, exp_real, exp_imag);
    acb_exp(exponential_component, exponent, precision);

    // Construct the real and imaginary arguments of the gamma component
    // and then calculate the value of the gamma component
    double gamma_arg_real = 1;
    double gamma_arg_imag = -(dimensionless_frequency/2);

    acb_t gamma_argument;
    acb_t gamma_component;

    acb_init(gamma_argument);
    acb_init(gamma_component);

    acb_set_d_d(gamma_argument, gamma_arg_real, gamma_arg_imag);

    acb_gamma(gamma_component, gamma_argument, precision);

    // For those of non-trivial value calculate the real and imaginary parts
    // of the arguments for the confluent hypergeometric function component
    // of the amplification factor
    double hyper_arg_a_imag = dimensionless_frequency/2;

    double impact_parameter_sq = impact_parameter * impact_parameter;
    double hyper_arg_z_imag = hyper_arg_a_imag * impact_parameter_sq;

    // Construct the arguments to the confluent hypergeometric function
    acb_t hyper_arg_a;
    acb_t hyper_arg_b;
    acb_t hyper_arg_z;

    acb_init(hyper_arg_a);
    acb_init(hyper_arg_b);
    acb_init(hyper_arg_z);

    acb_set_d_d(hyper_arg_a, 0, hyper_arg_a_imag);
    acb_one(hyper_arg_b);
    acb_set_d_d(hyper_arg_z, 0, hyper_arg_z_imag);

    // Calculate the confluent hypergeometric function value
    acb_t hyper_component;
    acb_init(hyper_component);

    acb_hypgeom_1f1(hyper_component,
                    hyper_arg_a,
                    hyper_arg_b,
                    hyper_arg_z,
                    0,
                    precision);

    // Construct the final value of the amplification factor by first
    // calculating the value of the exponential * gamma components and then
    // finally multiplying by the hyper component value
    acb_t exponential_gamma;

    acb_init(exponential_gamma);

    acb_mul(exponential_gamma,
            exponential_component,
            gamma_component,
            precision);
    acb_mul(amplification_factor,
            exponential_gamma,
            hyper_component,
            precision);

    // Clean up of declared acb_ts for memory management
    acb_clear(exponent);
    acb_clear(exponential_component);
    acb_clear(gamma_argument);
    acb_clear(gamma_component);
    acb_clear(hyper_arg_a);
    acb_clear(hyper_arg_b);
    acb_clear(hyper_arg_z);
    acb_clear(hyper_component);
    acb_clear(exponential_gamma);
}

// Function constructs two matrices containing the real and imaginary parts of
// the value of the amplification factor function based upon two vectors
// containing values of dimensionless frequency and impact parameters. It
// returns these inside of a pair object
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
    AmplificationFactorMatrices(std::vector<double> dimensionless_frequency,
                                std::vector<double> impact_parameter,
                                slong precision) {
    // Get the size of the vectors
    int impact_parameter_size = impact_parameter.size();
    int dimensionless_frequency_size = dimensionless_frequency.size();

    // Set up the vectors for the amplification factor real and imaginary parts
    std::vector<std::vector<double>> amplification_factor_real(
        impact_parameter_size,
        std::vector<double> (dimensionless_frequency_size));
    std::vector<std::vector<double>> amplification_factor_imag(
        impact_parameter_size,
        std::vector<double> (dimensionless_frequency_size));

    // Create the loop that goes through the structures to generate the
    // amplification factor values for the dimensionless frequency and impact
    // parameter combinations. This loop is performed in parallel and using
    // the dynamic method - this is because the amount of time that will be
    // spent on each calculation is not equal
    #pragma omp parallel for collapse(2) schedule(dynamic)

    for (int i=0; i < impact_parameter_size; i++) {
        for (int j=0; j < dimensionless_frequency_size; j++) {
            acb_t amplification_factor;
            acb_init(amplification_factor);
            AmplificationFactorCalculation(amplification_factor,
                                           dimensionless_frequency[j],
                                           impact_parameter[i],
                                           precision);

            amplification_factor_real[i][j] = arf_get_d(
                arb_midref(acb_realref(amplification_factor)), ARF_RND_NEAR);
            amplification_factor_imag[i][j] = arf_get_d(
                arb_midref(acb_imagref(amplification_factor)), ARF_RND_NEAR);
        }
    }

    // Wrap the amplification factor real and imaginary parts into a single
    // pair object and return it
    std::pair<std::vector<std::vector<double>>,
              std::vector<std::vector<double>>>
        amplification_factor_matrices(amplification_factor_real,
                                      amplification_factor_imag);

    return amplification_factor_matrices;
}

// Main Function - this takes in five arguments:
//    dimensionless_frequency_file - input file contianing dimensionless
//                                   frequency values
//    impact_parameter_file - input file containing impact parameter values
//    amplification_factor_real_file - output file containing the real parts of
//                                     the amplification factor values
//    amplification_factor_imag_file - output file containing the imaginary
//                                     parts of the amplification factor values
//    precision - integer value to be used as the arithmetic precision in the
//                amplification factor calculations
//
// Function takes in a dimensionless frequency and impact parmaeter file each
// containing a vector of numbers which it will then generate a matrix of
// amplification factor values of from these outputting these to the two
// amplification factor files
int main(int argc, char* argv[]) {
    // Read in all of the filenames
    std::string dimensionless_frequency_file = argv[1];
    std::string impact_parameter_file = argv[2];
    std::string amplification_factor_real_file = argv[3];
    std::string amplification_factor_imag_file = argv[4];

    // Read in the precision value
    slong precision = atoi(argv[5]);

    // For each of the dimensionless frequency and impact parmaeter files
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
        dimensionless_frequency, impact_parameter, precision);

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
