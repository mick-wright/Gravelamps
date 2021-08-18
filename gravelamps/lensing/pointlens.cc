// Program definition files for calculating the amplification factor for the
// point mass lens and the main loop for performing the calculation based upon
// these
//
// Mick Wright 2021

#include "pointlens.h"

// Function computes the value of xm - the source position divided by a length
// normalisation constant for the phase constant corresponding to the minimum
// time delay i.e. that of the image that travels the shortest path to the
// observer. Function is given by xm = (y + sqrt(y^2+4))/2
double StationaryPointMinimum(double source_position) {
    double source_position_sq = source_position * source_position;
    double xm = (source_position + sqrt(source_position_sq+4))/2;
    return xm;
}

// Function computes the value of phi - the phase constant used to obtain the
// minimum time delay induced by the lensing. Function is given by
// phi = (xm(y)-y)^2 / 2 - log(xm(y))
double MinTimeDelayPhaseConstant(double source_position) {
    double xm_y = StationaryPointMinimum(source_position);
    double x_minus_y = xm_y - source_position;
    double phi = (x_minus_y*x_minus_y)/2 - log(xm_y);
    return phi;
}

// Function computes the amplification factor for an axially symmetric point
// mass lens using full wave optics for given values of dimensionless frequency
// and source position with arithmetic precision given by precision
void AmplificationFactorCalculation(acb_t amplification_factor,
                                    double dimensionless_frequency,
                                    double source_position,
                                    slong precision) {
    // Calculate the real part of the exponential component of the
    // amplification factor
    double exp_real = (M_PI * dimensionless_frequency)/4;

    // Calculate the imaginary part of the exponential component of the
    // amplification factor
    double exp_log_part = log(dimensionless_frequency/2);
    double exp_phi_part = 2 * MinTimeDelayPhaseConstant(source_position);
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

    double source_position_sq = source_position * source_position;
    double hyper_arg_z_imag = hyper_arg_a_imag * source_position_sq;

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

// Function computes the mangification for the geometric optics approximation
// with the plus and minus images given by the state of mode. This is given by
// 1/2 +- (y^2 + 2)/(2y * sqrt(y^2 + 4)). Mode determines plus or minus
double Magnification(double source_position, int mode) {
    double magnification = 1./2.;
    double second_term_numerator = source_position*source_position + 2;
    double second_term_denominator = 2 * source_position * sqrt(
        source_position * source_position + 4);
    double second_term = second_term_numerator/second_term_denominator;

    if (mode == 1) {
        magnification += second_term;
    } else {
        magnification -= second_term;
    }

    return magnification;
}

// Function computes the time delay for the geometric optics approximation.
// This is given by:
// y * (sqrt(y^2 + 4)/2) + ln((sqrt(y^2+4)+y)/(sqrt(y^2+4)-y))
double TimeDelay(double source_position) {
    double first_term = source_position;
    double first_term_numerator = sqrt(source_position*source_position + 4);
    first_term *= first_term_numerator/2.;

    double second_term_numerator = first_term_numerator + source_position;
    double second_term_denominator = first_term_numerator - source_position;
    double second_term = log(second_term_numerator/second_term_denominator);

    double time_delay = first_term + second_term;

    return time_delay;
}

// Function computes the amplification factor for an axially symmetric point
// mass lens using the geometric optics approximation for given values of
// dimensionless frequency and source position
std::complex<double> AmplificationFactorGeometric(
    double dimensionless_frequency, double source_position) {
    // Using the complex_literals i operator
    using std::literals::complex_literals::operator""i;

    // First compute the magnifications
    double magnification_plus = Magnification(source_position, 1);
    double magnification_minus = Magnification(source_position, 0);

    // Compute the Time Delay
    double time_delay = TimeDelay(source_position);

    // Calculate the exponential term
    double exponent_imag  = dimensionless_frequency * time_delay;
    std::complex<double> exp_term = 1i * exponent_imag;
    exp_term = exp(exp_term);

    // Construct the final amplification factor
    double first_term = sqrt(abs(magnification_plus));
    std::complex<double> second_term =
        sqrt(abs(magnification_minus)) * exp_term;
    std::complex<double> geometric_factor = first_term - 1i * second_term;

    return geometric_factor;
}

// Function constructs two matrices containing the real and imaginary parts of
// the value of the amplification factor function based upon two vectors
// containing values of dimensionless frequency and source position. It
// returns these inside of a pair object
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
    AmplificationFactorMatrices(std::vector<double> dimensionless_frequency,
                                std::vector<double> source_position,
                                slong precision,
                                int approx_switch) {
    // Get the size of the vectors
    int source_position_size = source_position.size();
    int dimensionless_frequency_size = dimensionless_frequency.size();

    // Set up the vectors for the amplification factor real and imaginary parts
    std::vector<std::vector<double>> amplification_factor_real(
        source_position_size,
        std::vector<double> (dimensionless_frequency_size));
    std::vector<std::vector<double>> amplification_factor_imag(
        source_position_size,
        std::vector<double> (dimensionless_frequency_size));

    // Create the loop that goes through and calculates the amplification
    // factor values over which the amplification factor values will be
    // calculated through the full wave optics calculations. The schedule
    // here is dynamic due to the fact that the calculation takes differing
    // amounts of time at differing points
    #pragma omp parallel for collapse(2) schedule(dynamic)

    for (int i=0; i < source_position_size; i++) {
        for (int j=0; j < approx_switch; j++) {
            acb_t amplification_factor;
            acb_init(amplification_factor);
            AmplificationFactorCalculation(amplification_factor,
                                           dimensionless_frequency[j],
                                           source_position[i],
                                           precision);

            amplification_factor_real[i][j] = arf_get_d(
                arb_midref(acb_realref(amplification_factor)), ARF_RND_NEAR);
            amplification_factor_imag[i][j] = arf_get_d(
                arb_midref(acb_imagref(amplification_factor)), ARF_RND_NEAR);
        }
    }

    // Create the loop that goes through and calculates the amplification
    // factor values over which they will be calculated using the geometric
    // optics approximations
    if (approx_switch < dimensionless_frequency_size) {
        #pragma omp parallel for collapse(2) schedule(dynamic)
        for (int i=0; i < source_position_size; i++) {
            for (int j=approx_switch; j < dimensionless_frequency_size; j++) {
                std::complex<double> geometric_factor;
                geometric_factor =
                    AmplificationFactorGeometric(dimensionless_frequency[j],
                                                 source_position[i]);

                amplification_factor_real[i][j] = std::real(geometric_factor);
                amplification_factor_imag[i][j] = std::imag(geometric_factor);
            }
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

// Main Function - this takes in six arguments:
//    dimensionless_frequency_file - input file contianing dimensionless
//                                   frequency values
//    source_position_file - input file containing source position values
//    amplification_factor_real_file - output file containing the real parts of
//                                     the amplification factor values
//    amplification_factor_imag_file - output file containing the imaginary
//                                     parts of the amplification factor values
//    precision - integer value to be used as the arithmetic precision in the
//                amplification factor calculations
//    approx_switch - integer value which is the location in the dimensionless
//                    frequency array where the geometric optics location
//                    should take over
//
// Function takes in a dimensionless frequency and impact parmaeter file each
// containing a vector of numbers which it will then generate a matrix of
// amplification factor values of from these outputting these to the two
// amplification factor files
int main(int argc, char* argv[]) {
    // Read in all of the filenames
    std::string dimensionless_frequency_file = argv[1];
    std::string source_position_file = argv[2];
    std::string amplification_factor_real_file = argv[3];
    std::string amplification_factor_imag_file = argv[4];

    // Read in the precision value
    slong precision = atoi(argv[5]);

    // Read in the position where the geometric optics approximation
    // will take over
    slong approx_switch = atoi(argv[6]);

    // For each of the dimensionless frequency and impact parmaeter files
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
        dimensionless_frequency, source_position, precision, approx_switch);

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
