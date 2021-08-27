// Program definition files for calculating the amplification factor for the
// Singular Isothermal Sphere (SIS) lens and the main loop for performing the
// calculation based upon these
//
// Mick Wright

#include "sis.h"

// Function computes the amplification factor for an axially symmetric singular
// isothermal sphere (SIS) style lensing using full wave optics for given
// values of dimensionless frequency and source position. It does this using a
// summation method with the infinite sum approximated up to a given threshold
// and with a given arithmetic precision
void AmplificationFactorCalculation(acb_t amplification_factor,
                                    double dimensionless_frequency,
                                    double source_position,
                                    slong summation_upper_limit,
                                    slong precision) {
    // For the given value of impact parameter, calculate the phi, the phase
    // constant used to obtain a minimum time delay of zero
    double phi = source_position + 1./2.;

    // Calculate the prefactor term for the summation. This is given by
    // exp(i(w/2 * (y^2 + 2phi(y))))
    acb_t prefactor;
    acb_init(prefactor);

    double prefactor_dim_freq_term = dimensionless_frequency/2.;
    double prefactor_sour_pos_term = source_position * source_position;
    double prefactor_phi_term = 2. * phi;

    double prefactor_exponent_imag = prefactor_dim_freq_term
                                     * (prefactor_sour_pos_term
                                     + prefactor_phi_term);

    acb_set_d_d(prefactor, 0, prefactor_exponent_imag);
    acb_exp(prefactor, prefactor, precision);

    // Initialise the summation term
    acb_t summation_term;
    acb_init(summation_term);
    acb_zero(summation_term);

    // Loop through the summation from 0 up to the upper summation limit
    for (int n=0; n <= summation_upper_limit; n++) {
        // Calculate the gamma term, gamma(1+n/2)/n!
        double gamma_term_val = 1 + n/2.;

        acb_t gamma_term;
        acb_init(gamma_term);
        acb_set_d(gamma_term, gamma_term_val);
        acb_gamma(gamma_term, gamma_term, precision);

        arb_t factorial;
        arb_init(factorial);
        arb_fac_ui(factorial, n, precision);

        acb_div_arb(gamma_term, gamma_term, factorial, precision);

        // Calculate the power term (2w e^i3pi/2)^n/2
        double power_term_prefactor_val = 2 * dimensionless_frequency;
        double power_term_exponent_val = 3 * M_PI/2.;
        double power_term_power_val = n/2.;

        acb_t power_term;
        acb_t power_term_exponent;
        arb_t power_term_power;

        acb_init(power_term);
        acb_init(power_term_exponent);
        arb_init(power_term_power);

        acb_set_d(power_term, power_term_prefactor_val);
        acb_set_d_d(power_term_exponent, 0, power_term_exponent_val);
        arb_set_d(power_term_power, power_term_power_val);

        acb_exp(power_term_exponent, power_term_exponent, precision);
        acb_mul(power_term, power_term, power_term_exponent, precision);
        acb_pow_arb(power_term, power_term, power_term_power, precision);

        // Calculate the Hypergeometric function term:
        // 1f1(1 + n/2; 1; -i/2 wy^2)
        double hyper_arg_a_real = 1 + n/2.;
        double hyper_arg_z_imag = (-1./2.) * dimensionless_frequency
                                  * source_position * source_position;
	
        acb_t hyper_arg_a;
        acb_t hyper_arg_b;
        acb_t hyper_arg_z;

        acb_init(hyper_arg_a);
        acb_init(hyper_arg_b);
        acb_init(hyper_arg_z);

        acb_set_d(hyper_arg_a, hyper_arg_a_real);
        acb_one(hyper_arg_b);
        acb_set_d_d(hyper_arg_z, 0, hyper_arg_z_imag);

        acb_t hyper_term;
        acb_init(hyper_term);
        acb_hypgeom_1f1(
           hyper_term, hyper_arg_a, hyper_arg_b, hyper_arg_z, 0, precision);

        // Construct the value to add to the summation
        acb_t summation_temp;
        acb_init(summation_temp);
        acb_mul(summation_temp, gamma_term, power_term, precision);
        acb_mul(summation_temp, summation_temp, hyper_term, precision);

	acb_add(summation_term, summation_term, summation_temp, precision);

	// Memory Management - clear the declared acbs
	acb_clear(gamma_term);
	arb_clear(factorial);
	acb_clear(power_term);
	acb_clear(power_term_exponent);
	arb_clear(power_term_power);
	acb_clear(hyper_arg_a);
	acb_clear(hyper_arg_b);
	acb_clear(hyper_arg_z);
	acb_clear(hyper_term);
	acb_clear(summation_temp);
    }

    // Get the final result by multiplying the prefactor and the summation term
    acb_mul(amplification_factor, prefactor, summation_term, precision);

    // Memory Management - clear the remaining acbs
    acb_clear(prefactor);
    acb_clear(summation_term);
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
                                slong summation_upper_limit,
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
    for (int i=0; i < source_position_size; i++) {
        #pragma omp parallel for schedule(dynamic)
        for (int j=0; j < approx_switch; j++) {
            acb_t amplification_factor;
            acb_init(amplification_factor);
            AmplificationFactorCalculation(amplification_factor,
                                           dimensionless_frequency[j],
                                           source_position[i],
                                           summation_upper_limit,
                                           precision);
            amp_fac_real[i][j] = arf_get_d(
                arb_midref(acb_realref(amplification_factor)), ARF_RND_NEAR);
            amp_fac_imag[i][j] = arf_get_d(
                arb_midref(acb_imagref(amplification_factor)), ARF_RND_NEAR);
        }
        std::cout << "Completed column " << i+1 << " of "
                  << source_position_size << std::endl;
    }

    // Create the loop that goes through and calculates the amplification
    // factor values over which they will be calculated using the geometric
    // optics approximation
    if (approx_switch < dimensionless_frequency_size) {
        for (int i=0; i < source_position_size; i++) {
            #pragma omp parallel for schedule(dynamic)
            for (int j=approx_switch; j < dimensionless_frequency_size; j++) {
                std::complex<double> geometric_factor;
                geometric_factor =
                    AmplificationFactorGeometric(dimensionless_frequency[j],
                                                 source_position[i]);

                amp_fac_real[i][j] = std::real(geometric_factor);
                amp_fac_imag[i][j] = std::imag(geometric_factor);
            }
            std::cout << "Completed column " << i+1 << " of "
                      << source_position_size << std::endl;
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
//     summation_upper_limit - value to calculate the summation up to
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
    slong summation_upper_limit = atoi(argv[5]);
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
        summation_upper_limit, precision, approx_switch);

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
