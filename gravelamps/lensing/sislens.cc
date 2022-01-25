// Program definition files for generating interpolation data for the Singular
// Isothermal Sphere (SIS) lens model in both wave and geometric optics for
// user given data 
//
// Mick Wright

#include "src/sis.h"

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

    // For each of the dimensionless frequency and impact parmaeter files
    // open the file for reading, generate an iterator object which will get
    // all of the enclosed values and put them inside of a vector
    std::ifstream dim_freq_fstream(dimensionless_frequency_file);
    std::istream_iterator<double> dim_freq_start(dim_freq_fstream),
                                  dim_freq_end;
    std::vector<double> dimensionless_frequency(dim_freq_start, dim_freq_end);

    std::ifstream sour_pos_fstream(source_position_file);
    std::istream_iterator<double> sour_pos_start(sour_pos_fstream),
                                  sour_pos_end;
    std::vector<double> source_position(sour_pos_start, sour_pos_end);

    // Store the lengths of these vectors for matrix generation
    int source_position_size = source_position.size();
    int dimensionless_frequency_size = dimensionless_frequency.size();

    // Set up the matrices for the amplification factor real and imaginary
    // parts
    std::vector<std::vector<double>> amplification_factor_real;
    std::vector<std::vector<double>> amplification_factor_imag;

    // If the amplification factor files exist already, read them in to avoid
    // repeating calculation
    std::ifstream amp_fac_real(amplification_factor_real_file);
    std::ifstream amp_fac_imag(amplification_factor_imag_file);

    if (amp_fac_real.is_open() && amp_fac_imag.is_open()) {
        // Load into the matrix from the file
        std::string line;
        double value;

        while (std::getline(amp_fac_real, line)) {
            std::istringstream iss(line);
            std::vector<double> tmp;
            while (iss >> value) {
               tmp.push_back(value);
            }
            amplification_factor_real.push_back(tmp);
        }

        while (std::getline(amp_fac_imag, line)) {
            std::istringstream iss(line);
            std::vector<double> tmp;
            while (iss >> value) {
               tmp.push_back(value);
            }
            amplification_factor_imag.push_back(tmp);
        }
    }

    // Now resize the matrix to the correct size
    amplification_factor_real.resize(
        source_position_size,
        std::vector<double>(dimensionless_frequency_size));
    amplification_factor_imag.resize(
        source_position_size,
        std::vector<double>(dimensionless_frequency_size));

    // Close the istreams
    amp_fac_real.close();
    amp_fac_imag.close();

    // Get outstreams for the amplification factor matrices
    std::ofstream amp_fac_realout(amplification_factor_real_file);
    std::ofstream amp_fac_imagout(amplification_factor_imag_file);

    // Loop through the dimensionless frequency and source position values
    // if the amplification factor value is zero (i.e. not read in from file
    // already) calculate the value in either wave or geometric optics as
    // determined by the user's switchover value. This is done with parallel
    // threading
    for (int i=0; i < source_position_size; i++) {
        #pragma omp parallel for ordered schedule(dynamic)
        for (int j=0; j < dimensionless_frequency_size; j++) {
            if (amplification_factor_real[i][j] != 0.0
                && amplification_factor_imag[i][j] != 0.0) {
                {}
            } else {
                if (j >= approx_switch) {
                    std::complex<double> geometric_factor;
                    geometric_factor =
                        AmplificationFactorGeometric(
                            dimensionless_frequency[j],
                            source_position[i]);

                    amplification_factor_real[i][j] =
                        std::real(geometric_factor);
                    amplification_factor_imag[i][j] =
                        std::imag(geometric_factor);
                } else {
                    acb_t amplification_factor;
                    acb_init(amplification_factor);

                    AmplificationFactorCalculation(amplification_factor,
                                                   dimensionless_frequency[j],
                                                   source_position[i],
                                                   summation_upper_limit,
                                                   precision);

                    amplification_factor_real[i][j] = arf_get_d(
                        arb_midref(acb_realref(amplification_factor)),
                        ARF_RND_NEAR);
                    amplification_factor_imag[i][j] = arf_get_d(
                        arb_midref(acb_imagref(amplification_factor)),
                        ARF_RND_NEAR);
                }
            }
            #pragma omp ordered
            {
                amp_fac_realout << amplification_factor_real[i][j] << " ";
                amp_fac_imagout << amplification_factor_imag[i][j] << " ";
            }
        }
        amp_fac_realout << "\n";
        amp_fac_imagout << "\n";

        std::cout << "Completed source position value " << i+1 << " of " <<
           source_position_size << std::endl;
    }

    return 0;
}
