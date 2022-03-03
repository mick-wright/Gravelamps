// Program definition files for generating interpolation data for the point
// mass lens model in both the wave and geometric optics regimes for user
// given parameters
//
// Mick Wright 2021

#include "src/pointlens.h"

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
    // Read in all of the filenames from the arguments
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

    //Close the streams
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
                amp_fac_realout << amplification_factor_real[i][j] << "\t";
                amp_fac_imagout << amplification_factor_imag[i][j] << "\t";
            }
        }
        #pragma omp critical
        amp_fac_realout << std::endl;
        amp_fac_imagout << std::endl;

        std::cout << "Completed source position value " << i+1 << " of " <<
           source_position_size << std::endl;
    }

    return 0;
}
