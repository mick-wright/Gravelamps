// Function definition for the lensing data generation for the isolated point
// mass lens model
//
// Mick Wright 2022

#include "src/point.h"

volatile sig_atomic_t checkpoint_flag = 0;

// Checkpoint signal handler
void checkpoint_sig_handler(int signum) {
    checkpoint_flag = 85;
}

// Function takes in values of dimensionless frequency and source position and
// calculates the amplification factor for the model using either wave or
// geometric optics depending upon whether the dimensionless frequency is above
// the specified switch. It returns in either case a pair of doubles containing
// the real and imaginary components
//
// Input:
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//      slong precision : arithmetic precision to use in arbitrary precision
//                        calculations
//      slong geo_switch : value of dimensionless frequency at which to switch
//                         from full wave optics to the geometric optics
//                         approximation
//
// Output:
//      double* amplification_factor : Array of two doubles containing the real
//                                     and imaginary components of the
//                                     amplification factor
double* AmplificationFactor(double dimensionless_frequency,
                            double source_position,
                            slong precision,
                            slong geo_switch) {
    double* amplification_factor;
    if (dimensionless_frequency > geo_switch) {
        amplification_factor =\
            PyAmplificationFactorGeometric(dimensionless_frequency,
                                           source_position);
    } else {
        acb_t amplification;
        acb_init(amplification);
        AmplificationFactorWave(amplification,
                                dimensionless_frequency,
                                source_position,
                                precision);
        amplification_factor = new double[2];
        amplification_factor[0] =
            arf_get_d(arb_midref(acb_realref(amplification)), ARF_RND_NEAR);
        amplification_factor[1] =
            arf_get_d(arb_midref(acb_imagref(amplification)), ARF_RND_NEAR);
    }
    return amplification_factor;
}

// Function takes in a pair of files containing a grid of values of
// dimensionless frequency and source position and generates a pair of files
// containing the corresponding grid of the real and imaginary part sof the
// amplification factor. If these files exist already it will read them in and
// use them to avoid repeating calcualtions, allowing for the process to be
// interrupted.
//
// Input:
//      char* dimensionless_frequency_file : path of file containing
//                                           dimensionless frequency values
//                                           over which to calculate the
//                                           amplification factor
//      char* source_position_file : path of file containing source position
//                                   values over which to caclulate the
//                                   amplification factor
//      char* amplification_factor_real_file : path of file over which to
//                                             output the real component of the
//                                             calculated amplification factor
//                                             values
//      char* amplification_factor_imag_file : path of file over which to
//                                             output the imaginary component
//                                             of the calculated amplification
//                                             factor values
//      slong precision : arithmetic precision to use in arbitrary precision
//                        calculations
//      slong geo_swtich : value of dimensionless frequency at which to switch
//                         from full wave optics to the geometric optics
//                         approximation
int GenerateLensData(char* dimensionless_frequency_file,
                     char* source_position_file,
                     char* amplification_factor_real_file,
                     char* amplification_factor_imag_file,
                     int precision_int,
                     int geo_switch_int) {
    // Set up Signal handlers
    std::signal(SIGALRM, checkpoint_sig_handler);
    std::signal(SIGINT, checkpoint_sig_handler);

    slong precision = precision_int;
    slong geo_switch = geo_switch_int;

    // Get vectors of the dimensionless frequency and source position values
    // from the input files
    std::vector<double> dimensionless_frequency =\
        GetVector(dimensionless_frequency_file);
    std::vector<double> source_position = GetVector(source_position_file);

    // Store the lengths of these vectors for matrix generation
    int dimensionless_frequency_size = dimensionless_frequency.size();
    int source_position_size = source_position.size();

    // Generate the empty amplification factor matrices and populate with any
    // already extant values from the output file locations
    std::vector<std::vector<double>> amplification_factor_real =\
        GetMatrix(amplification_factor_real_file,
                  source_position_size,
                  dimensionless_frequency_size);

    std::vector<std::vector<double>> amplification_factor_imag =\
        GetMatrix(amplification_factor_imag_file,
                  source_position_size,
                  dimensionless_frequency_size);

    // Open outstreams for the amplification factor matrices
    std::ofstream real_outstream(amplification_factor_real_file);
    std::ofstream imag_outstream(amplification_factor_imag_file);

    // Loop through the created matrices. If the value is zero at a given point
    // it is assumed to have not been read from the file and is calculated
    // using either wave or geometric optics as determined by the geo_switch
    // value. This is done with parallel threading.
    for (int i=0; i < source_position_size; i++) {
        if (checkpoint_flag == 85) {
             std::cout << "Received checkpoint instruction" << std::endl;
             return checkpoint_flag;
        }
        #pragma omp parallel for ordered schedule(dynamic) shared(checkpoint_flag)
        for (int j=0; j < dimensionless_frequency_size; j++) {
            if (checkpoint_flag == 85) {
                continue;
            }
            if (amplification_factor_real[i][j] != 0.0
                && amplification_factor_imag[i][j] != 0.0) {
                    {}
            } else {
                double* amplification_factor =\
                    AmplificationFactor(dimensionless_frequency[j],
                                        source_position[i],
                                        precision,
                                        geo_switch);

                amplification_factor_real[i][j] = amplification_factor[0];
                amplification_factor_imag[i][j] = amplification_factor[1];
            }
            #pragma omp ordered
            {
            real_outstream << amplification_factor_real[i][j] << " ";
            imag_outstream << amplification_factor_imag[i][j] << " ";

            std::cout << "Completed dimensionless frequency " <<
                std::setfill('0') << std::setw(5) << j+1 << " of " <<
                dimensionless_frequency_size << "\r";

            real_outstream.flush();
            imag_outstream.flush();
            }
        }
        real_outstream << std::endl;
        imag_outstream << std::endl;

        std::cout << std::endl << "Completed source position value " << i+1 <<
            " of " << source_position_size << std::endl;
    }

    if (checkpoint_flag == 85) {
        return checkpoint_flag;
    } else {
        return 0;
    }
}
