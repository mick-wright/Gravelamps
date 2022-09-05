// Generic utility function definitions for the assistance in the calculation
// of the amplification factor by any of the models
//
// Mick Wright 2022

#include "src/utils.h"

// Function takes a filename and reads it into a vector of values
//
// Input:
//      char* filename : filename to read vector from
//
// Output:
//      std::vector<double> value_vector : vector containing values from file
std::vector<double> GetVector(char* filename) {
    std::ifstream file_stream(filename);
    std::istream_iterator<double> file_start(file_stream), file_end;
    std::vector<double> value_vector(file_start, file_end);

    return value_vector;
}

// Function takes a filename as well as matrix size to generate a matrix of the
// specified size and populate it with any values that are already present
// in the specified file name.
//
// Input:
//      char* filename : filename to read any completed values from
//      int row_num : number of rows for the full matrix
//      int col_num : number of columns for the full matrix
//
// Output:
//      std::vector<std::vector<double>> matrix : matrix of specified size with
//                                                any values from the file
//                                                placed within
std::vector<std::vector<double>> GetMatrix(char* filename,
                                           int row_num,
                                           int col_num) {
    // Create temporary matrix to fill with any values from the file
    std::vector<std::vector<double>> tmp_matrix;

    // Open the file for reading
    std::ifstream file_stream(filename);

    // Read in any values from the file if it's open
    if (file_stream.is_open()) {
        std::string line;
        double value;

        while (std::getline(file_stream, line)) {
            std::istringstream iss(line);
            std::vector<double> tmp;
            while (iss >> value) {
                tmp.push_back(value);
            }
            tmp_matrix.push_back(tmp);
        }
    }

    // Close the file
    file_stream.close();

    // Create the full specified matrix
    std::vector<std::vector<double>> matrix(row_num,
                                            std::vector<double>(col_num));

    // Populate the matrix with any values from the file
    int row_count = 0;
    for (auto i : tmp_matrix) {
        int column_count = 0;
        for (auto j : i) {
            matrix[row_count][column_count] = j;
            column_count++;
        }
        row_count++;
    }

    return matrix;
}

// Function destroys object placed within, used for deallocating the memory
// of resultant objects sent to the python frontend.
//
// Input:
//      double* object : object to be destroyed
void destroyObj(double* object) {
    delete object;
    object = NULL;
}

