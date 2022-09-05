// Header file containing generic utility functions that can be used in any of
// the models.
//
// Mick Wright 2022

#ifndef GRAVELAMPS_MODEL_SRC_UTILS_H_
#define GRAVELAMPS_MODEL_SRC_UTILS_H_

#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

// Function takes a filename and reads it into a vector of values.
std::vector<double> GetVector(char* filename);

// Function takes a filename as well as matrix size to generate a matrix of the
// specified size and populate it with any values that are already present in
// the specified file name.
std::vector<std::vector<double>> GetMatrix(char* filename,
                                           int row_num,
                                           int col_num);

// Function destroys objects placed within, used for deallocating the memory
// of resultant objects sent to the python frontend
extern "C" {
    void destroyObj(double* object);
}

#endif  // GRAVELAMPS_MODEL_SRC_UTILS_H_
