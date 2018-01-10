
#ifndef IOMATRIX_H
#define IOMATRIX_H

#include <iostream>
#include <fstream>

// devuelve true si en caso de exito, false caso contrario
template <class Mat>
bool readMatrixFromFile(const char* filename, Mat& matrix, unsigned int nrows, unsigned int ncols)
{    
    std::ifstream input(filename,  std::ifstream::in); // abro el archivo
    if (!input.is_open()) return false; // error al abrir el archivo
    
    for (unsigned int i = 0; i < nrows; ++i){
        for (unsigned int j = 0; j < ncols; ++j){
            input >> matrix[i][j]; // la matrix debe tener el tamaño adecuado
            if (!input.good()) return false; // error al leer el archivo
        }
    }
    
    input.close(); // cierro el archivo
    return true;   // no ocurrio ningun error
}

#endif
