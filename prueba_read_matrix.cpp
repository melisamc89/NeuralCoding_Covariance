#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <string>
#include <iostream>

#include "nrutilc.h"
#include "iomatrix.h" // aca esta definida readMatrixFromFile mas generica

using namespace std;

// sobrecargada para no tener que pasar el tamaño. si puede poner en un .h aparte
template <class T>
bool readMatrixFromFile(const char* filename, NRMat<T>& matrix)
{
    return readMatrixFromFile(filename, matrix, matrix.nrows(), matrix.ncols());
}


int main()
{
	const unsigned int n=30001, m=4;
    string name("gs_15.dat"); 
    
	NRMat<double> mat(n, m);  // cuidado que la libreria no cheque n, m > 0
    if(!readMatrixFromFile(name.c_str(), mat, n, m)) {
        cout << "Error al leer " << name << endl;
        return -1;
    }
    
	for(unsigned int i = 0; i < n; ++i){
		for(unsigned int j = 0; j < m; ++j){
			printf("%lf\t", mat[i][j]);
		}
		printf("\n");
	}
    
    return 0;	
}
