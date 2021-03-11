#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include "src/Matrix.h"

typedef std::complex<double> dcomplex;
typedef std::complex<float> fcomplex;

int main() {
    std::vector<fcomplex> data = {{8,1}, 1.2, 1.2, 2};
    Matrix<fcomplex> M(2,2, data);  // the matrix A to be diagonalized

    // assert(M.IsHermitian());
    // M.fillRand();
    M.print();

    std::vector<float> evals;

    diag(M, evals, 'V');

    for(const auto& x : evals) {
        std::cout << x << " ";
    }
    M.print();



}
