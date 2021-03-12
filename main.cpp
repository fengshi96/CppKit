#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include "src/Matrix.h"

typedef std::complex<double> dcomplex;
typedef std::complex<float> fcomplex;

int main() {
    std::vector<fcomplex> data = {{2,1}, 1.2, 1.4, 2};
    Matrix<fcomplex> M(2,2, data);  // the matrix A to be diagonalized

    // assert(M.IsHermitian());
    // M.fillRand();
    M.print();
    M.ajoint();
    M.print();

    std::vector<fcomplex> evals;
    std::vector<fcomplex> evecs;

    diag(M, evals, evecs, 'V');

    for(const auto& x : evals) {
        std::cout << x << " ";
    }
    std::cout << std:: endl << std:: endl << std:: endl;


    // test transpose
    std::vector<fcomplex> vec = {{1,2}, {2,5}, {1,9}, {11,12}, {223, 11}, {0,99}};
    Matrix<fcomplex> N(2, 3, vec);
    N.print();
    std::cout << std::endl;

    N.ajoint();
    N.print();




}
