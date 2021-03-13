#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include "src/Vector.h"
#include "src/Matrix.h"

typedef std::complex<double> dcomplex;
typedef std::complex<float> fcomplex;

int main() {

//    // test vector vector product
//    std::vector<double> X = {2,2};
//    std:: vector<double> Y = {4,5};
//    int R = vxv(X, Y);
//    std::cout << R << std::endl;
//
//    std::cout << (typeid(dcomplex).name() == typeid(std::complex<double>).name()) <<std::endl;
//    std::cout << typeid(std::complex<double>).name() <<std::endl;

    // test mv prod abstraction
    std::vector<dcomplex> X = {1,2};
    std::vector<dcomplex> matini = {1,2,2,4};
    Matrix<dcomplex> A(2,2,matini);

    std::vector<dcomplex> R0 = prod(A, X);
    std::vector<dcomplex> R = A.prod(X);
    R0 *= 1;
    for (const auto& x : R ) {
        std::cout << x << " ";
    }






}
