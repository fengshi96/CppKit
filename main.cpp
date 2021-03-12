#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <typeinfo>
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

//    // test mv prod abstraction
//    std::vector<fcomplex> X = {{2,1},2};
//
//    std::vector<fcomplex> matini = {1,2,2,4};
//    Matrix<fcomplex> A(2,2,matini);
//
//    std::vector<fcomplex> R = A.dot(X);
//    for (const auto& x : R ) {
//        std::cout << x << " ";
//    }

        // test mm prod abstraction

        std::vector<dcomplex> matini1 = {1,2,2,4};
        Matrix<dcomplex> A(2,2,matini1);

        std::vector<dcomplex> matini2 = {1,1,2,0};
        Matrix<dcomplex> B(2,2,matini2);

        Matrix<dcomplex> R = A.dot(B, 'g');
        R.print();

//      // test vector vector multiplication
//      std::vector<fcomplex> X = {{1,1},2};
//      std::vector<fcomplex> Y = {4,5};
//      dcomplex R = dot(X,Y);
//      std::cout << R << " ";




}
