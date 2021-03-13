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
//    std::vector<dcomplex> X = {{2,1},2};
//
//    std::vector<dcomplex> matini = {1,2,2,4};
//    Matrix<dcomplex> A(2,2,matini);
//
//    std::vector<dcomplex> R = A.dot(X);
//    dcomplex a = 1.0;
//    vscal(a, R);
//    for (const auto& x : R ) {
//        std::cout << x << " ";
//    }

//        // test mm prod abstraction
//
//        std::vector<double> matini1 = {1,2,2,4};
//        Matrix<double> A(2,2,matini1);
//
//        std::vector<double> matini2 = {1,1,2,0};
//        Matrix<double> B(2,2,matini2);
//
//        Matrix<double> R = A.dot(B, 'g');
//        double a = 3.0;
//        vscal(a, R);
//        R.print();

//      // test vector vector multiplication
//      std::vector<fcomplex> X = {{1,1},2};
//      std::vector<fcomplex> Y = {4,5};
//      dcomplex R = dot(X,Y);
//      std::cout << R << " ";

//        // test m +- m
//        std::vector<double> matini1 = {1,2,2,4};
//        Matrix<double> A(2,2,matini1);
//        A.print();
//        std::cout << std::endl;
//
//        std::vector<double> matini2 = {1,3,2,22};
//        Matrix<double> B(2,2,matini2);
//        B.print();
//        std::cout << std::endl;
//
//        A -= B;
//        A.print();

//        // test a*matrix
//        std::vector<double> matini1 = {1,2,2,4};
//        Matrix<double> A(2,2,matini1);
//        A.print();
//        std::cout << std::endl;
//        A *= 4;
//        A.print();
        //test norm
        std::vector<dcomplex> v = {{1,1},1};
        std::cout << norm(v) << std::endl;






}
