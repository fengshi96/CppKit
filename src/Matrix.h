#ifndef SYMMETRY_CPP_MATRIX_H
#define SYMMETRY_CPP_MATRIX_H

#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include <iterator>
#include <functional>
#include <complex>

typedef std::complex<double> dcomplex;
typedef std::complex<float> fcomplex;

template<class T>
class Matrix {
public:

    //set all elements to zero
    Matrix() : nrow(0), ncol(0)
    {}

    //allocate number of row col and elements
    Matrix(int nrow, int ncol) : nrow(nrow), ncol(ncol), data_(nrow * ncol)
    {}

    // initialize with matrix elements in col major
    Matrix(int nrow, int ncol, std::vector<T>& v) : nrow(nrow), ncol(ncol), data_(nrow * ncol)
    {data_ = v;}

    // copy constructor
    Matrix(const Matrix<T>& m) {
        nrow=m.nrow;
        ncol=m.ncol;
        data_=m.data_;
    }

    const T& operator()(int i,int j) const {
        //    std::cout << i << " " << j << std::endl;
        assert(i < nrow && j < ncol);
        assert(i+ j * nrow < data_.size());
        return data_[i+ j * nrow];
    }

    T& operator()(int i, int j) {
        //    std::cout << i << " " << j << std::endl;
        assert(i < nrow && j < ncol);
        assert(i+ j * nrow < data_.size());
        return data_[i+ j * nrow];
    }  // caller operator. returns a reference of that element

    void print() {
            std::cout.precision(8);
            std::cout << "shape:= (" << nrow << "," << ncol << ")" << std::endl;
            for(int i=0; i < nrow; i++) {
                for(int j=0; j < ncol; j++) {
                    std::cout << data_[i+ j * nrow] << "\t";
                }
                std::cout << std::endl;
            }
    }

    void clear(){
        for(int i=0; i < nrow; i++) {
            for(int j=0; j < ncol; j++) {
                data_[i+ j * nrow] = 0.0;
            }
        }
    }

    void resize(int newrow, int newcol) {
        assert(newrow>0 && newcol>0);
        nrow=newrow;
        ncol=newcol;
        data_.clear();
        data_.resize(newrow*newcol);
    }

    inline int rows() {
        return nrow;
    }

    inline int cols() {
        return ncol;
    }

    void fill(T val) {
        std::fill(data_.begin(),data_.end(),val);
    }

    void fillRand() {
        // First create an instance of an engine.
        std::random_device rnd_device;
        // Specify the engine and distribution.
        std::mt19937 mersenne_engine {rnd_device()};  // Generates random
        std::uniform_real_distribution<double> dist {0, 1};

        auto gen = [&dist, &mersenne_engine](){
            return dist(mersenne_engine);
        };
        generate(data_.begin(), data_.end(), gen);
    } // https://stackoverflow.com/a/23143753/14853469

    int numNonZeros(){
        int counter=0;
        for(int i=0; i<data_.size(); i++)
            if(data_[i]!=0.0) counter++;
        return counter;
    }
    void del(){
        nrow=0; ncol=0;
        data_.resize(0);
    }

    bool IsSquare() {
        bool out = true;
        if (nrow != ncol) {
            out = false;
            return out;
        }
        return out;
    }

    bool IsHermitian() {
        bool out=true;
        if (not IsSquare()) {out = false; return out;}
        for(int i=0; i < nrow; i++) {
            for(int j=0; j < ncol; j++) {
                std::complex<T> Hij = data_[i+ j * nrow];
                std::complex<T> Hji = data_[j+ i * nrow];
                if(Hij != std::conj(Hji)) {
                    std::string tmp = "Hij != Hji " + std::to_string(i)+"-"+ std::to_string(j)+" \n";
                    std::cout << i << " \t " << j << " \t " << Hij << " \t " << Hji << std::endl;
                    out=false;
                    return out;
                }
            }
        }
        return out;
    }

    bool IsSymmetric() {
        bool out = true;
        if (not IsSquare()) {out = false; return out;}
        for(int i=0; i < nrow; i++) {
            for(int j=0; j < ncol; j++) {
                T Hij = data_[i+ j * nrow];
                T Hji = data_[j+ i * nrow];
                if(Hij != Hji) {
                    std::string tmp = "Hij != Hji " + std::to_string(i)+"-"+ std::to_string(j)+" \n";
                    std::cout << i << " \t " << j << " \t " << Hij << " \t " << Hji << std::endl;
                    out = false;
                    return out;
                }
            }
        }
        return out;
    }

    void conjugate()
    {
        int n = data_.size();
        for (int i = 0; i < n; ++i)
            data_[i] = std::conj(data_[i]);
    }

    void transpose(Matrix<T>& m2, const Matrix<T>& m);
    void transpose();
    void ajoint(Matrix<T>& m2, const Matrix<T>& m);
    void ajoint();

private:
    int nrow, ncol;
    std::vector<T> data_;
};

// functions
template<class T>
void Matrix<T>::transpose(Matrix<T>& m2,const Matrix<T>& m)
{
    m2.resize(m.cols(),m.rows());
    for (int i=0;i<m2.rows();++i)
        for (int j=0;j<m2.cols();++j)
            m2(i,j) = m(j,i);
}

template<class T>
void Matrix<T>::transpose() {
    if (IsSquare()) {
        for (int i = 0; i < nrow; ++i) {
            for (int j = 0; j < i; ++j) {
                T tmp = data_[i+ j * nrow];
                data_[i+ j * nrow] = data_[j+ i * nrow];
                data_[j+ i * nrow] = tmp;
            }
        }
    }
    else {
        int newnrow = ncol; int newncol = nrow;
        Matrix<T> m2(newnrow, newncol);
        for (int i=0; i<newnrow; ++i)
            for (int j=0; j<newncol; ++j)
                m2(i,j) = data_[j+ i * nrow];

        this->resize(newnrow, newncol);
        for (int i=0; i<newnrow; ++i)
            for (int j=0; j<newncol; ++j)
                data_[i+ j * nrow] = m2(i, j);
    }
}

template<class T>
void Matrix<T>::ajoint(Matrix<T>& m2, const Matrix<T>& m){
    m2(m);
    m2.conjugate();
    m2.transpose();
}

template<class T>
void Matrix<T>::ajoint(){
    this->conjugate();
    this->transpose();
}

// These go to Matrix.cpp
// real symmetric and Hermitian matrix
void diag(Matrix<dcomplex> &A, std::vector<double> &evals, char option);
void diag(Matrix<fcomplex> &A, std::vector<float> &evals, char option);
void diag(Matrix<double> &A, std::vector<double> &evals, char option);
void diag(Matrix<float> &A, std::vector<float> &evals, char option);

// non-symmetric real and complex matrix
void diag(Matrix<dcomplex> &m, std::vector<dcomplex>& evals, std::vector<dcomplex>& vr, char option);
void diag(Matrix<fcomplex> &m, std::vector<fcomplex>& evals, std::vector<fcomplex>& vr, char option);
void diag(Matrix<double> &m, std::vector<double>& evalsRe, std::vector<double>& evalsIm,
          std::vector<double> vr, char option);
void diag(Matrix<float> &m, std::vector<float>& evalsRe, std::vector<float>& evalsIm,
          std::vector<float> vr, char option);

#endif
