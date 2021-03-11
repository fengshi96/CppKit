#include "Matrix.h"

typedef std::complex<double> dcomplex;
typedef std::complex<float> fcomplex;

extern "C" {
    void zheev_(char *,char *,int *, dcomplex *, int *, double *, dcomplex *,int *, double *, int *);
    void cheev_(char *,char *,int *, fcomplex *, int *, float *, fcomplex *,int *, float *, int *);
    void dsyev_(char *,char *,int *,double *, int *, double *, double *, int *, int *);
    void ssyev_(char *,char *,int *,float *, int *, float *, float *, int *, int *);
    void dgeev_(char*, char*, int*, double *, int*, double* , double*, double*, int*, double*, int*, double*, int*, int*);
    void sgeev_(char*, char*, int*, float *, int*, float* , float*, float*, int*, float*, int*, float*, int*, int*);
    void zgeev_(char*, char*, int*, dcomplex *, int*, dcomplex* , dcomplex*, int*, dcomplex*,
                int*, dcomplex*, int*, double*, int*);
    void cgeev_(char*, char*, int*, fcomplex *, int*, fcomplex* , fcomplex*, int*, fcomplex*,
            int*, fcomplex*, int*, float*, int*);
}

// diag double precision Hermitian matrix m
void diag(Matrix<dcomplex>& m, std::vector<double>& evals, char option){
    char jobz=option;  // 'N':  Compute eigenvalues only; 'V':  Compute eigenvalues and eigenvectors.
    char uplo='U';  // 'U':  Upper triangle of A is stored;
    int n=m.rows();
    int lda=m.cols();
    assert(n==lda);
    assert(m.IsHermitian());

    evals.resize(n);
    std::vector<dcomplex> work(3);
    std::vector<double> rwork(3*n-2);
    int info,lwork= -1;  // If LWORK = -1, then a workspace query is assumed

    fill(evals.begin(),evals.end(),0);

    // query:
    zheev_(&jobz,&uplo,&n,&(m(0,0)),&lda,&(evals[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<" during query \n";
        perror("diag: zheev: failed with info!=0.\n");
    }
    lwork = static_cast<int>(std::real(work[0]))+1;
    work.resize(lwork+1);

    // real work:
    zheev_(&jobz,&uplo,&n,&(m(0,0)),&lda,&(evals[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"work info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }
    //for(int i=0;i<n;i++) { cout << eval[i] << " \t ";} cout << endl;
    //sort(eigs_.begin(),eigs_.end()); // sort Eigenvalues and Hamiltonian
}

// diag single precision Hermitian matrix m
void diag(Matrix<std::complex<float>> &m, std::vector<float>& evals, char option)
{
    char jobz=option;
    char uplo='U';
    int n=m.rows();
    int lda=m.cols();
    assert(n==lda);
    assert(m.IsHermitian());

    std::vector<std::complex<float> > work(3);
    std::vector<float> rwork(3*n);
    int info,lwork= -1;

    evals.resize(n);

    // query:
    cheev_(&jobz,&uplo,&n,&(m(0,0)),&lda,&(evals[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<" during query\n";
        throw std::runtime_error("diag: cheev_: failed with info!=0.\n");
    }

    const int NB = 256;
    lwork = std::max(1 + static_cast<int>(std::real(work[0])), (NB + 2)*n);
    work.resize(lwork);

    // real work:
    cheev_(&jobz,&uplo,&n,&(m(0,0)),&lda,&(evals[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"work info="<<info<<"\n";
        throw std::runtime_error("diag: cheev: failed with info!=0.\n");
    }
}


// diag double precision real symmetric matrix
void diag(Matrix<double> &m, std::vector<double>& evals, char option)
{
    char jobz=option;
    char uplo='U';
    int n=m.rows();
    int lda=m.cols();
    assert(n==lda);
    assert(m.IsSymmetric());

    std::vector<double> work(3);
    int info;
    int lwork= -1;

    if (lda<=0) throw std::runtime_error("lda<=0\n");

    evals.resize(n);

    // query:
    dsyev_(&jobz,&uplo,&n,&(m(0,0)),&lda, &(evals[0]),&(work[0]),&lwork, &info);
    if (info!=0) {
        std::cerr<<"info="<<info<<" during query\n";
        throw std::runtime_error("diag: dsyev_: failed with info!=0.\n");
    }

    const int NB = 256;
    lwork = std::max(1 + static_cast<int>(work[0]), (NB + 2)*n);
    work.resize(lwork);
    // real work:
    dsyev_(&jobz,&uplo,&n,&(m(0,0)),&lda, &(evals[0]),&(work[0]),&lwork, &info);
    if (info!=0) {
        std::cerr<<"work info="<<info<<"\n";
        throw std::runtime_error("diag: dsyev_: failed with info!=0.\n");
    }
}

// diag single precision real symmetric matrix
void diag(Matrix<float> &m, std::vector<float>& evals,char option)
{
    char jobz=option;
    char uplo='U';
    int n=m.rows();
    int lda=m.cols();
    assert(n==lda);
    assert(m.IsSymmetric());

    std::vector<float> work(3);
    int info;
    int lwork= -1;

    if (lda<=0) throw std::runtime_error("lda<=0\n");

    evals.resize(n);

    // query:
    ssyev_(&jobz,&uplo,&n,&(m(0,0)),&lda, &(evals[0]),&(work[0]),&lwork, &info);
    if (info!=0) {
        std::cerr<<"info="<<info<<" during query\n";
        throw std::runtime_error("diag: ssyev_: failed with info!=0.\n");
    }

    const int NB = 256;
    lwork = std::max(1 + static_cast<int>(work[0]), (NB + 2)*n);
    work.resize(lwork);

    // real work:
    ssyev_(&jobz,&uplo,&n,&(m(0,0)),&lda, &(evals[0]),&(work[0]),&lwork, &info);
    if (info!=0) {
        std::cerr<<"work info="<<info<<"\n";
        throw std::runtime_error("diag: ssyev_: failed with info!=0.\n");
    }
}

// diag double precision real non-symmetric matrix
void diag(Matrix<double> &m, std::vector<double>& evalsRe, std::vector<double>& evalsIm,
          std::vector<double> vr, char option){
    char jobvl = 'N';
    char jobvr = option;
    int n=m.rows();
    int lda=m.cols();
    assert(n==lda);

    std::vector<double> vl;
    vr.resize(n * n);

    int ldvl = 1;
    int ldvr = n;
    std::vector<double> work(3);
    int lwork= -1;
    int info;

    evalsRe.resize(n);
    evalsIm.resize(n);

    // query:
    dgeev_(&jobvl, &jobvr, &n, &(m(0,0)), &lda, &(evalsRe[0]), &(evalsIm[0]),
           &(vl[0]), &ldvl, &(vr[0]), &ldvr, &(work[0]), &lwork, &info);
    if (info!=0) {
        std::cerr<<"info="<<info<<" during query\n";
        throw std::runtime_error("diag: dgeev_: failed with info!=0.\n");
    }
    lwork = static_cast<int>(work[0])+1;
    work.resize(lwork);

    // real work:
    dgeev_(&jobvl, &jobvr, &n, &(m(0,0)), &lda, &(evalsRe[0]), &(evalsIm[0]),
           &(vl[0]), &ldvl, &(vr[0]), &ldvr, &(work[0]), &lwork, &info);
    if (info!=0) {
        std::cerr<<"work info="<<info<<"\n";
        throw std::runtime_error("diag: dgeev_: failed with info!=0.\n");
    }
}

// diag single precision real non-symmetric matrix
void diag(Matrix<float> &m, std::vector<float>& evalsRe, std::vector<float>& evalsIm,
          std::vector<float> vr, char option){
    char jobvl = 'N';
    char jobvr = option;
    int n=m.rows();
    int lda=m.cols();
    assert(n==lda);

    std::vector<float> vl;
    vr.resize(n * n);

    int ldvl = 1;
    int ldvr = n;
    std::vector<float> work(3);
    int lwork= -1;
    int info;

    evalsRe.resize(n);
    evalsIm.resize(n);

    // query:
    sgeev_(&jobvl, &jobvr, &n, &(m(0,0)), &lda, &(evalsRe[0]), &(evalsIm[0]),
           &(vl[0]), &ldvl, &(vr[0]), &ldvr, &(work[0]), &lwork, &info);
    if (info!=0) {
        std::cerr<<"info="<<info<<" during query\n";
        throw std::runtime_error("diag: sgeev_: failed with info!=0.\n");
    }
    lwork = static_cast<int>(work[0])+1;
    work.resize(lwork);

    // real work:
    sgeev_(&jobvl, &jobvr, &n, &(m(0,0)), &lda, &(evalsRe[0]), &(evalsIm[0]),
           &(vl[0]), &ldvl, &(vr[0]), &ldvr, &(work[0]), &lwork, &info);
    if (info!=0) {
        std::cerr<<"work info="<<info<<"\n";
        throw std::runtime_error("diag: sgeev_: failed with info!=0.\n");
    }
}

// diag double precision complex non-symmetric matrix
void diag(Matrix<dcomplex> &m, std::vector<dcomplex>& evals, std::vector<dcomplex>& vr, char option){
    char jobvl = 'N';
    char jobvr = option;
    int n=m.rows();
    int lda=m.cols();
    assert(n==lda);

    int ldvl = 1;
    int ldvr = n;
    std::vector<dcomplex> vl;
    evals.resize(n);
    vr.resize(ldvr * n);
    std::vector<dcomplex> work(3);
    std::vector<double> rwork(2 * n);
    int lwork= -1;
    int info;

    // query:
    zgeev_(&jobvl, &jobvr, &n, &(m(0,0)), &lda, &(evals[0]), &(vl[0]), &ldvl, &(vr[0]), &ldvr,
           &(work[0]), &lwork, &(rwork[0]), &info);
    if (info!=0) {
        std::cerr<<"info="<<info<<" during query\n";
        throw std::runtime_error("diag: zgeev_: failed with info!=0.\n");
    }
    const int NB = 256;
    lwork = std::max(1 + static_cast<int>(std::real(work[0])), (NB + 2)*n);
    work.resize(lwork);

    // real work:
    zgeev_(&jobvl, &jobvr, &n, &(m(0,0)), &lda, &(evals[0]), &(vl[0]), &ldvl, &(vr[0]), &ldvr,
           &(work[0]), &lwork, &(rwork[0]), &info);
    if (info!=0) {
        std::cerr<<"work info="<<info<<"\n";
        throw std::runtime_error("diag: zgeev_: failed with info!=0.\n");
    }
}


// diag single precision complex non-symmetric matrix
void diag(Matrix<fcomplex> &m, std::vector<fcomplex>& evals, std::vector<fcomplex>& vr, char option){
    char jobvl = 'N';
    char jobvr = option;
    int n=m.rows();
    int lda=m.cols();
    assert(n==lda);

    int ldvl = 1;
    int ldvr = n;
    std::vector<fcomplex> vl;
    evals.resize(n);
    vr.resize(ldvr * n);
    std::vector<fcomplex> work(3);
    std::vector<float> rwork(2 * n);
    int lwork= -1;
    int info;

    // query:
    cgeev_(&jobvl, &jobvr, &n, &(m(0,0)), &lda, &(evals[0]), &(vl[0]), &ldvl, &(vr[0]), &ldvr,
           &(work[0]), &lwork, &(rwork[0]), &info);
    if (info!=0) {
        std::cerr<<"info="<<info<<" during query\n";
        throw std::runtime_error("diag: cgeev_: failed with info!=0.\n");
    }
    const int NB = 256;
    lwork = std::max(1 + static_cast<int>(std::real(work[0])), (NB + 2)*n);
    work.resize(lwork);

    // real work:
    cgeev_(&jobvl, &jobvr, &n, &(m(0,0)), &lda, &(evals[0]), &(vl[0]), &ldvl, &(vr[0]), &ldvr,
           &(work[0]), &lwork, &(rwork[0]), &info);
    if (info!=0) {
        std::cerr<<"work info="<<info<<"\n";
        throw std::runtime_error("diag: cgeev_: failed with info!=0.\n");
    }
}