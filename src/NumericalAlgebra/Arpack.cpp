//
// Created by huangkai on 25-1-15.
//

#include "Arpack.hpp"

ArpackSolver::ArpackSolver(int N, int nev)
: N_(N), nev_(nev), ncv_(2 * nev + 1), ldv_(N), ldz_(N),
        lworkl_(ncv_ * (ncv_ + 8)), resid_(N), V_(ldv_ * ncv_), z_(ldz_ * nev_),
        d_(nev_), workd_(3 * N), workl_(lworkl_), select_(ncv_), iparam_{},
        ipntr_{}, tol_(0.0)
{
  if (ncv_ <= nev_) {
             throw std::invalid_argument("NCV must be greater than NEV.");
      }
        // Initialize iparam
        iparam_[0] = 1;      // Use exact shifts
        iparam_[2] = 10 * N; // Maximum iterations
        iparam_[3] = 1;      // Block size (fixed to 1)

    ritz_option=arpack::which::smallest_algebraic;
}



void ArpackSolver::SetOption(arpack::which RitzOption){
    ritz_option=RitzOption;
}



void ArpackSolver::runRealSymmetric(const SparseMatrix &OP)
{
    iparam_[6] = 1; // Standard eigenproblem
    int info = 0, ido = 0;
    do {
        arpack::saupd(ido, arpack::bmat::identity, N_, ritz_option, nev_, tol_,
                      resid_.givePointer(), ncv_, V_.givePointer(), ldv_, iparam_, ipntr_,
                      workd_.givePointer(), workl_.givePointer(), lworkl_, info);

        if (ido == 1 || ido == -1) {
          OP.Mult(&(workd_[ipntr_[0] - 1]),&(workd_[ipntr_[1] - 1]));
        }
    } while (ido == 1 || ido == -1);
    if (info < 0) {
        std::cerr << "Error in ARPACK saupd: info = " << info << std::endl;
        throw std::runtime_error("ARPACK saupd failed.");
    }

    // Compute eigenvalues and vectors
    arpack::seupd(1, arpack::howmny::ritz_vectors, select_.givePointer(), d_.givePointer(),
                  z_.givePointer(), ldz_, 0.0, arpack::bmat::identity, N_, ritz_option,
                  nev_, tol_, resid_.givePointer(), ncv_, V_.givePointer(), ldv_, iparam_,
                  ipntr_, workd_.givePointer(), workl_.givePointer(), lworkl_, info);

    if (info < 0) {
        std::cerr << "Error in ARPACK seupd: info = " << info << std::endl;
        throw std::runtime_error("ARPACK seupd failed.");
    }
}

void ArpackSolver::runGeneralizedSymmetric(const PardisoSolver &solver, const SparseMatrix &B){
    iparam_[6] = 3; // Generalized eigenproblem
    FloatArray x1(N_), x2(N_);
    int info = 0, ido = 0;
    do {
        arpack::saupd(ido, arpack::bmat::generalized, N_, ritz_option, nev_, tol_,
                      resid_.givePointer(), ncv_, V_.givePointer(), ldv_, iparam_, ipntr_,
                      workd_.givePointer(), workl_.givePointer(), lworkl_, info);
        if (ido == -1) {
            cblas_dcopy(N_, &(workd_[ipntr_[0] - 1]), 1, x1.givePointer(), 1);
            B.Mult(x1, x2);
            solver.Mult(x1, x2);
            cblas_dcopy(N_, x1.givePointer(), 1, &(workd_[ipntr_[1] - 1]), 1);
        } else if (ido == 1) {
            cblas_dcopy(N_, &(workd_[ipntr_[2] - 1]), 1, &(workd_[ipntr_[1] - 1]),
                        1);
            cblas_dcopy(N_, &(workd_[ipntr_[1] - 1]), 1, x2.givePointer(), 1);
            solver.Mult(x2, x1);
            cblas_dcopy(N_, x1.givePointer(), 1, &(workd_[ipntr_[1] - 1]), 1);
        } else if (ido == 2) {
            cblas_dcopy(N_, &(workd_[ipntr_[0] - 1]), 1, x1.givePointer(), 1);
            B.Mult(x1, x2);
            cblas_dcopy(N_, x2.givePointer(), 1, &(workd_[ipntr_[1] - 1]), 1);
        }
    } while (ido == 1 || ido == -1 || ido == 2);
    if (info < 0) {
        std::cerr << "Error in ARPACK saupd: info = " << info << std::endl;
        throw std::runtime_error("ARPACK saupd failed.");
    }
    // Compute eigenvalues and vectors
    arpack::seupd(1, arpack::howmny::ritz_vectors, select_.givePointer(), d_.givePointer(),
                  z_.givePointer(), ldz_, 0.0, arpack::bmat::generalized, N_,
                  ritz_option, nev_, tol_, resid_.givePointer(), ncv_, V_.givePointer(), ldv_,
                  iparam_, ipntr_, workd_.givePointer(), workl_.givePointer(), lworkl_, info);
    if (info < 0) {
        std::cerr << "Error in ARPACK seupd: info = " << info << std::endl;
        throw std::runtime_error("ARPACK seupd failed.");
    }
}


void ArpackSolver::runGeneralizedNonSymmertricRight(const PardisoSolver &solver, const SparseMatrix &Stif,const SparseMatrix &Mass,
    const SparseMatrix &Damp) {
    ncv_=2*nev_+1;
    ldv_=N_*2;
    lworkl_=ncv_*(3*ncv_+6);
    resid_.resize(N_*2);
    V_.resize(ldv_*ncv_);
    z_.resize(N_*4*nev_);
    workd_.resize(6*N_);
    workl_.resize(lworkl_);

    select_.resize(ncv_);
    FloatArray workev_(3*ncv_);
    FloatArray dr(nev_*2);
    FloatArray di(nev_*2);

    resid_[0] = 1;
    iparam_[2] = 20*N_;
    iparam_[6] = 3; // mode

    EigenvaluesComplex_.resize(nev_);
    EigenvectorsComplex_.resize(N_,nev_);


    int iter=0;
    int info = 0, ido = 0;
    do {
        arpack::naupd(ido, arpack::bmat::generalized, N_* 2, ritz_option, nev_, tol_,
                      resid_.givePointer(), ncv_, V_.givePointer(), ldv_, iparam_, ipntr_,
                      workd_.givePointer(), workl_.givePointer(), lworkl_, info);
        if (ido == -1||ido ==1) {
            cblas_dcopy(N_, &(workd_[ipntr_[0]  - 1]), 1, &(workd_[ipntr_[1] +N_ -1]), 1);
            FloatArray work1(N_), work2(N_),temp1(N_),temp2(N_);
            cblas_dcopy(N_, &(workd_[ipntr_[0] - 1]), 1, work1.givePointer(), 1);
            cblas_dcopy(N_, &(workd_[ipntr_[0] +N_ -1]), 1, work2.givePointer(), 1);
            Damp.Mult(work1,temp1);
            Mass.Mult(work2,temp2);
            temp1=-1.0*temp1-temp2;
            solver.Mult(temp1,temp2);
            cblas_dcopy(N_,temp2.givePointer(),1, &(workd_[ipntr_[1]-1]), 1);
        }
        else if (ido == 2) {
            cblas_dcopy(N_, &(workd_[ipntr_[0] +N_ - 1]), 1, &(workd_[ipntr_[1] +N_-1]), 1);
            FloatArray work1(N_), work2(N_),temp1(N_),temp2(N_);
            cblas_dcopy(N_, &(workd_[ipntr_[0] - 1]), 1, temp1.givePointer(), 1);
            Stif.Mult(temp1,temp2);
            cblas_dcopy(N_,temp2.givePointer(),1, &(workd_[ipntr_[1]-1]), 1);
        }
        iter++;
    } while (ido == 1 || ido == -1 || ido == 2);


    if (info < 0) {
        std::cerr << "Error in ARPACK saupd: info = " << info << std::endl;
        throw std::runtime_error("ARPACK saupd failed.");
    }

    arpack::neupd(1, arpack::howmny::ritz_vectors, select_.givePointer(), dr.givePointer(), di.givePointer(), z_.givePointer(),
        ldv_, 0.0, 0.0, workev_.givePointer(), arpack::bmat::generalized, N_ * 2,
            ritz_option, nev_, tol_, resid_.givePointer(), ncv_, V_.givePointer(), ldv_,
            iparam_, ipntr_, workd_.givePointer(), workl_.givePointer(), lworkl_, info);
    if (info < 0) {
        std::cerr << "Error in ARPACK seupd: info = " << info << std::endl;
        throw std::runtime_error("ARPACK seupd failed.");
    }


    EigenvaluesComplex_.resize(nev_);
    for (int i = 0; i < nev_; i++) {
        ComplexScalar temp(dr[i],di[i]);
        EigenvaluesComplex_.Set(i,temp);
    }

    for (int j=0;j<nev_-nev_%2;j=j+2) {
        for (int i=0;i<N_;++i) {
            EigenvectorsComplex_.Set(i,j,z_[j*N_*2+i],z_[j*N_*2+2*N_+i]);
            EigenvectorsComplex_.Set(i,j+1,z_[j*N_*2+i],-z_[j*N_*2+2*N_+i]);
        }
    }

    if (nev_ %2 ==1) {
        for (int i=0;i<N_;++i) {
            EigenvectorsComplex_.Set(i,nev_-1,z_[(nev_-1)*N_*2+i],z_[(nev_-1)*N_*2+2*N_+i]);
        }
    }
}


void ArpackSolver::runGeneralizedNonSymmertricLeft(const PardisoSolver &solver, const SparseMatrix &Stif,const SparseMatrix &Mass,
    const SparseMatrix &Damp) {
        ncv_=2*nev_+1;
    ldv_=N_*2;
    lworkl_=ncv_*(3*ncv_+6);
    resid_.resize(N_*2);
    V_.resize(ldv_*ncv_);
    z_.resize(N_*4*nev_);
    workd_.resize(6*N_);
    workl_.resize(lworkl_);
    tol_=1e-18;

    select_.resize(ncv_);
    FloatArray workev_(3*ncv_);
    FloatArray dr(nev_*2);
    FloatArray di(nev_*2);

    resid_[0] = 1;
    iparam_[2] = 20*N_;
    iparam_[6] = 3; // mode

    EigenvaluesComplex_.resize(nev_);
    EigenvectorsComplex_.resize(N_,nev_);

    matrix_descr descr = { SPARSE_MATRIX_TYPE_GENERAL,SPARSE_FILL_MODE_FULL,SPARSE_DIAG_NON_UNIT };
    int iter=0;
    int info = 0, ido = 0;
    do {
        arpack::naupd(ido, arpack::bmat::generalized, N_* 2, ritz_option, nev_, tol_,
                      resid_.givePointer(), ncv_, V_.givePointer(), ldv_, iparam_, ipntr_,
                      workd_.givePointer(), workl_.givePointer(), lworkl_, info);
        if (ido == -1 || ido == 1) {
            FloatArray intermediate_vector(N_),temp(N_);
            mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, -1.0, *Damp.getMKLMatrix(),descr , &(workd_[ipntr_[0] - 1]), 0.0, intermediate_vector.givePointer());
            cblas_daxpy(N_, 1.0, &(workd_[ipntr_[0] +N_ - 1]), 1, intermediate_vector.givePointer(), 1);
            solver.Mult(intermediate_vector,temp);
            cblas_dcopy(N_, temp.givePointer(), 1, &workd_[ipntr_[1] - 1], 1);
            mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, -1.0, *Mass.getMKLMatrix(), descr, &workd_[ipntr_[0] - 1], 0.0, &workd_[ipntr_[1] + N_ - 1]);

        }
        else if (ido == 2) {
            cblas_dcopy(N_, &workd_[ ipntr_[0] + N_ - 1], 1, &workd_[ipntr_[1] + N_ - 1], 1);
            mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, 1.0, *Stif.getMKLMatrix(), descr, &workd_[ipntr_[0] - 1], 0.0, &workd_[ipntr_[1] - 1]);
        }
        iter++;
    } while (ido == 1 || ido == -1 || ido == 2);


    if (info < 0) {
        std::cerr << "Error in ARPACK saupd: info = " << info << std::endl;
        throw std::runtime_error("ARPACK saupd failed.");
    }

    arpack::neupd(1, arpack::howmny::ritz_vectors, select_.givePointer(), dr.givePointer(), di.givePointer(), z_.givePointer(),
        ldv_, 0.0, 0.0, workev_.givePointer(), arpack::bmat::generalized, N_ * 2,
            ritz_option, nev_, tol_, resid_.givePointer(), ncv_, V_.givePointer(), ldv_,
            iparam_, ipntr_, workd_.givePointer(), workl_.givePointer(), lworkl_, info);
    if (info < 0) {
        std::cerr << "Error in ARPACK seupd: info = " << info << std::endl;
        throw std::runtime_error("ARPACK seupd failed.");
    }


    EigenvaluesComplex_.resize(nev_);
    for (int i = 0; i < nev_; i++) {
        ComplexScalar temp(dr[i],di[i]);
        EigenvaluesComplex_.Set(i,temp);
    }

    for (int j=0;j<nev_-nev_%2;j=j+2) {
        for (int i=0;i<N_;++i) {
            EigenvectorsComplex_.Set(i,j,z_[j*N_*2+i],z_[j*N_*2+2*N_+i]);
            EigenvectorsComplex_.Set(i,j+1,z_[j*N_*2+i],-z_[j*N_*2+2*N_+i]);
        }
    }

    if (nev_ %2 ==1) {
        for (int i=0;i<N_;++i) {
            EigenvectorsComplex_.Set(i,nev_-1,z_[(nev_-1)*N_*2+i],z_[(nev_-1)*N_*2+2*N_+i]);
        }
    }
}





FloatArray ArpackSolver::getEigenvector(int i) const {
    FloatArray result(N_);
    for (int j = 0; j < N_; j++) {
        result[j]=z_[i*N_+j];
    }
    return result;
}

FloatMatrix ArpackSolver::getEigenvector(IntArray& idx) const {
    FloatMatrix result(N_, idx.giveSize());
    for (int i=0;i<idx.giveSize();i++) {
        for(int j=0;j<N_;j++) {
            result(j,i)=z_[idx[i]*N_+j];
        }
    }
    return result;
}


ComplexVector ArpackSolver::getComplexEigenvector(int i) const {
    ComplexVector result(N_);
    for (int j = 0; j < N_; j++) {
        result(j)=EigenvectorsComplex_.At(i,j);
    }
    return result;
}

ComplexMatrix ArpackSolver::getComplexEigenvector(IntArray& idx) const {
    ComplexMatrix result(N_, idx.giveSize());
    for (int i=0;i<idx.giveSize();i++) {
        for(int j=0;j<N_;j++) {
            result(j,i)=EigenvectorsComplex_.At(j,idx[i]);
        }
    }
    return result;
}

ComplexVector ArpackSolver::getComplexEigenvalue(IntArray& idx) const {
    ComplexVector result(idx.giveSize());
    for (int i=0;i<idx.giveSize();i++) {
        result(i)=EigenvaluesComplex_(idx(i));
    }
    return result;
}






void ArpackSolver::SetEigenvalue(int index, double value) {
    for (int i=0;i<nev_;i++) {
        for(int j=0;j<N_;j++) {
            z_[i*N_+index]=value;
        }
    }
}

void ArpackSolver::DevideValue(int index, double value) {
    for (int j=0;j<N_;j++) {
        z_[index*N_+j]=z_[index*N_+j]/value;
    }
}



FloatMatrix ArpackSolver::getEigenvector(std::vector<int>& idx)const {
    FloatMatrix result(N_, idx.size());
    for (int i=0;i<idx.size();i++) {
        for(int j=0;j<N_;j++) {
            result(j,i)=z_[idx[i]*N_+j];
        }
    }
    return result;
}

void ArpackSolver::changeEigenvector(IntArray index, FloatMatrix values) {
    for (int idx=0;idx<index.giveSize();idx++) {
        for (int j=0;j<N_;j++) {
            z_[index(idx)*N_+j]=values(j,idx);
        }
    }
}

