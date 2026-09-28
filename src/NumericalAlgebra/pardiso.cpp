//
// Created by huangkai on 25-1-15.
//

#include "pardiso.hpp"
#include "sparsemat.hpp"
#include "error.h"

PardisoSolver::PardisoSolver()
{
    // Indicate that default parameters are changed
    iparm[0] = 1;
    // Use METIS for fill-in reordering
    iparm[1] = 2;
    // Write the solution into the x vector data
    iparm[5] = 0;
    // Maximum number of iterative refinement steps
    iparm[7] = 2;
    // Perturb the pivot elements with 1E-13
    iparm[9] = 13;
    // Use nonsymmetric permutation
    iparm[10] = 1;
    // Perform a check on the input data
    iparm[26] = 1;
#ifdef PMTOP_USE_SINGLE
    // Single precision
    iparm[27] = 1;
#endif
    // 0-based indexing in CSR data structure
    iparm[34] = 1;
    // Maximum number of numerical factorizations
    maxfct = 1;
    // Which factorization to use. This parameter is ignored and always assumed
    // to be equal to 1. See MKL documentation.
    mnum = 1;
    // Print statistical information in file
    msglvl = 0;
    // Initialize error flag
    error = 0;
    // Real nonsymmetric matrix
    mtype = MatType::REAL_NONSYMMETRIC;
    // Number of right hand sides
    nrhs = 1;
}

void PardisoSolver::SetOperator(const SparseMatrixComplex &op) {

    auto mat = const_cast<SparseMatrixComplex *>(&op);
    width = mat->Cols();
    m = mat->Rows();
    nnz = mat->NumNonZeroElems();
    const int *Ap = mat->GetI().givePointer();
    const int *Ai = mat->GetJ().givePointer();
    const MKL_Complex16 *Ax = mat->GetA().Data();
    csr_rowptr = new int[m + 1];
    reordered_csr_colind = new int[nnz];
    reordered_csr_val = new MKL_Complex16[nnz];
    for (int i = 0; i <= m; i++)
    {
        csr_rowptr[i] = Ap[i];
    }
    mat->SortColumnIndices();
    for (int i = 0; i < nnz; i++)
    {
        reordered_csr_colind[i] = Ai[i];
        reordered_csr_val[i] = Ax[i];
    }
    // Analyze inputs
    phase = 11;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &m, reordered_csr_val, csr_rowptr,
            reordered_csr_colind, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);

    if (error != 0) {
        throw std::runtime_error("Error during symbolic factorization");
    }

    // Numerical factorization
    phase = 22;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &m, reordered_csr_val, csr_rowptr,
            reordered_csr_colind, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);

    if (error != 0) {
        throw std::runtime_error("Error during numerical factorization");
    }

}



void PardisoSolver::SetOperator(const Operator &op)
{
    auto mat = const_cast<SparseMatrix *>(dynamic_cast<const SparseMatrix *>(&op));

    width = mat->Width();

    m = mat->Size();

    nnz = mat->NumNonZeroElems();

    const int *Ap = mat->GetI();
    const int *Ai = mat->GetJ();
    const double *Ax = mat->GetData();

    csr_rowptr = new int[m + 1];
    reordered_csr_colind = new int[nnz];
    reordered_csr_nzval = new double[nnz];

    for (int i = 0; i <= m; i++)
    {
        csr_rowptr[i] = Ap[i];
    }
    mat->SortColumnIndices();


    for (int i = 0; i < nnz; i++)
    {
        reordered_csr_colind[i] = Ai[i];
        reordered_csr_nzval[i] = Ax[i];
    }

    // Analyze inputs
    phase = 11;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &m, reordered_csr_nzval, csr_rowptr,
            reordered_csr_colind, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);

    if (error != 0) {
        throw std::runtime_error("Error during symbolic factorization");
    }

    // Numerical factorization
    phase = 22;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &m, reordered_csr_nzval, csr_rowptr,
            reordered_csr_colind, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);

    if (error != 0) {
        throw std::runtime_error("Error during numerical factorization");
    }

}

void PardisoSolver::Mult(const FloatArray &b, FloatArray &x) const
{
    if (b.giveSize() != static_cast<size_t>(m)) {
        throw std::invalid_argument(
            "Size of b must match the number of rows in the matrix");
    }
    // Solve
    phase = 33;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &m, reordered_csr_nzval, csr_rowptr,
            reordered_csr_colind, &idum, &nrhs,
            iparm, &msglvl, const_cast<double *>(b.givePointer()), x.givePointer(), &error);

    if (error != 0) {
        throw std::runtime_error("Error during solution");
    }

}

void PardisoSolver::Mult(double *b, double *x) const {
    // Solve
    phase = 33;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &m, reordered_csr_nzval, csr_rowptr,
            reordered_csr_colind, &idum, &nrhs,
            iparm, &msglvl, b, x, &error);

    if (error != 0) {
        throw std::runtime_error("Error during solution");
    }
}

void PardisoSolver::Mult(ComplexVector &b, ComplexVector &x) const {
    // Solve
    phase = 33;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &m, reordered_csr_val, csr_rowptr,
            reordered_csr_colind, &idum, &nrhs,
            iparm, &msglvl, b.Data(), x.Data(), &error);

    if (error != 0) {
        throw std::runtime_error("Error during solution");
    }
}




void PardisoSolver::SetPrintLevel(int print_level)
{
    msglvl = print_level;
}

void PardisoSolver::SetMatrixType(MatType mat_type)
{
    mtype = mat_type;
}

PardisoSolver::~PardisoSolver()
{
    // Release all internal memory
    phase = -1;
    if (reordered_csr_nzval != nullptr) {
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &m, reordered_csr_nzval, csr_rowptr,
                reordered_csr_colind, &idum, &nrhs,
                iparm, &msglvl, &ddum, &ddum, &error);
        delete[] reordered_csr_nzval;
    }
    if (reordered_csr_val!=nullptr) {
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &m, reordered_csr_val, csr_rowptr,
                reordered_csr_colind, &idum, &nrhs,
                iparm, &msglvl, &ddum, &ddum, &error);
        delete[] reordered_csr_val;
    }

    delete[] csr_rowptr;
    delete[] reordered_csr_colind;
}