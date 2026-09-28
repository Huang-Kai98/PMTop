// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

//
// Created by huangkai on 25-1-15.
//

#ifndef PARDISO_H
#define PARDISO_H
#include "floatarray.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "operator.hpp"
#include <iostream>
#include <mkl.h>

#include "sparsemat.hpp"


/**
 * @brief MKL Parallel Direct Sparse Solver PARDISO
 *
 * Interface to MKL PARDISO: the direct sparse solver based on PARDISO
 */

class PardisoSolver
{
public:
    enum MatType
    {
        REAL_STRUCTURE_SYMMETRIC = 1,
        REAL_SYMMETRIC_POSITIVE_DEFINITE = 2,
        REAL_SYMMETRIC_INDEFINITE = -2,
        REAL_NONSYMMETRIC = 11,
        COMPLEX_NONSYMMETRIC=13
     };

    /**
    * @brief Construct a new PardisoSolver object
    *
    */
    PardisoSolver();

    /**
    * @brief Set the Operator object and perform factorization
    *
    * @a op needs to be of type SparseMatrix.
    *
    * @param op Operator to use in factorization and solve
    */
    void SetOperator(const Operator &op);

    void SetOperator(const SparseMatrixComplex& op);

    /**
    * @brief Solve
    *
    * @param b RHS vector
    * @param x Solution vector
    */
    void Mult(const FloatArray &b, FloatArray &x) const;

    void Mult(double *b, double *x) const;

    void Mult(ComplexVector &b, ComplexVector &x) const;

    /**
    * @brief Set the print level for MKL Pardiso
    *
    * Prints statistics after the factorization and after each solve.
    *
    * @param print_lvl Print level
    */
    void SetPrintLevel(int print_lvl);

    /**
    * @brief Set the matrix type
    *
    * The matrix type supported is either real and symmetric or real and
    * non-symmetric.
    *
    * @param mat_type Matrix type
    */
    void SetMatrixType(MatType mat_type);

    ~PardisoSolver();
private:
    int width;  ///< Dimension of the input / number of columns in the matrix.
    // Global number of rows
    int m;

    // Number of nonzero entries
    int nnz;

    // CSR data structure for the copy data of the local CSR matrix
    int *csr_rowptr = nullptr;
    double *reordered_csr_nzval = nullptr;
    int *reordered_csr_colind = nullptr;

    MKL_Complex16 *reordered_csr_val = nullptr;

    // Internal solver memory pointer pt,
    // 32-bit: int pt[64]
    // 64-bit: long int pt[64] or void *pt[64] should be OK on both architectures
    mutable void *pt[64] = {0};

    // Solver control parameters, detailed description can be found in the
    // constructor.
    mutable int iparm[64] = {0};
    mutable int maxfct, mnum, msglvl, phase, error;
    int mtype;
    int nrhs;

    // Dummy variables
    mutable int idum;
    mutable double ddum;

};

#endif //PARDISO_H
