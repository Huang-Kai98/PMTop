// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

//
// Created by huangkai on 25-1-15.
//

#ifndef ARPACK_H
#define ARPACK_H

#include "pardiso.hpp"
#include "arpack.hpp"
#include "floatarray.h"
#include "floatmatrix.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <mkl_spblas.h>
#include <stdexcept>
#include <vector>

#include "ComplexMatrix.h"
#include "ComplexVector.hpp"
#include "sparsemat.hpp"
#include "Logger.h"

/**
 * @brief Arpack
 *
 * ARPACK is a collection of Fortran77 subroutines designed to solve large scale eigenvalue problems.
 */

class ArpackSolver {
public:
    ArpackSolver(int N, int nev);

    /**
    * @brief Set the RitzOption
    *
    * @a RitzOption needs to be of type arpack::which.
    * largest_algebraic,smallest_algebraic,largest_magnitude,
    * smallest_magnitude,largest_real,smallest_real,largest_imaginary,
    * smallest_imaginary,both_ends
    *
    */
    void SetOption(arpack::which RitzOption);

    /**
    * @brief Compute Standard Eigen Problem
    *
    * @a solver need Pardiso factorization , OP SparseMatrix
    */
    void runRealSymmetric(const SparseMatrix &OP);

    /**
    * @brief Compute Generalized Eigen Problem
    *
    * @a solver need Pardiso factorization , OP SparseMatrix
    */
    void runGeneralizedSymmetric(const PardisoSolver &solver, const SparseMatrix &B);

    void runGeneralizedNonSymmertricRight(const PardisoSolver &solver, const SparseMatrix &Stif, const SparseMatrix &Mass,
        const SparseMatrix &Damp);

    void runGeneralizedNonSymmertricLeft(const PardisoSolver &solver, const SparseMatrix &Stif, const SparseMatrix &Mass,
        const SparseMatrix &Damp);


    /**
    * @brief Compute Eigen Problem
    *
    * @a solver need Pardiso factorization , OP SparseMatrix
    */
    void Mult(const SparseMatrix &OP) {
        runRealSymmetric(OP);
    }
    void Mult(const PardisoSolver &solver, const SparseMatrix &B) {
        runGeneralizedSymmetric(solver, B);
    };

    void Mult(const PardisoSolver &solver, const SparseMatrix &A, const SparseMatrix &B, const SparseMatrix &C, std::string SolutionOpt) {
        if (SolutionOpt == "Right") {
            runGeneralizedNonSymmertricRight(solver,A,B,C);
        }
        else if (SolutionOpt == "Left") {
            runGeneralizedNonSymmertricLeft(solver,A,B,C);
        }
        else {
            LOG_ERROR("Nukown SolutionOpt");
        }
    }


    FloatArray getEigenvalues() const { return d_; }
    FloatArray getEigenvectors() const { return z_; }
    FloatArray getEigenvector(int i) const;
    FloatMatrix getEigenvector(IntArray &idx) const;
    ComplexVector getComplexEigenvalues() const {return EigenvaluesComplex_;};
    ComplexMatrix getComplexEigenvectors() const {return EigenvectorsComplex_;};
    ComplexVector getComplexEigenvector(int i) const;
    ComplexMatrix getComplexEigenvector(IntArray &idx) const;
    ComplexVector getComplexEigenvalue(IntArray &idx) const;

    void SetEigenvalue(int index,double value);
    void DevideValue(int index,double value);
    int GetNev(){return nev_;};
    FloatMatrix getEigenvector(std::vector<int> &idx) const;


    void changeEigenvector(IntArray index, FloatMatrix values);

private:
    int N_;
    int nev_;
    int ncv_;
    int ldv_;
    int ldz_;
    int lworkl_;

    FloatArray resid_;
    FloatArray V_;
    FloatArray z_;
    FloatArray d_;
    FloatArray workd_;
    FloatArray workl_;
    IntArray select_;

    ComplexVector EigenvaluesComplex_;
    ComplexMatrix EigenvectorsComplex_;

    int iparam_[11];
    int ipntr_[14];
    double tol_;

    arpack::which ritz_option;
};



#endif //ARPACK_H
