#pragma once
#include "floatarray.h"
#include "matrix.hpp"

/** Conjugate gradient method. Given Matrix A, vector b and initial guess
    x, iteratively solve A x = b. When the default arguments are used
    CG doesn't print current residual and number of iterations, maximum
    number of iterations is 1000, the relative tolerance is 10e-12 and
    the absolute tolerance is 10e-24. */
void CG(const Operator &A, const FloatArray &b, FloatArray &x,
        int print_iter = 0, int max_num_iter = 1000, double RTOLERANCE = 10e-12,
        double ATOLERANCE = 10e-24);

/** Preconditioned conjugate gradient method. Given Matrix A, preconditioner
    Matrix B, vector b and initial guess x, iteratively solve A x = b.
    When the default arguments are used PCG doesn't print current residuals
    and number of iterations, maximum number of iterations is 1000, the
    relative tolerance is 10e-12 and the absolute tolerance is 10e-24.
    Remark : if no better initial guess is available, the user may set
    it as B b (since not done in PCG routine). */
void PCG(const Operator &A, const Operator &B, const FloatArray &b,
         FloatArray &x, int print_iter = 0, int max_num_iter = 1000,
         double RTOLERANCE = 10e-12, double ATOLERANCE = 10e-24, int save = 0);

/// A GMRES solver
int GMRES(const Operator &A, FloatArray &x, const FloatArray &b,
          const Operator &M, int &max_iter, int m, double &tol, double &atol,
          int printit);

/** Adaptive restarted GMRES.
    m_max and m_min(=1) are the maximal and minimal restart parameters.
    m_step(=1) is the step to use for going from m_max and m_min.
    cf(=0.4) is a desired convergance factor. */
int aGMRES(const Operator &A, FloatArray &x, const FloatArray &b,
           const Operator &M, int &max_iter, int m_max, int m_min, int m_step,
           double cf, double &tol, double &atol, int printit);

/// A BiCG-Stab solver
int BiCGSTAB(const Operator &A, FloatArray &x, const FloatArray &b,
             const Operator &M, int &max_iter, double &tol, double atol,
             int printit);

/// Stationary linear iteration: x <- x + B (b - A x)
void SLI(const Operator &A, const Operator &B, const FloatArray &b,
         FloatArray &x, int print_iter = 0, int max_num_iter = 1000,
         double RTOLERANCE = 10e-12, double ATOLERANCE = 10e-24);