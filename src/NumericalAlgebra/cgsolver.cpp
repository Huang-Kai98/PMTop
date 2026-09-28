#include "floatarray.h"
#include "matrix.hpp"
#include "sparsemat.hpp"
#include <iomanip>
#include <iostream>

// Conjugate Gradient solver

void CG(const Operator &A, const FloatArray &b, FloatArray &x,
        int print_iter = 0, int max_num_iter = 1000, double RTOLERANCE = 10e-12,
        double ATOLERANCE = 10e-24) {

  int i, dim = x.giveSize();
  double den, nom, nom0, betanom, alpha, beta, r0;
  FloatArray r(dim), d(dim), Ad(dim);

  A.Mult(x, r);
  r = b - r; // r = b - A x
  // subtract(b, r, r); // r = b - A x
  d = r;
  nom0 = nom = dot(r, r);

  if (print_iter == 1)
    std::cout << "   Iteration : " << std::setw(3) << 0 << "  (r, r) = " << nom
              << std::endl;

  if ((r0 = nom * RTOLERANCE) < ATOLERANCE)
    r0 = ATOLERANCE;
  if (nom < r0)
    return;

  A.Mult(d, Ad);
  den = d * Ad;

  if (den <= 0.0) {
    if (nom0 > 0.0)
      std::cout << "Operator A is not postive definite. (Ar,r) = " << den
                << std::endl;
    return;
  }

  // start iteration                          //  d = r, Ad = A r
  for (i = 1; i < max_num_iter; i++) {
    alpha = nom / den; // alpha = (r_o,r_o)/(Ar_o,r_o)
    x = x + alpha * d;
    r = r - alpha * Ad; // r_n = r_o - alpha * Ad
    // add(x, alpha, d, x);   //   x =   x + alpha * d
    // add(r, -alpha, Ad, r); // r_n = r_o - alpha * Ad
    betanom = r * r; // betanom = (r_o, r_o)

    if (print_iter == 1)
      std::cout << "   Iteration : " << std::setw(3) << i
                << "  (r, r) = " << betanom << std::endl;

    if (betanom < r0) {
      if (print_iter == 2)
        std::cout << "Number of CG iterations: " << i << std::endl;
      else if (print_iter == 3)
        std::cout << "(r_0, r_0) = " << nom0 << std::endl
                  << "(r_N, r_N) = " << betanom << std::endl
                  << "Number of CG iterations: " << i << std::endl;
      break;
    }

    beta = betanom / nom; // beta = (r_n,r_n)/(r_o,r_o)
    r = r + beta * d;     //    r = r_n + beta * d
    // add(r, beta, d, d);   //    d = r_n + beta * d
    A.Mult(d, Ad); //   Ad = A d
    den = d * Ad;  //  den = (d , A d)
    if (den <= 0.0) {
      if (d * d > 0.0)
        std::cout << "Operator A is not postive definite. (Ad,d) = " << den
                  << std::endl;
    }
    nom = betanom; //  nom = (r_n, r_n)
  }
  if (i == max_num_iter && print_iter >= 0) {
    std::cerr << "CG: No convergence!" << std::endl;
    std::cout << "(r_0, r_0) = " << nom0 << std::endl
              << "(r_N, r_N) = " << betanom << std::endl
              << "Number of CG iterations: " << i << std::endl;
  }
}
