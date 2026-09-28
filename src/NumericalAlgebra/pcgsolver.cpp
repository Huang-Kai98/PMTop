#include "floatarray.h"
#include "matrix.hpp"
#include "sparsemat.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>

// Preconditined conjugate gradient

void PCG(const Operator &A, const Operator &B, const FloatArray &b,
         FloatArray &x, int print_iter = 0, int max_num_iter = 1000,
         double RTOLERANCE = 10e-12, double ATOLERANCE = 10e-24, int save = 0) {
  int i, dim = x.giveSize();
  double r0, den, nom, nom0, betanom, alpha, beta;
  FloatArray r(dim), d(dim), z(dim);

  A.Mult(x, r); //    r = A x
  r = b - r;    //    r = b  - A x
  // subtract(b, r, r); //    r = b  - r
  B.Mult(r, z); //    z = B r
  d = z;
  nom0 = nom = z * r;

  if (print_iter == 1)
    std::cout << "   Iteration : " << std::setw(3) << 0
              << "  (B r, r) = " << nom << std::endl;

  if ((r0 = nom * RTOLERANCE) < ATOLERANCE)
    r0 = ATOLERANCE;
  if (nom < r0)
    return;

  A.Mult(d, z);
  den = z * d;

  if (den < 0.0) {
    std::cout << "Negative denominator in step 0 of PCG: ";
    std::cout << den << std::endl;
    //    return;
  }

  if (den == 0.0)
    return;

  // start iteration
  for (i = 1; i <= max_num_iter; i++) {

    if (save)
      if (i % save == 0) {
        std::cout << "saving the solution vector on iteration " << i
                  << std::endl;
        std::ofstream out("pcg.x");
        out << x;
      }

    alpha = nom / den;
    x += alpha * d; //  x = x + alpha d
    r -= alpha * z; //  r = r - alpha z
    // add(x, alpha, d, x);  //  x = x + alpha d
    // add(r, -alpha, z, r); //  r = r - alpha z

    B.Mult(r, z); //  z = B r
    betanom = r * z;

    if (print_iter == 1)
      std::cout << "   Iteration : " << std::setw(3) << i
                << "  (B r, r) = " << betanom << std::endl;

    if (betanom < r0) {
      if (print_iter == 2)
        std::cout << "Number of PCG iterations: " << i << std::endl;
      else if (print_iter == 3)
        std::cout << "(B r_0, r_0) = " << nom0 << std::endl
                  << "(B r_N, r_N) = " << betanom << std::endl
                  << "Number of PCG iterations: " << i << std::endl;
      break;
    }

    beta = betanom / nom;
    d = z + beta * d; //  d = z + beta d
    // add(z, beta, d, d); //  d = z + beta d
    A.Mult(d, z);
    den = d * z;
    nom = betanom;
  }
  if (i > max_num_iter) {
    std::cerr << "PCG: No convergence!" << std::endl;
    std::cout << "(B r_0, r_0) = " << nom0 << std::endl
              << "(B r_N, r_N) = " << betanom << std::endl
              << "Number of PCG iterations: " << (i - 1) << std::endl;
  }
  if (print_iter >= 1 || i > max_num_iter) {
    if (i > max_num_iter)
      i--;
    std::cout << "Average reduction factor = " << pow(betanom / nom0, 0.5 / i)
              << std::endl;
  }
}

// Preconditioned stationary linear iteration
void SLI(const Operator &A, const Operator &B, const FloatArray &b,
         FloatArray &x, int print_iter = 0, int max_num_iter = 1000,
         double RTOLERANCE = 10e-12, double ATOLERANCE = 10e-24) {
  int i, dim = x.giveSize();
  double r0, nom, nomold = 1, nom0, cf;
  FloatArray r(dim), z(dim);

  r0 = -1.0;

  for (i = 1; i < max_num_iter; i++) {
    A.Mult(x, r); //    r = A x
    r = b - r;    //    r = b  - A x
    // subtract(b, r, r); //    r = b  - A x
    B.Mult(r, z); //    z = B r

    nom = z * r;

    if (r0 == -1.0) {
      nom0 = nom;
      r0 = nom * RTOLERANCE;
      if (r0 < ATOLERANCE)
        r0 = ATOLERANCE;
    }

    cf = sqrt(nom / nomold);
    if (print_iter == 1) {
      std::cout << "   Iteration : " << std::setw(3) << i
                << "  (B r, r) = " << nom;
      if (i > 1)
        std::cout << "\tConv. rate: " << cf;
      std::cout << std::endl;
    }
    nomold = nom;

    if (nom < r0) {
      if (print_iter == 2)
        std::cout << "Number of iterations: " << i << std::endl
                  << "Conv. rate: " << cf << std::endl;
      else if (print_iter == 3)
        std::cout << "(B r_0, r_0) = " << nom0 << std::endl
                  << "(B r_N, r_N) = " << nom << std::endl
                  << "Number of iterations: " << i << std::endl;
      break;
    }
    x += z; //  x = x + z
    // add(x, 1.0, z, x); //  x = x + B (b - A x)
  }

  if (i == max_num_iter) {
    std::cerr << "No convergence!" << std::endl;
    std::cout << "(B r_0, r_0) = " << nom0 << std::endl
              << "(B r_N, r_N) = " << nom << std::endl
              << "Number of iterations: " << i << std::endl;
  }
}
