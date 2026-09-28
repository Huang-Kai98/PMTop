// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#include "floatarray.h"
#include "matrix.hpp"
#include <iomanip>
#include <iostream>

//*****************************************************************
// Iterative template routine -- BiCGSTAB
//
// BiCGSTAB solves the unsymmetric linear system Ax = b
// using the Preconditioned BiConjugate Gradient Stabilized method
//
// BiCGSTAB follows the algorithm described on p. 27 of the
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
//*****************************************************************

int BiCGSTAB(const Operator &A, FloatArray &x, const FloatArray &b,
             const Operator &M, int &max_iter, double &tol, double atol,
             int printit) {
  int i, n = A.Size();
  double resid;
  double rho_1, rho_2 = 1.0, alpha = 1.0, beta, omega = 1.0;
  FloatArray p(n), phat(n), s(n), shat(n), t(n), v(n), r(n), rtilde(n);

  A.Mult(x, r); //  r = A * x
  r = b - r;
  // subtract(b, r, r); //  r = b - r
  rtilde = r;

  resid = dot(r, r);
  if (printit)
    std::cout << "   iter " << 0 << ",   (r, r) = " << resid << std::endl;
  tol *= resid;
  tol = (atol > tol) ? atol : tol;

  if (resid <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  for (i = 1; i <= max_iter; i++) {
    rho_1 = dot(rtilde, r);
    // rho_1 = rtilde * r;
    if (rho_1 == 0) {
      tol = resid;
      if (printit)
        std::cout << "   iter " << i << ",   (r, r) = " << resid << std::endl;
      return 2;
    }
    if (i == 1)
      p = r;
    else {
      beta = (rho_1 / rho_2) * (alpha / omega);
      p = p - omega * v;
      p = r + beta * p;
      // add(p, -omega, v, p); //  p = p - omega * v
      // add(r, beta, p, p); //  p = r + beta * p
    }
    M.Mult(p, phat); //  phat = M^{-1} * p
    A.Mult(phat, v); //  v = A * phat
    alpha = rho_1 / dot(rtilde, v);
    // alpha = rho_1 / (rtilde * v);
    s = r - alpha * v;
    // add(r, -alpha, v, s); //  s = r - alpha * v
    resid = dot(s, s);
    if (resid < tol) {
      x = x + alpha * phat;
      // x.Add(alpha, phat); //  x = x + alpha * phat
      tol = resid;
      if (printit)
        std::cout << "   iter " << i << ",   (s, s) = " << resid << std::endl;
      return 0;
    }
    if (printit)
      std::cout << "   iter " << i << ",   (s, s) = " << resid << ";   ";
    M.Mult(s, shat); //  shat = M^{-1} * s
    A.Mult(shat, t); //  t = A * shat
    omega = dot(t, s) / dot(t, t);
    // omega = (t * s) / (t * t);
    x += alpha * phat;
    x += omega * shat;
    r = s - omega * t;
    // x.Add(alpha, phat);   //  x += alpha * phat
    // x.Add(omega, shat);   //  x += omega * shat
    // add(s, -omega, t, r); //  r = s - omega * t

    rho_2 = rho_1;
    resid = dot(r, r);
    // resid = (r * r);
    if (printit)
      std::cout << "(r, r) = " << resid << std::endl;
    if (resid < tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
    if (omega == 0) {
      tol = resid;
      return 3;
    }
  }

  tol = resid;
  return 1;
}
