
#include "NumericalAlgebra.hpp"
#include "pardiso.hpp"
#include "Arpack.hpp"
#include <fstream>
#include <iostream>
#include <mkl.h>
#include <mkl_cblas.h>
#include <mkl_service.h>
#include <mkl_spblas.h>
#include <mkl_types.h>
int main() {

  const int N = 8;   // Matrix size
  const int nev = 3; // Number of eigenvalues

  IntArray ia = { 0, 4, 7, 9, 11, 14, 16, 17, 18};
  IntArray ja =
  { 0,   2,       5, 6,
      1, 2,    4,
         2,             7,
            3,       6,
               4, 5, 6,
                  5,    7,
                     6,
                        7
  };
  FloatArray a =
  { 1.0,      1.0,           2.0, 7.0,
        1.0, 8.0,      2.0,
              1.0,                     5.0,
                   1.0,           9.0,
                        1.0, 1.0, 5.0,
                            1.0,      5.0,
                                 1.0,
                                       1.0
  };
  IntArray iab = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  IntArray jab = {0, 1, 2, 3, 4, 5, 6, 7};
  FloatArray ab = {0.5, 1, 1, 1, 1, 1, 1, 1};
  FloatArray b = {1, 1, 1, 1, 1, 1, 1, 1};
  FloatArray b2 = {1, 1, 1, 1, 1, 1, 1, 1};


  SparseMatrix B(ia.givePointer(), ja.givePointer(), a.givePointer(), N, N,
                 true);
  SparseMatrix A(iab.givePointer(), jab.givePointer(), ab.givePointer(), N, N,
                 true);

    SparseMatrix C(iab.givePointer(), jab.givePointer(), b.givePointer(), N, N,
                 true);

  B.covertToMKLCSR();
    C.covertToMKLCSR();


    FloatArray res(N);
    PardisoSolver M;
    M.SetMatrixType(PardisoSolver::REAL_NONSYMMETRIC);
    M.SetOperator(B);
    M.Mult(b,res);
    res.printYourself();

  FloatArray res2(N);
  A.Mult(res,res2);
  res2.printYourself();

  ArpackSolver M1(N,nev);
    M1.Mult(C);
    FloatArray res3=M1.getEigenvalues();
    FloatArray res4=M1.getEigenvectors();
    res3.printYourself();
    res4.printYourself();


    PardisoSolver MA;
    MA.SetOperator(A);
    M1.SetOption(arpack::which::largest_magnitude);
    M1.Mult(MA,C);

    FloatArray res5=M1.getEigenvalues();
    FloatArray res6=M1.getEigenvectors();
    res5.printYourself();
    res6.printYourself();


  return 0;
}