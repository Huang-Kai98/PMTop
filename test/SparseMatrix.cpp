// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#include "EigenSolver.hpp"
#include "NumericalAlgebra.hpp"
#include "PardisoSolver.hpp"
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

    IntArray ia = {0, 4, 7, 9, 11, 14, 16, 17, 18};
    IntArray ja = {0, 2, 5, 6, 1, 2, 4, 2, 7, 3, 6, 4, 5, 6, 5, 7, 6, 7};
    FloatArray a = {8.0, 1.0, 2.0, 7.0, -4.0, 8.0,  2.0, 1.0,  5.0,
                    7.0, 9.0, 5.0, 1.0, 5.0,  -1.0, 5.0, 11.0, 5.0};
    IntArray iab = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    IntArray jab = {0, 1, 2, 3, 4, 5, 6, 7};
    FloatArray ab = {1, 1, 1, 1, 1, 1, 1, 1};
    FloatArray b = {1, 1, 1, 1, 1, 1, 1, 1};
    FloatArray b2 = {1, 1, 1, 1, 1, 1, 1, 1};

    SparseMatrix B(ia.givePointer(), ja.givePointer(), a.givePointer(), N, N,
                   true);
    // B.Print();

    SparseMatrix A(N, N);
    A.Add(0, 0, 1);
    A.Add(1, 1, 1);
    A.Add(2, 2, 1);
    A.Add(3, 3, 1);
    A.Add(4, 4, 1);
    A.Add(5, 5, 1);
    A.Add(6, 6, 1);
    A.Add(7, 7, 1);

    // A.Print();
    FloatArray y(N);

    IntArray i = {1, 2};
    IntArray j = {1, 2};
    FloatMatrix subm(2, 2);
    subm(0, 0) = 1;
    subm(0, 1) = 2;
    subm(1, 0) = 3;
    subm(1, 1) = 4;
    // A.AddSubMatrix(i, j, subm);
    A.Finalize();
    A.covertToMKLCSR();
    // A.Print();
    B.covertToMKLCSR();
    // std::cout << B.InnerProduct(b, b2);

    FloatMatrix F(N, N);
    F(0, 0) = 1;
    F(1, 1) = 1;
    F(2, 2) = 1;
    F(3, 3) = 1;
    F(4, 4) = 1;
    F(5, 5) = 1;
    F(6, 6) = 1;
    F(7, 7) = 1;

    FloatMatrix FR = B.RAP(F);
    std::cout << FR << std::endl;
    SparseMatrix *C = Mult(A, B);

    B.Print();

    SparseMatrix *C2 = RAP(B, A);
    C2->Print();

    delete C2;
    // C->Print();
    SparseMatrix *D = Transpose(B);

    // D->Print();

    delete D;

    delete C;

    return 0;
}