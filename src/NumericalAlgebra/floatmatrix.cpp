// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#include "floatmatrix.h"
#include "error.hpp"
#include "floatarray.h"
#include "mathfem.h"
#include "operator.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mkl_service.h>
#include <mkl_types.h>
#include <numeric>
#include <ostream>
#ifdef __MKL_MODULE
#include "mkl.h"
#endif

#define RESIZE(nr, nc)                                                         \
  {                                                                            \
    this->nRows = nr;                                                          \
    this->nColumns = nc;                                                       \
    int nsize = this->nRows * this->nColumns;                                  \
    if (nsize < (int)this->values.size()) {                                    \
      this->values.resize(nsize);                                              \
    } else if (nsize > (int)this->values.size()) {                             \
      this->values.assign(nsize, 0.);                                          \
    }                                                                          \
  }

FloatMatrix ::FloatMatrix(const FloatArray &vector, bool transpose)
//
// constructor : creates (vector->giveSize(),1) FloatMatrix
// if transpose = 1 creates (1,vector->giveSize()) FloatMatrix
//
{
  if (transpose) {
    nRows = 1; // column vector
    nColumns = vector.giveSize();
  } else {
    nRows = vector.giveSize(); // row vector- default
    nColumns = 1;
  }

  values = vector.values;
}

FloatMatrix FloatMatrix::transpose() const {
  FloatMatrix result(nColumns, nRows);
  for (int i = 0; i < nRows; ++i) {
    for (int j = 0; j < nColumns; ++j) {
      result(j, i) = (*this)(i, j);
    }
  }
  return result;
}
void FloatMatrix::GetColumnReference(int c, FloatArray &col) {
  col.SetDataAndSize(values.data() + c * nRows, nRows);
}

void FloatMatrix::GetColumn(int c, FloatArray &col) {
  col.resize(nRows);
  for (int i = 0; i < nRows; i++) {
    col[i] = (*this)(i, c);
  }
}

FloatMatrix::FloatMatrix(const FloatMatrix &mat, char ch) {
  if (ch == 't') {
    (*this) = mat.transpose();
  } else {
    (*this) = mat;
  }
}
void CalcAdjugate(const FloatMatrix &a, FloatMatrix &adja) {
#ifdef MFEM_DEBUG
  if (a.isSquare() || adjat.isSquare() || a.giveNumberOfRows() < 2 ||
      a.giveNumberOfRows() > 3) {
    pmtop_error("FloatMatrix::CalcAdjugateTranspose (...)");
  }
#endif
  if (a.giveNumberOfRows() == 2) {
    adja(0, 0) = a(1, 1);
    adja(0, 1) = -a(0, 1);
    adja(1, 0) = -a(1, 0);
    adja(1, 1) = a(0, 0);
  } else {
    adja(0, 0) = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1);
    adja(0, 1) = a(0, 2) * a(2, 1) - a(0, 1) * a(2, 2);
    adja(0, 2) = a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1);

    adja(1, 0) = a(1, 2) * a(2, 0) - a(1, 0) * a(2, 2);
    adja(1, 1) = a(0, 0) * a(2, 2) - a(0, 2) * a(2, 0);
    adja(1, 2) = a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2);

    adja(2, 0) = a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0);
    adja(2, 1) = a(0, 1) * a(2, 0) - a(0, 0) * a(2, 1);
    adja(2, 2) = a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0);
  }
};

void CalcAdjugateTranspose(const FloatMatrix &a, FloatMatrix &adjat) {

#ifdef NDEBUG
  if (a.isSquare() || adjat.isSquare() || a.giveNumberOfRows() < 2 ||
      a.giveNumberOfRows() > 3) {
    pmtop_error("FloatMatrix::CalcAdjugateTranspose (...)");
  }
#endif
  if (a.giveNumberOfRows() == 2) {
    adjat(0, 0) = a(1, 1);
    adjat(1, 0) = -a(0, 1);
    adjat(0, 1) = -a(1, 0);
    adjat(1, 1) = a(0, 0);
  } else {
    adjat(0, 0) = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1);
    adjat(1, 0) = a(0, 2) * a(2, 1) - a(0, 1) * a(2, 2);
    adjat(2, 0) = a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1);

    adjat(0, 1) = a(1, 2) * a(2, 0) - a(1, 0) * a(2, 2);
    adjat(1, 1) = a(0, 0) * a(2, 2) - a(0, 2) * a(2, 0);
    adjat(2, 1) = a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2);

    adjat(0, 2) = a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0);
    adjat(1, 2) = a(0, 1) * a(2, 0) - a(0, 0) * a(2, 1);
    adjat(2, 2) = a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0);
  }
}

void FloatMatrix::GradToDiv(FloatArray &div) {

#ifdef NDEBUG
  if (nColumns * nRows != div.giveSize())
    pmtop_error("DenseMatrix::GradToDiv (...)");
#endif

  // div(dof*j+i) <-- (*this)(i,j)

  int n = nRows * nColumns;
  double *ddata = div.givePointer();

  for (int i = 0; i < n; i++)
    ddata[i] = values[i];
}

FloatMatrix ::FloatMatrix(
    std ::initializer_list<std ::initializer_list<double>> mat) {
  RESIZE(mat.size(), mat.begin()->size())
  auto p = this->values.begin();
  for (auto col : mat) {
#ifndef NDEBUG
    if (this->nRows != (int)col.size()) {
      std::cerr << "Initializer list has inconsistent column sizes."
                << std::endl;
    }
#endif
    for (auto x : col) {
      *p = x;
      p++;
    }
  }
}

FloatMatrix &FloatMatrix ::operator=(
    std ::initializer_list<std ::initializer_list<double>> mat) {
  RESIZE((int)mat.begin()->size(), (int)mat.size());
  auto p = this->values.begin();
  for (auto col : mat) {
#ifndef NDEBUG
    if (this->nRows != (int)col.size()) {
      std::cerr << "Initializer list has inconsistent column sizes."
                << std::endl;
    }
#endif
    for (auto x : col) {
      *p = x;
      p++;
    }
  }

  return *this;
}
FloatMatrix &FloatMatrix ::operator=(std ::initializer_list<FloatArray> mat) {
  RESIZE(mat.begin()->giveSize(), (int)mat.size());
  auto p = this->values.begin();
  for (auto col : mat) {
#ifndef NDEBUG
    if (this->nRows != col.giveSize()) {
      std::cerr << "Initializer list has inconsistent column sizes."
                << std::endl;
    }
#endif
    for (auto x : col) {
      *p = x;
      p++;
    }
  }

  return *this;
}

void FloatMatrix ::checkBounds(int i, int j) const
// Checks that the receiver includes a position (i,j).
{
  if (i <= 0) {
    std::cerr << "matrix error on rows : " << i << " <= 0" << std::endl;
  }
  if (j <= 0) {
    std::cerr << "matrix error on columns : " << j << " <= 0" << std::endl;
  }

  if (i > nRows) {
    std::cerr << "matrix error on rows : " << i << " > " << nRows << std::endl;
  }

  if (j > nColumns) {
    printf("APA \n");
    std::cerr << "matrix error on columns : " << j << " > " << nColumns
              << std::endl;
  }
}

bool FloatMatrix ::isFinite() const {
  for (double val : values) {
    if (!std::isfinite(val)) {
      return false;
    }
  }

  return true;
}

FloatMatrix FloatMatrix::HadamardProduct(FloatArray &x) const {
  FloatMatrix Result(this->giveNumberOfRows(),this->giveNumberOfColumns());
  #ifndef NDEBUG
  if (x.giveSize()!=this->nRows) {
    std::cerr << "dimensions of 'x' and 'matrix' mismatch" << std::endl;
  }
#endif
  for (int i = 0; i < x.giveSize(); i++) {
    for (int j = 0; j < this->giveNumberOfColumns(); j++) {
      Result(i,j) = (*this)(i,j)*x(i);
    }
  }
  return Result;
}


void FloatMatrix ::assemble(const FloatMatrix &src, const IntArray &loc) {
  int ii, jj, size = src.giveNumberOfRows();

#ifndef NDEBUG
  if (size != loc.giveSize()) {
    std::cerr << "dimensions of 'src' and 'loc' mismatch" << std::endl;
  }

  if (!src.isSquare()) {
    std::cerr << "'src' is not sqaure matrix" << std::endl;
  }
#endif

  for (int i = 1; i <= size; i++) {
    if ((ii = loc.at(i))) {
      for (int j = 1; j <= size; j++) {
        if ((jj = loc.at(j))) {
          this->at(ii, jj) += src.at(i, j);
        }
      }
    }
  }
}

void FloatMatrix ::assemble(const FloatMatrix &src, const IntArray &rowind,
                            const IntArray &colind) {
  int ii, jj;
  int nr = src.giveNumberOfRows();
  int nc = src.giveNumberOfColumns();

#ifndef NDEBUG
  if (nr != rowind.giveSize()) {
    std::cerr << "row dimensions of 'src' and 'rowind' mismatch" << std::endl;
  }

  if (nc != colind.giveSize()) {
    std::cerr << "column dimensions of 'src' and 'colind' mismatch"
              << std::endl;
  }
#endif

  for (int i = 1; i <= nr; i++) {
    if ((ii = rowind.at(i))) {
      for (int j = 1; j <= nc; j++) {
        if ((jj = colind.at(j))) {
          this->at(ii, jj) += src.at(i, j);
        }
      }
    }
  }
}

void FloatMatrix ::assembleT(const FloatMatrix &src, const IntArray &rowind,
                             const IntArray &colind) {
  int ii, jj;
  int nr = src.giveNumberOfRows();
  int nc = src.giveNumberOfColumns();

#ifndef NDEBUG
  if (nr != rowind.giveSize()) {
    std::cerr << "row dimensions of 'src' and 'rowind' mismatch" << std::endl;
  }

  if (nc != colind.giveSize()) {
    std::cerr << "column dimensions of 'src' and 'colind' mismatch"
              << std::endl;
  }
#endif

  for (int i = 1; i <= nr; i++) {
    if ((ii = rowind.at(i))) {
      for (int j = 1; j <= nc; j++) {
        if ((jj = colind.at(j))) {
          this->at(jj, ii) += src.at(i, j);
        }
      }
    }
  }
}

void FloatMatrix ::assemble(const FloatMatrix &src, const int *rowind,
                            const int *colind) {
  int ii, jj;
  int nr = src.giveNumberOfRows();
  int nc = src.giveNumberOfColumns();

  for (int i = 1; i <= nr; i++) {
    if ((ii = rowind[i - 1])) {
      for (int j = 1; j <= nc; j++) {
        if ((jj = colind[j - 1])) {
          this->at(ii, jj) += src.at(i, j);
        }
      }
    }
  }
}
void FloatMatrix ::beTranspositionOf(const FloatMatrix &src) {
  // receiver becomes a transposition of src
  int nrows = src.giveNumberOfColumns(), ncols = src.giveNumberOfRows();
  RESIZE(nrows, ncols);

  for (int i = 1; i <= nrows; i++) {
    for (int j = 1; j <= ncols; j++) {
      this->at(i, j) = src.at(j, i);
    }
  }
}

void FloatMatrix::Lump() {
  for (int i = 0; i < nRows; i++) {
    double L = 0.0;
    for (int j = 0; j < nColumns; j++) {
      L += (*this)(i, j);
      (*this)(i, j) = 0.0;
    }
    (*this)(i, i) = L;
  }
}
void FloatMatrix ::beProductOf(const FloatMatrix &aMatrix,
                               const FloatMatrix &bMatrix)
// Receiver = aMatrix * bMatrix
{
#ifndef NDEBUG
  if (aMatrix.nColumns != bMatrix.nRows) {
    std::cerr << "error in product A*B : dimensions do not match, A(*,"
              << aMatrix.nColumns << "), B(" << bMatrix.nRows << ",*)"
              << std::endl;
  }
#endif
  RESIZE(aMatrix.nRows, bMatrix.nColumns);
#ifdef __MKL_MODULE
  double alpha = 1., beta = 0.;
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, this->nRows,
              this->nColumns, aMatrix.nColumns, alpha, aMatrix.givePointer(),
              aMatrix.nRows, bMatrix.givePointer(), bMatrix.nRows, beta,
              this->givePointer(), this->nRows);
#else
  for (int i = 1; i <= aMatrix.nRows; i++) {
    for (int j = 1; j <= bMatrix.nColumns; j++) {
      double coeff = 0.;
      for (int k = 1; k <= aMatrix.nColumns; k++) {
        coeff += aMatrix.at(i, k) * bMatrix.at(k, j);
      }

      this->at(i, j) = coeff;
    }
  }
#endif
}

void FloatMatrix ::beTProductOf(const FloatMatrix &aMatrix,
                                const FloatMatrix &bMatrix)
// Receiver = aMatrix^T * bMatrix
{
#ifndef NDEBUG
  if (aMatrix.nRows != bMatrix.nRows) {
    std::cerr << "error in product A*B : dimensions do not match" << std::endl;
  }
#endif
  RESIZE(aMatrix.nColumns, bMatrix.nColumns);
#ifdef __MKL_MODULE
  double alpha = 1., beta = 0.;
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, this->nRows,
              this->nColumns, aMatrix.nRows, alpha, aMatrix.givePointer(),
              aMatrix.nRows, bMatrix.givePointer(), bMatrix.nRows, beta,
              this->givePointer(), this->nRows);
#else
  for (int i = 1; i <= aMatrix.nColumns; i++) {
    for (int j = 1; j <= bMatrix.nColumns; j++) {
      double coeff = 0.;
      for (int k = 1; k <= aMatrix.nRows; k++) {
        coeff += aMatrix.at(k, i) * bMatrix.at(k, j);
      }

      this->at(i, j) = coeff;
    }
  }
#endif
}
void FloatMatrix ::beProductTOf(const FloatMatrix &aMatrix,
                                const FloatMatrix &bMatrix)
// Receiver = aMatrix * bMatrix^T
{
#ifndef NDEBUG
  if (aMatrix.nColumns != bMatrix.nColumns) {
    std::cerr << "error in product A*B : dimensions do not match" << std::endl;
  }
#endif
  RESIZE(aMatrix.nRows, bMatrix.nRows);
#ifdef __MKL_MODULE
  double alpha = 1., beta = 0.;
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, this->nRows,
              this->nColumns, aMatrix.nColumns, alpha, aMatrix.givePointer(),
              aMatrix.nRows, bMatrix.givePointer(), bMatrix.nRows, beta,
              this->givePointer(), this->nRows);
#else
  for (int i = 1; i <= aMatrix.nRows; i++) {
    for (int j = 1; j <= bMatrix.nRows; j++) {
      double coeff = 0.;
      for (int k = 1; k <= aMatrix.nColumns; k++) {
        coeff += aMatrix.at(i, k) * bMatrix.at(j, k);
      }

      this->at(i, j) = coeff;
    }
  }
#endif
}

FloatArray FloatMatrix::ComputeEigenValueAndEigenVector() {
  FloatArray EigenValue(this->nRows);
  LAPACKE_dsyevd(CblasColMajor, 'V', 'U', this->nRows, this->givePointer(), this->nRows,
                 EigenValue.givePointer());
  return EigenValue;
}

void FloatMatrix ::addProductOf(const FloatMatrix &aMatrix,
                                const FloatMatrix &bMatrix)
// Receiver = aMatrix * bMatrix
{
#ifndef NDEBUG
  if (aMatrix.nColumns != bMatrix.nRows) {
    std::cerr << "error in product A*B : dimensions do not match" << std::endl;
  }
  if (aMatrix.nRows != this->nRows || bMatrix.nColumns != this->nColumns) {
    std::cerr << "error in product receiver : dimensions do not match"
              << std::endl;
  }
#endif

#ifdef __MKL_MODULE
  double alpha = 1., beta = 1.;
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, this->nRows,
              this->nColumns, aMatrix.nColumns, alpha, aMatrix.givePointer(),
              aMatrix.nRows, bMatrix.givePointer(), bMatrix.nRows, beta,
              this->givePointer(), this->nRows);
#else
  for (int i = 1; i <= aMatrix.nRows; i++) {
    for (int j = 1; j <= bMatrix.nColumns; j++) {
      double coeff = 0.;
      for (int k = 1; k <= aMatrix.nColumns; k++) {
        coeff += aMatrix.at(i, k) * bMatrix.at(k, j);
      }

      this->at(i, j) += coeff;
    }
  }
#endif
}
void FloatMatrix ::addTProductOf(const FloatMatrix &aMatrix,
                                 const FloatMatrix &bMatrix)
// Receiver += aMatrix^T * bMatrix
{
#ifndef NDEBUG
  if (aMatrix.nRows != bMatrix.nRows) {
    std::cerr << "error in product A*B : dimensions do not match" << std::endl;
  }
  if (aMatrix.nColumns != this->nColumns || bMatrix.nColumns != this->nRows) {
    std::cerr << "error in product receiver : dimensions do not match"
              << std::endl;
  }
#endif

#ifdef __MKL_MODULE
  double alpha = 1., beta = 1.;
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, this->nRows,
              this->nColumns, aMatrix.nRows, alpha, aMatrix.givePointer(),
              aMatrix.nRows, bMatrix.givePointer(), bMatrix.nRows, beta,
              this->givePointer(), this->nRows);
#else
  for (int i = 1; i <= aMatrix.nColumns; i++) {
    for (int j = 1; j <= bMatrix.nColumns; j++) {
      double coeff = 0.;
      for (int k = 1; k <= aMatrix.nRows; k++) {
        coeff += aMatrix.at(k, i) * bMatrix.at(k, j);
      }

      this->at(i, j) += coeff;
    }
  }
#endif
}
void FloatMatrix ::beDyadicProductOf(const FloatArray &vec1,
                                     const FloatArray &vec2)
// Receiver = vec1 * vec2^T
{
  int n1 = vec1.giveSize();
  int n2 = vec2.giveSize();
  RESIZE(n1, n2);
  for (int j = 1; j <= n2; j++) {
    for (int i = 1; i <= n1; i++) {
      this->at(i, j) = vec1.at(i) * vec2.at(j);
    }
  }
}
void FloatMatrix ::beNMatrixOf(const FloatArray &n, int nsd) {
  this->resize(nsd, n.giveSize() * nsd);
  for (int i = 0; i < n.giveSize(); ++i) {
    for (int j = 0; j < nsd; ++j) {
      (*this)(j, i *nsd + j) = n(i);
    }
  }
}
void FloatMatrix ::beLocalCoordSys(
    const FloatArray &normal) { // normal should be at the first position,
                                // easier for interface material models
  if (normal.giveSize() == 1) {
    this->resize(1, 1);
    this->at(1, 1) = normal(0);
  } else if (normal.giveSize() == 2) {
    this->resize(2, 2);
    this->at(1, 1) = normal(0);
    this->at(1, 2) = normal(1);

    this->at(2, 1) = normal(1);
    this->at(2, 2) = -normal(0);

  } else if (normal.giveSize() == 3) {
    // Create a permutated vector of n, *always* length 1 and significantly
    // different from n.
    FloatArray b,
        t = {normal(1), -normal(2), normal(0)}; // binormal and tangent

    // Construct orthogonal vector
    double npn = t.dotProduct(normal);
    t.add(-npn, normal);
    t.normalize();
    b.beVectorProductOf(t, normal);

    this->resize(3, 3);
    this->at(1, 1) = normal.at(1);
    this->at(1, 2) = normal.at(2);
    this->at(1, 3) = normal.at(3);

    this->at(2, 1) = b.at(1);
    this->at(2, 2) = b.at(2);
    this->at(2, 3) = b.at(3);

    this->at(3, 1) = t.at(1);
    this->at(3, 2) = t.at(2);
    this->at(3, 3) = t.at(3);
  } else {
    std::cerr << "Normal needs 1 to 3 components." << std::endl;
  }
}

void FloatMatrix::setSubMatrix(const FloatMatrix &src, IntArray indexRow, IntArray indexCol) {
  int srcRows = src.giveNumberOfRows(), srcCols = src.giveNumberOfColumns();
  for (int j = 0; j < srcCols; j++) {
    for (int i = 0; i < srcRows; i++) {
      (*this)(indexRow(i), indexCol(j)) = src(i, j);
    }
  }
}

void FloatMatrix ::setSubMatrix(const FloatMatrix &src, int sr, int sc) {
  sr--;
  sc--;

  int srcRows = src.giveNumberOfRows(), srcCols = src.giveNumberOfColumns();
#ifndef NDEBUG
  int nr = sr + srcRows;
  int nc = sc + srcCols;

  if ((this->giveNumberOfRows() < nr) || (this->giveNumberOfColumns() < nc)) {
    std::cerr << "Sub matrix doesn't fit inside allocated space." << std::endl;
  }
#endif

  // add sub-matrix
  for (int j = 0; j < srcCols; j++) {
    for (int i = 0; i < srcRows; i++) {
      (*this)(sr + i, sc + j) = src(i, j);
    }
  }
  // memcpy( &(*this)(sr + 0, sc + j), &src(0,j), srcRows * sizeof(int));
}
void FloatMatrix ::setTSubMatrix(const FloatMatrix &src, int sr, int sc) {
  sr--;
  sc--;

  int srcRows = src.giveNumberOfRows(), srcCols = src.giveNumberOfColumns();
#ifndef NDEBUG
  int nr = sr + srcCols;
  int nc = sc + srcRows;

  if ((this->giveNumberOfRows() < nr) || (this->giveNumberOfColumns() < nc)) {
    std::cerr << "Sub matrix doesn't fit inside allocated space." << std::endl;
  }
#endif

  // add sub-matrix
  for (int i = 0; i < srcCols; i++) {
    for (int j = 0; j < srcRows; j++) {
      (*this)(sr + i, sc + j) = src(j, i);
    }
  }
}
void FloatMatrix ::addSubVectorRow(const FloatArray &src, int sr, int sc) {
  sc--;

  int srcCols = src.giveSize();

  int nr = sr;
  int nc = sc + srcCols;

  if ((this->giveNumberOfRows() < nr) || (this->giveNumberOfColumns() < nc)) {
    this->resizeWithData(PMTop::max(this->giveNumberOfRows(), nr),
                         PMTop::max(this->giveNumberOfColumns(), nc));
  }

  // add sub-matrix
  for (int j = 1; j <= srcCols; j++) {
    this->at(sr, sc + j) += src.at(j);
  }
}
void FloatMatrix ::addSubVectorCol(const FloatArray &src, int sr, int sc) {
  sr--;

  int srcRows = src.giveSize();

  int nr = sr + srcRows;
  int nc = sc;

  if ((this->giveNumberOfRows() < nr) || (this->giveNumberOfColumns() < nc)) {
    this->resizeWithData(PMTop::max(this->giveNumberOfRows(), nr),
                         PMTop::max(this->giveNumberOfColumns(), nc));
  }

  // add sub-matrix
  for (int j = 1; j <= srcRows; j++) {
    this->at(sr + j, sc) += src.at(j);
  }
}
void FloatMatrix ::setColumn(const FloatArray &src, int c) {
  int nr = src.giveSize();
#ifndef NDEBUG
  if (this->giveNumberOfRows() != nr || c < 1 ||
      c > this->giveNumberOfColumns()) {
    std::cerr << "Size mismatch" << std::endl;
  }
#endif


#ifdef __MKL_MODULE
  cblas_dcopy(nr, src.givePointer(), 1, this->givePointer()+(c - 1) * nr, 1);
#else
  auto P = this->values.begin() + (c - 1) * nr;
  std ::copy(src.begin(), src.end(), P);
#endif
}

void FloatMatrix ::copyColumn(FloatArray &dest, int c) const {
  int nr = this->giveNumberOfRows();
#ifndef NDEBUG
  if (c < 1 || c > this->giveNumberOfColumns()) {
    std::cerr << "Column outside range (" << c << ")" << std::endl;
  }
#endif

  dest.resize(nr);
  auto P = this->values.begin() + (c - 1) * nr;
  std ::copy(P, P + nr, dest.begin());
}
void FloatMatrix ::copySubVectorRow(const FloatArray &src, int sr, int sc) {
  sc--;

  int srcCols = src.giveSize();

  int nr = sr;
  int nc = sc + srcCols;

  if ((this->giveNumberOfRows() < nr) || (this->giveNumberOfColumns() < nc)) {
    this->resizeWithData(PMTop::max(this->giveNumberOfRows(), nr),
                         PMTop::max(this->giveNumberOfColumns(), nc));
  }

  // add sub-matrix
  for (int j = 1; j <= srcCols; j++) {
    this->at(sr, sc + j) = src.at(j);
  }
}
void FloatMatrix ::plusProductSymmUpper(const FloatMatrix &a,
                                        const FloatMatrix &b, double dV)
// Adds to the receiver the product  a(transposed).b dV .
// The receiver size is adjusted, if necessary.
// This method assumes that both the receiver and the product above are
// symmetric matrices, and therefore computes only the upper half of the
// receiver ; the lower half is not modified. Other advantage : it does
// not compute the transposition of matrix a.
{
  if (!this->isNotEmpty()) {
    this->nRows = a.nColumns;
    this->nColumns = b.nColumns;
    this->values.assign(a.nColumns * b.nColumns, 0.);
  }

#ifdef __MKL_MODULE
  ///@todo We should determine which is the best choice overall. For large
  /// systems more block matrix operations is necessary.
  /// For smaller systems the overhead from function calls might be larger, but
  /// the overhead might be tiny, or using symmetry at all might be
  /// undesireable.
  double beta = 1.;
  if (this->nRows < 20) {
    // Split the matrix into 2 columns, s1 + s2 = n ( = nRows = nColumns ).
    MKL_INT s1 = this->nRows / 2;
    MKL_INT s2 = this->nRows - s1;
    // First column block, we only take the first s rows by only taking the
    // first s columns in the matrix a.
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, s1, s1, a.nColumns, dV,
                a.givePointer(), a.nRows, b.givePointer(), b.nRows, beta,
                this->givePointer(), this->nRows);
    // Second column block starting a memory position c * nRows
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, this->nRows, s2,
                a.nRows, dV, a.givePointer(), a.nRows,
                &b.givePointer()[s1 * b.nRows], b.nRows, beta,
                &this->givePointer()[s1 * this->nRows], this->nRows);
  } else {
    // Get suitable blocksize. Around 10 rows should be suitable (slightly
    // adjusted to minimize number of blocks):
    MKL_INT block = (this->nRows - 1) / (this->nRows / 10) + 1;
    MKL_INT start = 0;
    MKL_INT end = block;
    while (start < this->nRows) {
      MKL_INT s = end - start;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, end, s, a.nRows, dV,
                  a.givePointer(), a.nRows, &b.givePointer()[start + b.nRows],
                  b.nRows, beta, &this->givePointer()[start * this->nRows],
                  this->nRows);
      start = end;
      end += block;
      if (end > this->nRows) {
        end = this->nRows;
      }
    }
  }
#else
  for (int i = 1; i <= nRows; i++) {
    for (int j = i; j <= nColumns; j++) {
      double summ = 0.;
      for (int k = 1; k <= a.nRows; k++) {
        summ += a.at(k, i) * b.at(k, j);
      }

      this->at(i, j) += summ * dV;
    }
  }
#endif
}

void FloatMatrix ::plusDyadSymmUpper(const FloatArray &a, double dV) {
  if (!this->isNotEmpty()) {
    this->nRows = a.giveSize();
    this->nColumns = a.giveSize();
    this->values.assign(this->nRows * this->nColumns, 0.);
  }
#ifdef __MKL_MODULE
  MKL_INT inc = 1;
  MKL_INT sizeA = a.giveSize();
  cblas_dsyr(CblasColMajor, CblasUpper, sizeA, dV, a.givePointer(), inc,
             this->givePointer(), sizeA);
#else
  for (int i = 1; i <= nRows; i++) {
    for (int j = i; j <= nColumns; j++) {
      this->at(i, j) += a.at(i) * a.at(j) * dV;
    }
  }
#endif
}

void FloatMatrix ::plusProductUnsym(const FloatMatrix &a, const FloatMatrix &b,
                                    double dV)
// Adds to the receiver the product  a(transposed).b dV .
// If the receiver has a null size, it is expanded.
// Advantage : does not compute the transposition of matrix a.
{
  if (!this->isNotEmpty()) {
    this->nRows = a.nColumns;
    this->nColumns = b.nColumns;
    this->values.assign(this->nRows * this->nColumns, 0.);
  }
#ifdef __MKL_MODULE
  double beta = 1.;
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, this->nRows,
              this->nColumns, a.nRows, dV, a.givePointer(), a.nRows,
              b.givePointer(), b.nRows, beta, this->givePointer(), this->nRows);
#else
  for (int i = 1; i <= nRows; i++) {
    for (int j = 1; j <= nColumns; j++) {
      double summ = 0.;
      for (int k = 1; k <= a.nRows; k++) {
        summ += a.at(k, i) * b.at(k, j);
      }

      this->at(i, j) += summ * dV;
    }
  }
#endif
}
void FloatMatrix ::plusDyadUnsym(const FloatArray &a, const FloatArray &b,
                                 double dV) {
  if (!this->isNotEmpty()) {
    this->nRows = a.giveSize();
    this->nColumns = b.giveSize();
    this->values.assign(this->nRows * this->nColumns, 0.);
  }
#ifdef __MKL_MODULE
  MKL_INT inc = 1;
  MKL_INT sizeA = a.giveSize();
  MKL_INT sizeB = b.giveSize();
  cblas_dger(CblasColMajor, sizeA, sizeB, dV, a.givePointer(), inc,
             b.givePointer(), inc, this->givePointer(), this->nRows);
#else
  for (int i = 1; i <= nRows; i++) {
    for (int j = 1; j <= nColumns; j++) {
      this->at(i, j) += a.at(i) * b.at(j) * dV;
    }
  }
#endif
}

bool FloatMatrix ::beInverseOf(const FloatMatrix &src)
// Receiver becomes inverse of given parameter src. If necessary, size is
// adjusted.
{
  double det;

#ifndef NDEBUG
  if (!src.isSquare()) {
    std::cerr << "cannot inverse a " << src.nRows << " by " << src.nColumns
              << " matrix" << std::endl;
  }
#endif

  RESIZE(src.nRows, src.nColumns);

  if (nRows == 1) {
    if (fabs(src.at(1, 1)) > 1.e-30) {
      this->at(1, 1) = 1. / src.at(1, 1);
      return true;
    } else {
      return false;
    }
  } else if (nRows == 2) {
    det = src.at(1, 1) * src.at(2, 2) - src.at(1, 2) * src.at(2, 1);
    if (fabs(det) > 1.e-30) {
      this->at(1, 1) = src.at(2, 2) / det;
      this->at(2, 1) = -src.at(2, 1) / det;
      this->at(1, 2) = -src.at(1, 2) / det;
      this->at(2, 2) = src.at(1, 1) / det;
      return true;
    } else {
      return false;
    }
  } else if (nRows == 3) {
    det = src.at(1, 1) * src.at(2, 2) * src.at(3, 3) +
          src.at(1, 2) * src.at(2, 3) * src.at(3, 1) +
          src.at(1, 3) * src.at(2, 1) * src.at(3, 2) -
          src.at(1, 3) * src.at(2, 2) * src.at(3, 1) -
          src.at(2, 3) * src.at(3, 2) * src.at(1, 1) -
          src.at(3, 3) * src.at(1, 2) * src.at(2, 1);
    if (fabs(det) > 1.e-30) {
      this->at(1, 1) =
          (src.at(2, 2) * src.at(3, 3) - src.at(2, 3) * src.at(3, 2)) / det;
      this->at(2, 1) =
          (src.at(2, 3) * src.at(3, 1) - src.at(2, 1) * src.at(3, 3)) / det;
      this->at(3, 1) =
          (src.at(2, 1) * src.at(3, 2) - src.at(2, 2) * src.at(3, 1)) / det;
      this->at(1, 2) =
          (src.at(1, 3) * src.at(3, 2) - src.at(1, 2) * src.at(3, 3)) / det;
      this->at(2, 2) =
          (src.at(1, 1) * src.at(3, 3) - src.at(1, 3) * src.at(3, 1)) / det;
      this->at(3, 2) =
          (src.at(1, 2) * src.at(3, 1) - src.at(1, 1) * src.at(3, 2)) / det;
      this->at(1, 3) =
          (src.at(1, 2) * src.at(2, 3) - src.at(1, 3) * src.at(2, 2)) / det;
      this->at(2, 3) =
          (src.at(1, 3) * src.at(2, 1) - src.at(1, 1) * src.at(2, 3)) / det;
      this->at(3, 3) =
          (src.at(1, 1) * src.at(2, 2) - src.at(1, 2) * src.at(2, 1)) / det;
      return true;
    } else {
      return false;
    }
  } else {
#ifdef __MKL_MODULE
    MKL_INT n = this->nRows;

    MKL_INT *ipiv = (MKL_INT *)mkl_malloc(n * sizeof(MKL_INT), 64);

    MKL_INT lwork, info;
    *this = src;
    info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, this->givePointer(), n, ipiv);
    if (info != 0) {
      std::cerr << "dgetrf error " << info << std::endl;
      return false;
    }
    lwork = n * n;
    info = LAPACKE_dgetri(LAPACK_COL_MAJOR, n, this->givePointer(), n, ipiv);
    if (info > 0) {
      std::cerr << "Singular at " << info << std::endl;
      return false;
    } else if (info < 0) {
      std::cerr << "Error on input " << info << std::endl;
    }
    mkl_free(ipiv);
    return true;
#else
    // size >3 ... gaussian elimination - slow but safe
    //
    double piv, linkomb;
    FloatMatrix tmp = src;
    // initialize answer to be unity matrix;
    this->zero();
    for (int i = 1; i <= nRows; i++) {
      this->at(i, i) = 1.0;
    }

    // lower triangle elimination by columns
    for (int i = 1; i < nRows; i++) {
      piv = tmp.at(i, i);
      if (fabs(piv) < 1.e-30) {
        std::cerr << "pivot (" << i << "," << i
                  << ") to close to small (< 1.e-20)" << std::endl;
        return false;
      }

      for (int j = i + 1; j <= nRows; j++) {
        linkomb = tmp.at(j, i) / tmp.at(i, i);
        for (int k = i; k <= nRows; k++) {
          tmp.at(j, k) -= tmp.at(i, k) * linkomb;
        }

        for (int k = 1; k <= nRows; k++) {
          this->at(j, k) -= this->at(i, k) * linkomb;
        }
      }
    }

    // upper triangle elimination by columns
    for (int i = nRows; i > 1; i--) {
      piv = tmp.at(i, i);
      for (int j = i - 1; j > 0; j--) {
        linkomb = tmp.at(j, i) / piv;
        for (int k = i; k > 0; k--) {
          tmp.at(j, k) -= tmp.at(i, k) * linkomb;
        }

        for (int k = nRows; k > 0; k--) {
          // tmp -> at(j,k)-= tmp  ->at(i,k)*linkomb;
          this->at(j, k) -= this->at(i, k) * linkomb;
        }
      }
    }

    // diagonal scaling
    for (int i = 1; i <= nRows; i++) {
      for (int j = 1; j <= nRows; j++) {
        this->at(i, j) /= tmp.at(i, i);
      }
    }
    return true;
#endif
  }
}

void FloatMatrix ::beSubMatrixOf(const FloatMatrix &src, int topRow,
                                 int bottomRow, int topCol, int bottomCol)
/*
 * modifies receiver to be  submatrix of the src matrix
 * size of receiver  submatrix is determined from
 * input parameters
 */
{
#ifndef NDEBUG
  if ((topRow < 1) || (bottomRow < 1) || (topCol < 1) || (bottomCol < 1)) {
    std::cerr << "subindexes size mismatch" << std::endl;
  }

  if ((src.nRows < bottomRow) || (src.nColumns < bottomCol) ||
      ((bottomRow - topRow) > src.nRows) ||
      ((bottomCol - topCol) > src.nColumns)) {
    std::cerr << "subindexes size mismatch" << std::endl;
  }
#endif

  int topRm1, topCm1;
  topRm1 = topRow - 1;
  topCm1 = topCol - 1;

  // allocate return value
  this->resize(bottomRow - topRm1, bottomCol - topCm1);
  for (int i = topRow; i <= bottomRow; i++) {
    for (int j = topCol; j <= bottomCol; j++) {
      this->at(i - topRm1, j - topCm1) = src.at(i, j);
    }
  }
}

void FloatMatrix::AddMatrix(double a, FloatMatrix &A, int ro, int co) {
  int h, ah, aw;
  double *p, *ap;

  h = nRows;
  ah = A.giveNumberOfRows();
  aw = A.giveNumberOfColumns();

#ifdef MFEM_DEBUG
  if (co + aw > nColumns || ro + ah > h)
    mfem_error("DenseMatrix::AddMatrix (...) 2");
#endif

  p = values.data() + ro + co * h;
  ap = A.givePointer();

  for (int c = 0; c < aw; c++) {
    for (int r = 0; r < ah; r++)
      p[r] += a * ap[r];
    p += h;
    ap += ah;
  }
}

void FloatMatrix::AddMatrix(FloatMatrix &aMatrix, int ro, int co) {
  int h, ah, aw;
  double *p, *ap;

  h = nRows;
  ah = aMatrix.giveNumberOfRows();
  aw = aMatrix.giveNumberOfColumns();

#ifdef MFEM_DEBUG
  if (co + aw > nColumns || ro + ah > h)
    pmtop_error("DenseMatrix::AddMatrix (...) 1");
#endif

  p = values.data() + ro + co * h;
  ap = aMatrix.givePointer();

  for (int c = 0; c < aw; c++) {
    for (int r = 0; r < ah; r++)
      p[r] += ap[r];
    p += h;
    ap += ah;
  }
}

void FloatMatrix ::beSubMatrixOf(const FloatMatrix &src,
                                 const IntArray &indxRow,
                                 const IntArray &indxCol)
/*
 * Modifies receiver to be a sub-matrix of the src matrix.
 * sub-matrix has size(indxRow) x size(indxCol) with values given as
 * this(i,j) = src( indxRow(i), indxCol(j) )
 */
{
#ifndef NDEBUG
  if ((indxRow.maximum()+1) > src.giveNumberOfRows() ||
      (indxCol.maximum()+1) > src.giveNumberOfColumns() || (indxRow.minimum()+1) < 1 ||
      (indxCol.minimum()+1) < 1) {
    std::cerr << "index exceeds source dimensions" << std::endl;
  }
#endif

  int szRow = indxRow.giveSize();
  int szCol = indxCol.giveSize();
  this->resize(szRow, szCol);

  for (int i = 0; i < szRow; i++) {
    for (int j = 0; j < szCol; j++) {
      (*this)(i, j) = src(indxRow(i), indxCol(j));
    }
  }
}
void FloatMatrix ::add(const FloatMatrix &aMatrix)
// Adds aMatrix to the receiver. If the receiver has a null size,
// adjusts its size to that of aMatrix. Returns the modified receiver.
{
  if (aMatrix.nRows == 0 || aMatrix.nColumns == 0) {
    return;
  }

  if (!this->isNotEmpty()) {
    this->operator=(aMatrix);
    return;
  }
#ifndef NDEBUG
  if ((aMatrix.nRows != nRows || aMatrix.nColumns != nColumns) &&
      aMatrix.isNotEmpty()) {
    std::cerr << "dimensions mismatch : (r1,c1)+(r2,c2) : (" << nRows << ","
              << nColumns << ")+(" << aMatrix.nRows << "," << aMatrix.nColumns
              << ")" << std::endl;
  }
#endif

#ifdef __MKL_MODULE
  MKL_INT aSize = aMatrix.nRows * aMatrix.nColumns;
  MKL_INT inc = 1;
  double s = 1.;
  cblas_daxpy(CblasColMajor, s, aMatrix.givePointer(), inc, this->givePointer(),
              inc);
#else
  for (size_t i = 0; i < this->values.size(); i++) {
    this->values[i] += aMatrix.values[i];
  }
#endif
}

void FloatMatrix ::add(double s, const FloatMatrix &aMatrix)
// Adds aMatrix to the receiver. If the receiver has a null size,
// adjusts its size to that of aMatrix. Returns the modified receiver.
{
  if (aMatrix.nRows == 0 || aMatrix.nColumns == 0) {
    return;
  }

  if (!this->isNotEmpty()) {
    this->operator=(aMatrix);
    this->times(s);
    return;
  }
#ifndef NDEBUG
  if ((aMatrix.nRows != nRows || aMatrix.nColumns != nColumns) &&
      aMatrix.isNotEmpty()) {
    std::cerr << "dimensions mismatch : (r1,c1)+(r2,c2) : (" << nRows << ","
              << nColumns << ")+(" << aMatrix.nRows << "," << aMatrix.nColumns
              << ")" << std::endl;
  }
#endif

#ifdef __MKL_MODULE
  MKL_INT aSize = aMatrix.nRows * aMatrix.nColumns;
  MKL_INT inc = 1;
  cblas_daxpy(CblasColMajor, s, aMatrix.givePointer(), inc, this->givePointer(),
              inc);
#else
  for (size_t i = 0; i < this->values.size(); i++) {
    this->values[i] += s * aMatrix.values[i];
  }
#endif
}
void FloatMatrix ::subtract(const FloatMatrix &aMatrix)
// Adds aMatrix to the receiver. If the receiver has a null size,
// adjusts its size to that of aMatrix. Returns the modified receiver.
{
  if (!this->isNotEmpty()) {
    this->operator=(aMatrix);
    this->negated();
    return;
  }
#ifndef NDEBUG
  if ((aMatrix.nRows != nRows || aMatrix.nColumns != nColumns) &&
      aMatrix.isNotEmpty()) {
    std::cerr << "dimensions mismatch : (r1,c1)-(r2,c2) : (" << nRows << ","
              << nColumns << ")-(" << aMatrix.nRows << "," << aMatrix.nColumns
              << ")" << std::endl;
  }
#endif

#ifdef __MKL_MODULE
  MKL_INT aSize = aMatrix.nRows * aMatrix.nColumns;
  MKL_INT inc = 1;
  double s = -1.;
  cblas_daxpy(this->values.size(), s, aMatrix.givePointer(), inc, this->givePointer(),
              inc);
#else
  for (size_t i = 0; i < this->values.size(); i++) {
    this->values[i] -= aMatrix.values[i];
  }
#endif
}
bool FloatMatrix ::solveForRhs(const FloatArray &b, FloatArray &answer,
                               bool transpose)
// solves equation b = this * x
{
#ifndef NDEBUG
  if (!this->isSquare()) {
    std::cerr << "cannot solve a " << nRows << " by " << nColumns << " matrix"
              << std::endl;
  }

  if (nRows != b.giveSize()) {
    std::cerr << "dimension mismatch" << std::endl;
  }
#endif

#ifdef __MKL_MODULE
  MKL_INT info, nrhs = 1;

  MKL_INT *ipiv = (MKL_INT *)mkl_malloc(this->nRows * sizeof(MKL_INT), 64);
  answer = b;
  info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, this->nRows, this->nColumns,
                        this->givePointer(), this->nRows, ipiv);
  if (info != 0) {
    throw std::runtime_error("dgetrf error " + std::to_string(info));
  }
  if (info == 0) {
    if (transpose) {
      info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'T', this->nRows, nrhs,
                            this->givePointer(), this->nRows, ipiv,
                            answer.givePointer(), nrhs);
    } else {
      info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', this->nRows, nrhs,
                            this->givePointer(), this->nRows, ipiv,
                            answer.givePointer(), nrhs);
    }
  }
  mkl_free(ipiv);
#else
  int pivRow;
  double piv, linkomb, help;
  FloatMatrix *mtrx, trans;
  if (transpose) {
    trans.beTranspositionOf(*this);
    mtrx = &trans;
  } else {
    mtrx = this;
  }

  answer = b;

  // initialize answer to be unity matrix;
  // lower triangle elimination by columns
  for (int i = 1; i < nRows; i++) {
    // find the suitable row and pivot
    piv = fabs(mtrx->at(i, i));
    pivRow = i;
    for (int j = i + 1; j <= nRows; j++) {
      if (fabs(mtrx->at(j, i)) > piv) {
        pivRow = j;
        piv = fabs(mtrx->at(j, i));
      }
    }

    if (piv < 1.e-20) {
      return false;
    }

    // exchange rows
    if (pivRow != i) {
      for (int j = i; j <= nRows; j++) {
        help = mtrx->at(i, j);
        mtrx->at(i, j) = mtrx->at(pivRow, j);
        mtrx->at(pivRow, j) = help;
      }
      help = answer.at(i);
      answer.at(i) = answer.at(pivRow);
      answer.at(pivRow) = help;
    }

    for (int j = i + 1; j <= nRows; j++) {
      linkomb = mtrx->at(j, i) / mtrx->at(i, i);
      for (int k = i; k <= nRows; k++) {
        mtrx->at(j, k) -= mtrx->at(i, k) * linkomb;
      }

      answer.at(j) -= answer.at(i) * linkomb;
    }
  }

  // back substitution
  for (int i = nRows; i >= 1; i--) {
    help = 0.;
    for (int j = i + 1; j <= nRows; j++) {
      help += mtrx->at(i, j) * answer.at(j);
    }

    answer.at(i) = (answer.at(i) - help) / mtrx->at(i, i);
  }
#endif

  return true;
}

bool FloatMatrix ::solveForRhs(const FloatMatrix &b, FloatMatrix &answer,
                               bool transpose)
// solves equation b = this * x
// returns x. this and b are kept untouched
//
// gaussian elimination - slow but safe
//
{
#ifndef NDEBUG
  if (!this->isSquare()) {
    std::cerr << "cannot solve a " << nRows << " by " << nColumns << " matrix"
              << std::endl;
  }

  if (nRows != b.giveNumberOfRows()) {
    std::cerr << "dimension mismatch" << std::endl;
  }
#endif

#ifdef __MKL_MODULE
  MKL_INT info;
  MKL_INT *ipiv = (MKL_INT *)mkl_malloc(nRows * sizeof(MKL_INT), 64);
  answer = b;
  info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, this->nRows, this->nColumns,
                        this->givePointer(), this->nRows, ipiv);
  if (info != 0) {
    throw std::runtime_error("dgetrf error " + std::to_string(info));
  }
  if (info == 0) {
    if (transpose) {
      info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'T', this->nRows, answer.nColumns,
                            this->givePointer(), this->nRows, ipiv,
                            answer.givePointer(), this->nRows);
    } else {
      info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', this->nRows, answer.nColumns,
                            this->givePointer(), this->nRows, ipiv,
                            answer.givePointer(), this->nRows);
    }
  }
  mkl_free(ipiv);
  return true;
#else
  int pivRow, nPs;
  double piv, linkomb, help;
  FloatMatrix *mtrx, trans;
  if (transpose) {
    trans.beTranspositionOf(*this);
    mtrx = &trans;
  } else {
    mtrx = this;
  }

  nPs = b.giveNumberOfColumns();
  answer = b;
  // initialize answer to be unity matrix;
  // lower triangle elimination by columns
  for (int i = 1; i < nRows; i++) {
    // find the suitable row and pivot
    piv = fabs(mtrx->at(i, i));
    pivRow = i;
    for (int j = i + 1; j <= nRows; j++) {
      if (fabs(mtrx->at(j, i)) > piv) {
        pivRow = j;
        piv = fabs(mtrx->at(j, i));
      }
    }

    if (fabs(piv) < 1.e-20) {
      return false; // PMTOP_ERROR("pivot too small, cannot solve %d by %d
                    // matrix", nRows, nColumns);
    }

    // exchange rows
    if (pivRow != i) {
      for (int j = i; j <= nRows; j++) {
        help = mtrx->at(i, j);
        mtrx->at(i, j) = mtrx->at(pivRow, j);
        mtrx->at(pivRow, j) = help;
      }

      for (int j = 1; j <= nPs; j++) {
        help = answer.at(i, j);
        answer.at(i, j) = answer.at(pivRow, j);
        answer.at(pivRow, j) = help;
      }
    }

    for (int j = i + 1; j <= nRows; j++) {
      linkomb = mtrx->at(j, i) / mtrx->at(i, i);
      for (int k = i; k <= nRows; k++) {
        mtrx->at(j, k) -= mtrx->at(i, k) * linkomb;
      }

      for (int k = 1; k <= nPs; k++) {
        answer.at(j, k) -= answer.at(i, k) * linkomb;
      }
    }
  }

  // back substitution
  for (int i = nRows; i >= 1; i--) {
    for (int k = 1; k <= nPs; k++) {
      help = 0.;
      for (int j = i + 1; j <= nRows; j++) {
        help += mtrx->at(i, j) * answer.at(j, k);
      }

      answer.at(i, k) = (answer.at(i, k) - help) / mtrx->at(i, i);
    }
  }
  return true;
#endif
}
void FloatMatrix ::initFromVector(const FloatArray &vector, bool transposed)
//
// constructor : creates (vector->giveSize(),1) FloatMatrix
// if transpose = 1 creates (1,vector->giveSize()) FloatMatrix
//
{
  if (transposed) {
    this->nRows = 1;
    this->nColumns = vector.giveSize();
  } else {
    this->nRows = vector.giveSize();
    this->nColumns = 1;
  }

  this->values = vector.values;
}

void FloatMatrix ::zero() {
  std ::fill(this->values.begin(), this->values.end(), 0.);
}

void FloatMatrix ::beUnitMatrix() {
#ifndef NDEBUG
  if (!this->isSquare()) {
    std::cerr << "cannot make unit matrix of " << nRows << " by " << nColumns
              << " matrix" << std::endl;
  }
#endif

  this->zero();
  for (int i = 1; i <= nRows; i++) {
    this->at(i, i) = 1.0;
  }
}
void FloatMatrix ::bePinvID()
// this matrix is the product of the 6x6 deviatoric projection matrix ID
// and the inverse scaling matrix Pinv
{
  this->resize(6, 6);
  values[0] = values[7] = values[14] = 2. / 3.;
  values[1] = values[2] = values[6] = values[8] = values[12] = values[13] =
      -1. / 3.;
  values[21] = values[28] = values[35] = 0.5;
}
void FloatMatrix ::resize(int rows, int columns)
//
// resizes receiver, all data will be lost
//
{
  this->nRows = rows;
  this->nColumns = columns;
  this->values.assign(rows * columns, 0.);
}
void FloatMatrix ::resizeWithData(int rows, int columns)
//
// resizes receiver, all data kept
//
{
  // Check of resize if necessary at all.
  if (rows == this->nRows && columns == this->nColumns) {
    return;
  }

  FloatMatrix old(std ::move(*this));

  this->nRows = rows;
  this->nColumns = columns;
  this->values.resize(rows * columns);

  int ii = PMTop::min(rows, old.giveNumberOfRows());
  int jj = PMTop::min(columns, old.giveNumberOfColumns());
  // copy old values if possible
  for (int i = 1; i <= ii; i++) {
    for (int j = 1; j <= jj; j++) {
      this->at(i, j) = old.at(i, j);
    }
  }
}
void FloatMatrix ::hardResize(int rows, int columns)
//
// resizes receiver, all data will be lost
//
{
  this->nRows = rows;
  this->nColumns = columns;
  values.assign(rows * columns, 0.);
  this->values.shrink_to_fit();
}
double FloatMatrix ::giveDeterminant() const
// Returns the determinant of the receiver.
{
#ifndef NDEBUG
  if (!this->isSquare()) {
    std::cout << "cannot compute the determinant of a non-square " << nRows
              << " by " << nColumns << " matrix" << std::endl;
  }
#endif

  if (nRows == 1) {
    return values[0];
  } else if (nRows == 2) {
    return (values[0] * values[3] - values[1] * values[2]);
  } else if (nRows == 3) {
    return (
        values[0] * values[4] * values[8] + values[3] * values[7] * values[2] +
        values[6] * values[1] * values[5] - values[6] * values[4] * values[2] -
        values[7] * values[5] * values[0] - values[8] * values[3] * values[1]);
  } else {
    std::cout << "cannot compute the determinant of a matrix larger than 3x3"
              << std::endl;
  }

  return 0.;
}
void FloatMatrix ::beDiagonal(const FloatArray &diag) {
  int n = diag.giveSize();
  this->resize(n, n);
  for (int i = 0; i < n; ++i) {
    (*this)(i, i) = diag[i];
  }
}
double FloatMatrix ::giveTrace() const {
#ifndef NDEBUG
  if (!this->isSquare()) {
    std::cerr << "cannot compute the trace of a non-square " << nRows << " by "
              << nColumns << " matrix" << std::endl;
  }
#endif
  double answer = 0.;
  for (int k = 0; k < nRows; k++) {
    answer += values[k * (nRows + 1)];
  }
  return answer;
}
void FloatMatrix ::printYourself() const
// Prints the receiver on screen.
{
  printf("FloatMatrix with dimensions : %d %d\n", nRows, nColumns);
  if (nRows <= 250 && nColumns <= 250) {
    for (int i = 1; i <= nRows; ++i) {
      for (int j = 1; j <= nColumns && j <= 100; ++j) {
        printf("%10.3e  ", this->at(i, j));
      }

      printf("\n");
    }
  } else {
    printf("   large matrix : coefficients not printed \n");
  }
}
void FloatMatrix ::printYourselfToFile(const std::string filename,
                                       const bool showDimensions) const
// Prints the receiver to file.
{
  std ::ofstream matrixfile(filename);
  if (matrixfile.is_open()) {
    if (showDimensions)
      matrixfile << "FloatMatrix with dimensions : " << nRows << ", "
                 << nColumns << "\n";
    matrixfile << std::scientific << std::right << std::setprecision(3);
    for (int i = 1; i <= nRows; ++i) {
      for (int j = 1; j <= nColumns; ++j) {
        matrixfile << std::setw(10) << this->at(i, j) << "\t";
      }

      matrixfile << "\n";
    }
    matrixfile.close();
  } else {
    std::cerr << "Failed to write to file" << std::endl;
  }
}
void FloatMatrix ::printYourself(const std::string &name) const
// Prints the receiver on screen.
{
  printf("%s (%d x %d): \n", name.c_str(), nRows, nColumns);
  if (nRows <= 250 && nColumns <= 250) {
    for (int i = 1; i <= nRows; ++i) {
      for (int j = 1; j <= nColumns && j <= 100; ++j) {
        printf("%10.3e  ", this->at(i, j));
      }

      printf("\n");
    }
  } else {
    for (int i = 1; i <= nRows && i <= 20; ++i) {
      for (int j = 1; j <= nColumns && j <= 10; ++j) {
        printf("%10.3e  ", this->at(i, j));
      }
      if (nColumns > 10)
        printf(" ...");
      printf("\n");
    }
    if (nRows > 20)
      printf(" ...\n");
  }
}
void FloatMatrix ::pY() const
// Prints the receiver on screen with higher accuracy than printYourself.
{
  printf("[");
  for (int i = 1; i <= nRows; ++i) {
    for (int j = 1; j <= nColumns; ++j) {
      printf("%20.15e", this->at(i, j));
      if (j < nColumns) {
        printf(",");
      } else {
        printf(";");
      }
    }
  }

  printf("];\n");
}
void FloatMatrix ::writeCSV(const std ::string &name) const {
  FILE *file = fopen(name.c_str(), "w");
  for (int i = 1; i <= nRows; ++i) {
    for (int j = 1; j <= nColumns; ++j) {
      fprintf(file, "%10.3e, ", this->at(i, j));
    }

    fprintf(file, "\n");
  }
  fclose(file);
}
void FloatMatrix ::rotatedWith(const FloatMatrix &r, char mode)
// Returns the receiver 'a' rotated according the change-of-base matrix r.
// The method performs the operation  a = r^T . a . r . or the inverse
{
  FloatMatrix rta;

  if (mode == 'n') {
    rta.beTProductOf(r, *this); //  r^T . a
    this->beProductOf(rta, r);  //  r^T . a . r
  } else if (mode == 't') {
    rta.beProductOf(r, *this);  //  r . a
    this->beProductTOf(rta, r); //  r . a . r^T
  } else {
    std::cerr << "unsupported mode" << std::endl;
  }
}
void FloatMatrix ::symmetrized()
// Initializes the lower half of the receiver to the upper half.
{
#ifndef NDEBUG
  if (nRows != nColumns) {
    std::cerr << "cannot symmetrize a non-square matrix" << std::endl;
  }

#endif

  for (int i = 2; i <= nRows; i++) {
    for (int j = 1; j < i; j++) {
      this->at(i, j) = this->at(j, i);
    }
  }
}
void FloatMatrix ::times(double factor)
// Multiplies every coefficient of the receiver by factor. Answers the
// modified receiver.
{
#ifdef __MKL_MODULE
  cblas_dscal(this->values.size(),factor,this->givePointer(),1);
#else
  for (double &x : this->values) {
    x *= factor;
  }
#endif
  // dscal_ seemed to be slower for typical usage of this function.
}
void FloatMatrix ::negated() {
#ifdef __MKL_MODULE
  cblas_dscal(this->values.size(),-1.0,this->givePointer(),1);
#else
  for (double &x : this->values) {
    x = -x;
  }
#endif
}
double FloatMatrix ::computeFrobeniusNorm() const {
  return sqrt(std ::inner_product(this->values.begin(), this->values.end(),
                                  this->values.begin(), 0.));
}
double FloatMatrix ::computeNorm(char p) const {
#ifdef __MKL_MODULE
  FloatArray work(this->giveNumberOfRows());
  int lda = PMTop::max(this->nRows, 1);
  double norm = LAPACKE_dlange(LAPACK_COL_MAJOR, p, this->nRows, this->nColumns,
                               this->givePointer(), lda);
  return norm;

#else
  if (p == '1') { // Maximum absolute column sum.
    double col_sum, max_col = 0.0;
    for (int j = 1; j <= this->nColumns; j++) {
      col_sum = 0.0;
      for (int i = 1; i <= this->nRows; i++) {
        col_sum += fabs(this->at(i, j));
      }
      if (col_sum > max_col) {
        max_col = col_sum;
      }
    }
    return max_col;
  }
  ///@todo Use this when obtaining eigen values is implemented.
  /*else if (p == '2') {
   *  double lambda_max;
   *  FloatMatrix AtA;
   *  FloatArray eigs;
   *  AtA.beTProductOf(*this,this);
   *  Ata.eigenValues(eigs, 1);
   *  return sqrt(eigs(0));
   * } */
  else {
    std::cerr << "p == " << p << " not implemented." << std::endl;
    return 0.0;
  }
#endif
}

void FloatMatrix ::beMatrixForm(const FloatArray &aArray) {
  // Revrites the vector on matrix form (symmetrized matrix used if size is 6),
  // order: 11, 22, 33, 23, 13, 12
  // order: 11, 22, 33, 23, 13, 12, 32, 31, 21
#ifndef NDEBUG
  if (aArray.giveSize() != 6 && aArray.giveSize() != 9) {
    std::cerr << "array size mismatch" << std::endl;
  }
#endif
  this->resize(3, 3);
  if (aArray.giveSize() == 9) {
    this->at(1, 1) = aArray.at(1);
    this->at(2, 2) = aArray.at(2);
    this->at(3, 3) = aArray.at(3);
    this->at(2, 3) = aArray.at(4);
    this->at(1, 3) = aArray.at(5);
    this->at(1, 2) = aArray.at(6);
    this->at(3, 2) = aArray.at(7);
    this->at(3, 1) = aArray.at(8);
    this->at(2, 1) = aArray.at(9);
  } else if (aArray.giveSize() == 6) {
    this->at(1, 1) = aArray.at(1);
    this->at(2, 2) = aArray.at(2);
    this->at(3, 3) = aArray.at(3);
    this->at(2, 3) = aArray.at(4);
    this->at(1, 3) = aArray.at(5);
    this->at(1, 2) = aArray.at(6);
    this->at(3, 2) = aArray.at(4);
    this->at(3, 1) = aArray.at(5);
    this->at(2, 1) = aArray.at(6);
  }
}

void FloatMatrix::GetSubVector(const IntArray &row, int col, FloatArray &Result) const {
  Result.resize(row.giveSize());
  for (int i=0;i<row.giveSize();i++) {
    Result[i]=(*this)(row[i],col);
  }
}

void FloatMatrix::GetSubVector(const IntArray &row, IntArray &col, FloatMatrix &Result) const {
  Result.resize(row.giveSize(),col.giveSize());
  for (int i=0;i<row.giveSize();i++) {
    for (int j=0;j<col.giveSize();j++) {
      Result(i,j)=(*this)(row[i],col[j]);
    }
  }
}



void FloatMatrix ::changeComponentOrder() {
  // Changes index order between abaqus <-> PMTOP
  // #  ifndef NDEBUG
  //     if ( nRows != 6 || nColumns != 6 ) {
  //         PMTOP_ERROR("matrix dimension is not 6x6");
  //     }
  // #  endif

  if (nRows == 6 && nColumns == 6) {
    // This could probably be done more beautifully + efficiently.

    std ::swap(this->at(4, 1), this->at(6, 1));

    std ::swap(this->at(4, 2), this->at(6, 2));

    std ::swap(this->at(4, 3), this->at(6, 3));

    std ::swap(this->at(1, 4), this->at(1, 6));
    std ::swap(this->at(2, 4), this->at(2, 6));
    std ::swap(this->at(3, 4), this->at(3, 6));
    std ::swap(this->at(4, 4), this->at(6, 6));
    std ::swap(this->at(5, 4), this->at(5, 6));
    std ::swap(this->at(6, 4), this->at(4, 6));

    std ::swap(this->at(4, 5), this->at(6, 5));
  } else if (nRows == 9 && nColumns == 9) {
    // PMTOP:           11, 22, 33, 23, 13, 12, 32, 31, 21
    // UMAT:            11, 22, 33, 12, 13, 23, 32, 21, 31
    const int abq2oo[9] = {1, 2, 3, 6, 5, 4, 7, 9, 8};

    FloatMatrix tmp(9, 9);
    for (int i = 1; i <= 9; i++) {
      for (int j = 1; j <= 9; j++) {
        tmp.at(i, j) = this->at(abq2oo[i - 1], abq2oo[j - 1]);
      }
    }

    *this = tmp;
  }
}
double FloatMatrix ::computeReciprocalCondition(char p) const {
#ifndef NDEBUG
  if (!this->isSquare()) {
    std::cerr << "receiver must be square (is " << this->nRows << " by "
              << this->nColumns << ")" << std::endl;
  }
#endif
  double anorm = this->computeNorm(p);

#ifdef __MKL_MODULE
  MKL_INT n = this->nRows;

  MKL_INT info;
  double rcond;
  if (n > 3) {
    MKL_INT *iwork = (MKL_INT *)mkl_malloc(n * sizeof(MKL_INT), 64);
    FloatMatrix a_cpy = *this;
    info =
        LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, a_cpy.givePointer(), n, iwork);
    if (info < 0) {
      throw std::runtime_error("dgetrf error " + std::to_string(info));
    }
    info = LAPACKE_dgecon(LAPACK_COL_MAJOR, p, this->giveNumberOfRows(),
                          a_cpy.givePointer(), this->giveNumberOfRows(), anorm,
                          &rcond);
    if (info < 0) {
      throw std::runtime_error("dgecon error " + std::to_string(info));
    }
    mkl_free(iwork);
    return rcond;
  }

#endif
  if (this->giveDeterminant() <= 1e-6 * anorm) {
    return 0.0;
  }
  FloatMatrix inv;
  inv.beInverseOf(*this);
  return 1.0 / (inv.computeNorm(p) * anorm);
}
void FloatMatrix ::beMatrixFormOfStress(const FloatArray &aArray) {
  // Revrites the  matrix on vector form (symmetrized matrix used), order: 11,
  // 22, 33, 23, 13, 12
#ifndef NDEBUG
  if (aArray.giveSize() != 6 && aArray.giveSize() != 9) {
    std::cerr << "matrix dimension is not 3x3" << std::endl;
  }
#endif
  this->resize(3, 3);
  if (aArray.giveSize() == 9) {
    this->at(1, 1) = aArray.at(1);
    this->at(2, 2) = aArray.at(2);
    this->at(3, 3) = aArray.at(3);
    this->at(2, 3) = aArray.at(4);
    this->at(1, 3) = aArray.at(5);
    this->at(1, 2) = aArray.at(6);
    this->at(3, 2) = aArray.at(7);
    this->at(3, 1) = aArray.at(8);
    this->at(2, 1) = aArray.at(9);
  } else if (aArray.giveSize() == 6) {
    this->at(1, 1) = aArray.at(1);
    this->at(2, 2) = aArray.at(2);
    this->at(3, 3) = aArray.at(3);
    this->at(2, 3) = aArray.at(4);
    this->at(1, 3) = aArray.at(5);
    this->at(1, 2) = aArray.at(6);
    this->at(3, 2) = aArray.at(4);
    this->at(3, 1) = aArray.at(5);
    this->at(2, 1) = aArray.at(6);
  }
}

#if 0
bool FloatMatrix :: computeEigenValuesSymmetric(FloatArray &lambda, FloatMatrix &v, int neigs) const
{
#ifdef __MKL_MODULE
    double abstol = 1.0;
    MKL_INT lda, n, ldz, info, found;
    n = this->nRows;
    lda = n;
    ldz = n;
    FloatMatrix a;
    a = * this;
    if ( neigs == 0 )
    {
        neigs = n;
    }
    lambda.resize(neigs);
    v.resize(n, neigs);
    IntArray ifail(n);
    if ( neigs > 0 )
    {
        int one = 1;
        info=LAPACKE_dsyevx(LAPACK_COL_MAJOR,'N','I','U',n,a.givePointer(),
        n,0.,0.,one,neigs,abstol,&found,lambda.givePointer(),v.givePointer(),ldz,ifail.givePointer());
    }
    else
    {
        info=LAPACKE_dsyevx(LAPACK_COL_MAJOR,'A','A','U',n,a.givePointer(),n,0.,0.,0.,0.,abstol,&found,lambda.givePointer(),v.givePointer(),ldz,ifail.givePointer());
    }
    return info == 0;

#else
    PMTOP_ERROR("Requires MKL");
    return false;

#endif
}
#endif
/*
contextIOResultType FloatMatrix ::storeYourself(DataStream &stream) const
// writes receiver's binary image into stream
// use id to distinguish some instances
// return value >0 success
//              =0 file i/o error
{
  // write size
  if (!stream.write(nRows)) {
    return (CIO_IOERR);
  }

  if (!stream.write(nColumns)) {
    return (CIO_IOERR);
  }

  // write raw data
  if (!stream.write(this->givePointer(), nRows * nColumns)) {
    return (CIO_IOERR);
  }

  // return result back
  return CIO_OK;
}
contextIOResultType FloatMatrix ::restoreYourself(DataStream &stream)
// reads receiver from stream
// warning - overwrites existing data!
// returns 0 if file i/o error
//        -1 if id of class id is not correct
{
  // read size
  if (!stream.read(nRows)) {
    return (CIO_IOERR);
  }

  if (!stream.read(nColumns)) {
    return (CIO_IOERR);
  }

  this->values.resize(nRows * nColumns);

  // read raw data
  if (!stream.read(this->givePointer(), nRows * nColumns)) {
    return (CIO_IOERR);
  }

  // return result back
  return CIO_OK;
}
int FloatMatrix ::givePackSize(DataStream &buff) const {
  return buff.givePackSizeOfInt(1) + buff.givePackSizeOfInt(1) +
         buff.givePackSizeOfDouble(nRows * nColumns);
}
*/
bool FloatMatrix ::jaco_(FloatArray &eval, FloatMatrix &v, int nf) {
  /*
   * Solves the eigenvalues and eigenvectors of real
   * symmetric matrix by jacobi method.
   *  Written by bp. Inspired by ED WILSON jaco_ procedure.
   *
   * Parameters (input):
   * nf - number of significant figures
   *
   * Output params:
   * eval - eigen values (not sorted)
   * v    - eigenvectors (stored columvise)
   */

  /* Local variables */
  int neq = this->giveNumberOfRows();

  double c_b2 = .10;
  // double c_b27 = .01;

  /* Function Body */
#ifndef NDEBUG
  if (!isSquare()) {
    std::cerr << "Not square matrix" << std::endl;
  }
  // check for symmetry
  for (int i = 1; i <= neq; i++) {
    for (int j = i + 1; j <= neq; j++) {
      // if ( this->at(i, j) != this->at(j, i) ) {
      if (fabs(this->at(i, j) - this->at(j, i)) > 1.0e-6) {
        std::cerr << "Not Symmetric matrix" << std::endl;
      }
    }
  }

#endif

  eval.resize(neq);
  v.resize(neq, neq);

  for (int i = 1; i <= neq; i++) {
    eval.at(i) = this->at(i, i);
  }

  double tol = pow(c_b2, nf);
  double sum = 0.0;
  for (int i = 1; i <= neq; ++i) {
    for (int j = 1; j <= neq; ++j) {
      sum += fabs(this->at(i, j));
      v.at(i, j) = 0.0;
    }

    v.at(i, i) = 1.0;
  }

  if (sum <= 0.0) {
    return 0;
  }

  /* ---- REDUCE MATRIX TO DIAGONAL ---------------- */
  int ite = 0;
  double ssum;
  do {
    ssum = 0.0;
    for (int j = 2; j <= neq; ++j) {
      int ih = j - 1;
      for (int i = 1; i <= ih; ++i) {
        if ((fabs(this->at(i, j)) / sum) > tol) {
          ssum += fabs(this->at(i, j));
          /* ---- CALCULATE ROTATION ANGLE ----------------- */
          double aa =
              atan2(this->at(i, j) * 2.0, eval.at(i) - eval.at(j)) / 2.0;
          double si = sin(aa);
          double co = cos(aa);
          /*
           *   // ---- MODIFY "I" AND "J" COLUMNS OF "A" AND "V"
           *   for (k = 1; k <= neq; ++k) {
           *    tt = this->at(k, i);
           *    this->at(k, i) = co * tt + si * this->at(k, j);
           *    this->at(k, j) = -si * tt + co * this->at(k, j);
           *    tt = v.at(k, i);
           *    v.at(k, i) = co * tt + si * v.at(k, j);
           *    // L500:
           *    v.at(k, j) = -si * tt + co * v.at(k, j);
           *   }
           *   // ---- MODIFY DIAGONAL TERMS --------------------
           *   this->at(i, i) = co * this->at(i, i) + si * this->at(j, i);
           *   this->at(j, j) = -si * this->at(i, j) + co * this->at(j, j);
           *   this->at(i, j) = 0.0;
           *   // ---- MAKE "A" MATRIX SYMMETRICAL --------------
           *   for (k = 1; k <= neq; ++k) {
           *    this->at(i, k) = this->at(k, i);
           *    this->at(j, k) = this->at(k, j);
           *    // L600:
           *   }
           */
          // ---- MODIFY "I" AND "J" COLUMNS OF "A" AND "V"
          for (int k = 1; k < i; ++k) {
            double tt = this->at(k, i);
            this->at(k, i) = co * tt + si * this->at(k, j);
            this->at(k, j) = -si * tt + co * this->at(k, j);
            tt = v.at(k, i);
            v.at(k, i) = co * tt + si * v.at(k, j);
            v.at(k, j) = -si * tt + co * v.at(k, j);
          }

          // diagonal term (i,i)
          double tt = eval.at(i);
          eval.at(i) = co * tt + si * this->at(i, j);
          double aij = -si * tt + co * this->at(i, j);
          tt = v.at(i, i);
          v.at(i, i) = co * tt + si * v.at(i, j);
          v.at(i, j) = -si * tt + co * v.at(i, j);

          for (int k = i + 1; k < j; ++k) {
            double tt = this->at(i, k);
            this->at(i, k) = co * tt + si * this->at(k, j);
            this->at(k, j) = -si * tt + co * this->at(k, j);
            tt = v.at(k, i);
            v.at(k, i) = co * tt + si * v.at(k, j);
            v.at(k, j) = -si * tt + co * v.at(k, j);
          }

          // diagonal term (j,j)
          tt = this->at(i, j);
          double aji = co * tt + si * eval.at(j);
          eval.at(j) = -si * tt + co * eval.at(j);

          tt = v.at(j, i);
          v.at(j, i) = co * tt + si * v.at(j, j);
          v.at(j, j) = -si * tt + co * v.at(j, j);
          //
          for (int k = j + 1; k <= neq; ++k) {
            double tt = this->at(i, k);
            this->at(i, k) = co * tt + si * this->at(j, k);
            this->at(j, k) = -si * tt + co * this->at(j, k);
            tt = v.at(k, i);
            v.at(k, i) = co * tt + si * v.at(k, j);
            v.at(k, j) = -si * tt + co * v.at(k, j);
          }

          // ---- MODIFY DIAGONAL TERMS --------------------
          eval.at(i) = co * eval.at(i) + si * aji;
          eval.at(j) = -si * aij + co * eval.at(j);
          this->at(i, j) = 0.0;
        } else {
          /* ---- A(I,J) MADE ZERO BY ROTATION ------------- */
          ;
        }
      }
    }

    /* ---- CHECK FOR CONVERGENCE -------------------- */
    if (++ite > 50) {
      std::cerr << "too many iterations" << std::endl;
    }
  } while (fabs(ssum) / sum > tol);

  // restore original matrix
  for (int i = 1; i <= neq; i++) {
    for (int j = i; j <= neq; j++) {
      this->at(i, j) = this->at(j, i);
    }
  }

  return 0;
} /* jaco_ */

std ::ostream &operator<<(std ::ostream &out, const FloatMatrix &x) {
  out << x.nRows << " " << x.nColumns << " {";
  for (int i = 0; i < x.nRows; ++i) {
    for (int j = 0; j < x.nColumns; ++j) {
      out << " " << x(i, j);
    }
    out << ";";
  }
  out << "}";
  return out;
}

FloatMatrix &operator*=(FloatMatrix &x, const double &a) {
  x.times(a);
  return x;
}
FloatMatrix operator*(const FloatMatrix &a, const FloatMatrix &b) {
  FloatMatrix ans;
  ans.beProductOf(a, b);
  return ans;
}
FloatArray operator*(const FloatMatrix &a, const FloatArray &b) {
  FloatArray ans;
  ans.beProductOf(a, b);
  return ans;
}

FloatMatrix operator+(const FloatMatrix &a, const FloatMatrix &b) {
  FloatMatrix ans(a);
  ans.add(b);
  return ans;
}
FloatMatrix operator-(const FloatMatrix &a, const FloatMatrix &b) {
  FloatMatrix ans(a);
  ans.subtract(b);
  return ans;
}
FloatMatrix &operator+=(FloatMatrix &a, const FloatMatrix &b) {
  a.add(b);
  return a;
}
FloatMatrix &operator-=(FloatMatrix &a, const FloatMatrix &b) {
  a.subtract(b);
  return a;
}


bool FloatMatrix::readTXT(const std::string &filename) {
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    return false;
  }
  std::vector<std::vector<double>> rowData;
  std::string line;
  while (std::getline(infile, line)) {
    if (line.empty())
      continue;
    std::istringstream iss(line);
    std::vector<double> row;
    double value;
    while (iss >> value) {
      row.push_back(value);
    }
    rowData.push_back(row);
  }
  if (rowData.empty())
    return false;

  // 检查所有行的列数是否一致
  int numRows = static_cast<int>(rowData.size());
  int numCols = static_cast<int>(rowData[0].size());
  for (const auto &row : rowData) {
    if (static_cast<int>(row.size()) != numCols)
      return false;
  }
  // 设置矩阵的维度，并分配空间
  this->nRows = numRows;
  this->nColumns = numCols;
  values.resize(numRows * numCols);

  // 将数据转换为列优先顺序存储
  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      values[j * numRows + i] = rowData[i][j];
    }
  }
  return true;
}