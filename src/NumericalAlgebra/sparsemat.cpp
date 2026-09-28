// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#include "sparsemat.hpp"
#include "error.hpp"
#include "floatarray.h"
#include "floatmatrix.h"
#include "intarray.h"
#include <iomanip>
#include <mkl_spblas.h>
#include <array.hpp>
#include <assert.h>
#include <sort_pairs.hpp>
#include "array.cpp"

SparseMatrix::SparseMatrix(int nrows, int ncols) : Matrix(nrows) {
  RowIndex = NULL;
  J = NULL;
  A = NULL;

  Rows = new RowNode *[nrows];
  width = (ncols) ? (ncols) : (nrows);
  for (int i = 0; i < nrows; i++)
    Rows[i] = NULL;
  mkl_sparse_matrix = NULL;
}

int SparseMatrix::RowSize(int i) {
  if (RowIndex)
    return RowIndex[i + 1] - RowIndex[i];

  int s = 0;
  RowNode *row = Rows[i];
  for (; row != NULL; row = row->Prev)
    if (row->Value != 0.0)
      s++;
  return s;
}

double &SparseMatrix::Elem(int i, int j) { return operator()(i, j); }

const double &SparseMatrix::Elem(int i, int j) const {
  return operator()(i, j);
}

double &SparseMatrix::operator()(int i, int j) {
  int k, end;

#ifdef NDEBUG
  if (i >= size || i < 0 || j >= width || j < 0)
    pmtop_error("SparseMatrix::operator() #1");
#endif

  if (A == NULL)
    pmtop_error("SparseMatrix::operator() #2");

  end = RowIndex[i + 1];
  for (k = RowIndex[i]; k < end; k++)
    if (J[k] == j)
      return A[k];

  pmtop_error("SparseMatrix::operator() #3");
  return A[0];
}

const double &SparseMatrix::operator()(int i, int j) const {
  int k, end;
  static const double zero = 0.0;

#ifdef NDEBUG
  if (i >= size || i < 0 || j >= width || j < 0)
    pmtop_error("SparseMatrix::operator() const #1");
#endif

  if (A == NULL)
    pmtop_error("SparseMatrix::operator() const #2");
  end = RowIndex[i + 1];
  for (k = RowIndex[i]; k < end; k++)
    if (J[k] == j)
      return A[k];

  return zero;
}

void SparseMatrix::Mult(const FloatArray &x, FloatArray &y) const {
  y = 0.0;
  AddMult(x, y);
}

void SparseMatrix::Mult(const double * x, double * y) const {
#ifdef __MKL_MODULE
  if (Use_MKL == true)
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, mkl_sparse_matrix,
                    mkl_descr, x, 0.0, y);
#endif
  if (Use_MKL == false) {
    pmtop_error("SparseMatrix::Mult");
  }
}

void SparseMatrix::AddMult(const FloatArray &x, FloatArray &y,
                           const double a) const {
#ifdef NDEBUG
  if ((width != x.giveSize()) || (size != y.giveSize()))
    pmtop_error("SparseMatrix::AddMult() #1");
#endif

#ifdef __MKL_MODULE
  if (Use_MKL == true)
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, a, mkl_sparse_matrix,
                    mkl_descr, x.givePointer(), 0.0, y.givePointer());
#endif
  if (Use_MKL == false) {
    int i, j, end;
    double *Ap = A, *yp = y.givePointer();
    const double *xp = x.givePointer();

    if (Ap == NULL) {
      //  The matrix is not finalized, but multiplication is still possible
      for (i = 0; i < size; i++) {
        RowNode *row = Rows[i];
        double b = 0.0;
        for (; row != NULL; row = row->Prev)
          b += row->Value * xp[row->Column];
        *yp += a * b;
        yp++;
      }
      return;
    }

    int *Jp = J, *Ip = RowIndex;

    j = *Ip;
    if (a == 1.0)
      for (i = 0; i < size; i++) {
        double d;
        d = 0.0;
        Ip++;
        end = (*Ip);
        for (; j < end; j++) {
          d += (*Ap) * xp[*Jp];
          Ap++;
          Jp++;
        }
        *yp += d;
        yp++;
      }
    else
      for (i = 0; i < size; i++) {
        double d;
        d = 0.0;
        Ip++;
        end = (*Ip);
        for (; j < end; j++) {
          d += (*Ap) * xp[*Jp];
          Ap++;
          Jp++;
        }
        *yp += a * d;
        yp++;
      }
  }
}

void SparseMatrix::MultTranspose(const FloatArray &x, FloatArray &y) const {
  y = 0.0;
  AddMultTranspose(x, y);
}

void SparseMatrix::AddMultTranspose(const FloatArray &x, FloatArray &y,
                                    const double a) const {
#ifdef NDEBUG
  if ((size != x.giveSize()) || (width != y.giveSize()))
    pmtop_error("SparseMatrix::AddMultTranspose() #1");
#endif

#ifdef __MKL_MODULE
  if (Use_MKL == true)
    mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, a, mkl_sparse_matrix, mkl_descr,
                    x.givePointer(), 0.0, y.givePointer());
#endif

  if (Use_MKL == false) {
    int i, j, end;
    double *yp = y.givePointer();

    if (A == NULL) {
      // The matrix is not finalized, but multiplication is still possible
      for (i = 0; i < size; i++) {
        RowNode *row = Rows[i];
        double b = a * x(i);
        for (; row != NULL; row = row->Prev)
          yp[row->Column] += row->Value * b;
      }
      return;
    }

    for (i = 0; i < size; i++) {
      double xi = a * x(i);
      end = RowIndex[i + 1];
      for (j = RowIndex[i]; j < end; j++) {
        yp[J[j]] += A[j] * xi;
      }
    }
  }
}

void SparseMatrix::PartMult(const IntArray &rows, const FloatArray &x,
                            FloatArray &y) {
  if (A) {
    for (int i = 0; i < rows.giveSize(); i++) {
      int r = rows[i];
      int end = RowIndex[r + 1];
      double a = 0.0;
      for (int j = RowIndex[r]; j < end; j++)
        a += A[j] * x(J[j]);
      y(r) = a;
    }
  } else {
    pmtop_error("SparseMatrix::PartMult");
  }
}

double SparseMatrix::InnerProduct(const FloatArray &x,
                                  const FloatArray &y) const {

  double prod = 0.0;
#ifdef __MKL_MODULE
  if (Use_MKL == true) {
    FloatArray a(size);
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, mkl_sparse_matrix,
                    mkl_descr, x.givePointer(), 0.0, a.givePointer());
    prod = cblas_ddot(size, a.givePointer(), 1, y.givePointer(), 1);
  }
#endif
  if (Use_MKL == false) {
    for (int i = 0; i < size; i++) {
      double a = 0.0;
      if (A)
        for (int j = RowIndex[i], end = RowIndex[i + 1]; j < end; j++)
          a += A[j] * x(J[j]);
      else
        for (RowNode *node_p = Rows[i]; node_p != NULL; node_p = node_p->Prev)
          a += node_p->Value * x(node_p->Column);
      prod += a * y(i);
    }
  }
  return prod;
}

MatrixInverse *SparseMatrix::Inverse() const { return NULL; }

// 该函数用于对稀疏矩阵的某一行进行消元操作，同时更新右端项向量 rhs。
void SparseMatrix::EliminateRow(int row, const double sol, FloatArray &rhs) {
  RowNode *aux;

#ifdef NDEBUG
  if (row >= size || row < 0)
    pmtop_error("SparseMatrix::EliminateRow () #1");
#endif

  if (Rows == NULL)
    pmtop_error("SparseMatrix::EliminateRow () #2");

  for (aux = Rows[row]; aux != NULL; aux = aux->Prev) {
    rhs(aux->Column) -= sol * aux->Value;
    aux->Value = 0.0;
  }
}

void SparseMatrix::EliminateRow(int row) {
  RowNode *aux;

#ifdef NDEBUG
  if (row >= size || row < 0)
    pmtop_error("SparseMatrix::EliminateRow () #1");
#endif

  if (Rows == NULL)
    pmtop_error("SparseMatrix::EliminateRow () #2");

  for (aux = Rows[row]; aux != NULL; aux = aux->Prev)
    aux->Value = 0.0;
}

void SparseMatrix::EliminateCol(int col) {
  RowNode *aux;

  if (Rows == NULL)
    pmtop_error("SparseMatrix::EliminateCol () #1");

  for (int i = 0; i < size; i++)
    for (aux = Rows[i]; aux != NULL; aux = aux->Prev)
      if (aux->Column == col)
        aux->Value = 0.0;
}

void SparseMatrix::EliminateCols(IntArray &cols, FloatArray *x, FloatArray *b) {
  RowNode *aux;

  if (Rows == NULL)
    pmtop_error("SparseMatrix::EliminateCols () #1");

  for (int i = 0; i < size; i++)
    for (aux = Rows[i]; aux != NULL; aux = aux->Prev)
      if (cols[aux->Column]) {
        if (x && b)
          (*b)(i) -= aux->Value * (*x)(aux->Column);
        aux->Value = 0.0;
      }
}

// 该函数在稀疏矩阵中消除指定行列（rc）的非对角元素，同时更新右端项向量
// rhs。对于矩阵的两种存储格式（CSR 或链表），函数分别进行了处理。
// 消除稀疏矩阵的指定行 rc 和列 rc 的所有非对角元素。
// 根据参数 d，设置对角元素值： 如果 d ==1，对角元素更新为当前值乘以解 sol。
// 如果 d ==0，对角元素设置为 1.0，并更新 rhs 为解 sol。 更新
// rhs以反映消元操作对右端项向量的影响。
void SparseMatrix::EliminateRowCol(int rc, const double sol, FloatArray &rhs,
                                   int d) {
  int col;

#ifdef NDEBUG
  if (rc >= size || rc < 0)
    pmtop_error("SparseMatrix::EliminateRowCol () #1");
#endif

  if (Rows == NULL)
    for (int j = RowIndex[rc]; j < RowIndex[rc + 1]; j++)
      if ((col = J[j]) == rc)
        if (d) {
          rhs(rc) = A[j] * sol;
        } else {
          A[j] = 1.0;
          rhs(rc) = sol;
        }
      else {
        A[j] = 0.0;
        for (int k = RowIndex[col]; 1; k++)
          if (k == RowIndex[col + 1]) {
            pmtop_error("SparseMatrix::EliminateRowCol () #2");
          } else if (J[k] == rc) {
            rhs(col) -= sol * A[k];
            A[k] = 0.0;
            break;
          }
      }
  else
    for (RowNode *aux = Rows[rc]; aux != NULL; aux = aux->Prev)
      if ((col = aux->Column) == rc)
        if (d) {
          rhs(rc) = aux->Value * sol;
        } else {
          aux->Value = 1.0;
          rhs(rc) = sol;
        }
      else {
        aux->Value = 0.0;
        for (RowNode *node = Rows[col]; 1; node = node->Prev)
          if (node == NULL) {
            pmtop_error("SparseMatrix::EliminateRowCol () #3");
          } else if (node->Column == rc) {
            rhs(col) -= sol * node->Value;
            node->Value = 0.0;
            break;
          }
      }
}

void SparseMatrix::EliminateRowCol(int rc, int d) {
  int col;
  RowNode *aux, *node;

#ifdef NDEBUG
  if (rc >= size || rc < 0)
    pmtop_error("SparseMatrix::EliminateRowCol () #1");
#endif

  if (Rows == NULL)
    pmtop_error("SparseMatrix::EliminateRowCol () #2");

  for (aux = Rows[rc]; aux != NULL; aux = aux->Prev) {
    if ((col = aux->Column) == rc) {
      if (d == 0)
        aux->Value = 1e12;
    } else {
      aux->Value = 0.0;
      for (node = Rows[col]; 1; node = node->Prev)
        if (node == NULL) {
          pmtop_error("SparseMatrix::EliminateRowCol () #3");
        } else if (node->Column == rc) {
          node->Value = 0.0;
          break;
        }
    }
  }
}

void SparseMatrix::EliminateRowCol(int rc, SparseMatrix &Ae, int d) {
  int col;

  if (Rows) {
    RowNode *nd, *nd2;
    for (nd = Rows[rc]; nd != NULL; nd = nd->Prev) {
      if ((col = nd->Column) == rc) {
        if (d == 0) {
          Ae.Add(rc, rc, nd->Value - 1.0);
          nd->Value = 1.0;
        }
      } else {
        Ae.Add(rc, col, nd->Value);
        nd->Value = 0.0;
        for (nd2 = Rows[col]; 1; nd2 = nd2->Prev) {
          if (nd2 == NULL) {
            pmtop_error("SparseMatrix::EliminateRowCol");
          } else if (nd2->Column == rc) {
            Ae.Add(col, rc, nd2->Value);
            nd2->Value = 0.0;
            break;
          }
        }
      }
    }
  } else {
    pmtop_error("SparseMatrix::EliminateRowCol");
  }
}

void SparseMatrix::SetDiagIdentity() {
  for (int i = 0; i < size; i++)
    if (RowIndex[i + 1] == RowIndex[i] + 1 && fabs(A[RowIndex[i]]) < 1e-16)
      A[RowIndex[i]] = 1.0;
}

void SparseMatrix::EliminateZeroRows() {
  int i, j;
  double zero;

  for (i = 0; i < size; i++) {
    zero = 0.0;
    for (j = RowIndex[i]; j < RowIndex[i + 1]; j++)
      zero += fabs(A[j]);
    if (zero < 1e-12) {
      for (j = RowIndex[i]; j < RowIndex[i + 1]; j++)
        if (J[j] == i)
          A[j] = 1.0;
        else
          A[j] = 0.0;
    }
  }
}

// 该函数实现了稀疏矩阵的高斯-赛德尔前向迭代算法，用于求解线性方程组
// Ax = b。前向迭代意味着更新结果y 是按行逐一计算并立即使用最新的值。
// 利用稀疏矩阵的结构优化高斯-赛德尔迭代，减少计算成本。
// 遍历矩阵每一行，逐步更新解向量y
// 确保主对角线元素非零，如果为零则触发错误。
void SparseMatrix::Gauss_Seidel_forw(const FloatArray &x, FloatArray &y) const {
  int c, i, j, end, d, s = size, *Ip = RowIndex, *Jp = J;
  double sum, *Ap = A, *yp = y.givePointer();
  const double *xp = x.givePointer();

  if (A == NULL)
    pmtop_error("SparseMatrix::Gauss_Seidel_forw ()");

  j = Ip[0];
  for (i = 0; i < s; i++) {
    end = Ip[i + 1];
    sum = 0.0;
    d = -1;
    for (; j < end; j++)
      if ((c = Jp[j]) == i)
        d = j;
      else
        sum += Ap[j] * yp[c];

    if (d >= 0 && Ap[d] != 0.0)
      yp[i] = (xp[i] - sum) / Ap[d];
    else if (xp[i] == sum)
      yp[i] = sum;
    else
      pmtop_error("SparseMatrix::Gauss_Seidel_forw (...) #2");
  }
}

// 该函数实现了稀疏矩阵的高斯-赛德尔后向迭代算法，用于求解线性方程组
// Ax = b。后向迭代意味着按行从最后一行开始逆序更新解向量 y。
void SparseMatrix::Gauss_Seidel_back(const FloatArray &x, FloatArray &y) const {
  int i, j, beg, c, d;
  double sum, *Ap = A, *yp = y.givePointer();
  const double *xp = x.givePointer();
  int *Ip = RowIndex, *Jp = J;

  if (A == NULL)
    pmtop_error("SparseMatrix::Gauss_Seidel_back ()");

  j = Ip[size] - 1;
  for (i = size - 1; i >= 0; i--) {
    beg = Ip[i];
    sum = 0.;
    d = -1;
    for (; j >= beg; j--)
      if ((c = Jp[j]) == i)
        d = j;
      else
        sum += Ap[j] * yp[c];

    if (d >= 0 && Ap[d] != 0.0)
      yp[i] = (xp[i] - sum) / Ap[d];
    else if (xp[i] == sum)
      yp[i] = sum;
    else
      pmtop_error("SparseMatrix::Gauss_Seidel_back (...) #2");
  }
}

// 该函数计算稀疏矩阵的雅可比预处理因子（Jacobi
// Scaling）。雅可比预处理是一种对角线预处理，用于改善迭代解法（如高斯-赛德尔法、共轭梯度法）的收敛性能。
double SparseMatrix::GetJacobiScaling() const {
#ifdef NDEBUG
  if (A == NULL)
    pmtop_error("SparseMatrix::GetJacobiScaling()");
#endif

  double sc = 1.0;
  for (int i = 0; i < size; i++) {
    int d = -1;
    double norm = 0.0;
    for (int j = RowIndex[i]; j < RowIndex[i + 1]; j++) {
      if (J[j] == i)
        d = j;
      norm += fabs(A[j]);
    }
    if (d >= 0 && A[d] != 0.0) {
      double a = 1.8 * fabs(A[d]) / norm;
      if (a < sc)
        sc = a;
    } else
      pmtop_error("SparseMatrix::GetJacobiScaling() #2");
  }
  return sc;
}

// 该函数实现了稀疏矩阵的雅可比迭代法（Jacobi Iteration），用于求解线性方程组
// Ax=b。雅可比迭代法是一种基本的迭代方法，它将矩阵分解为对角线部分与其他部分，并逐步逼近解向量。
void SparseMatrix::Jacobi(const FloatArray &b, const FloatArray &x0,
                          FloatArray &x1, double sc) const {
#ifdef NDEBUG
  if (A == NULL)
    pmtop_error("SparseMatrix::Jacobi(...)");
#endif

  for (int i = 0; i < size; i++) {
    int d = -1;
    double sum = b(i);
    for (int j = RowIndex[i]; j < RowIndex[i + 1]; j++) {
      if (J[j] == i)
        d = j;
      else
        sum -= A[j] * x0(J[j]);
    }
    if (d >= 0 && A[d] != 0.0)
      x1(i) = sc * (sum / A[d]) + (1.0 - sc) * x0(i);
    else
      pmtop_error("SparseMatrix::Jacobi(...) #2");
  }
}

void SparseMatrix::Jacobi2(const FloatArray &b, const FloatArray &x0,
                           FloatArray &x1, double sc) const {
#ifdef NDEBUG
  if (A == NULL)
    pmtop_error("SparseMatrix::Jacobi2(...)");
#endif

  for (int i = 0; i < size; i++) {
    double resi = b(i), norm = 0.0;
    for (int j = RowIndex[i]; j < RowIndex[i + 1]; j++) {
      resi -= A[j] * x0(J[j]);
      norm += fabs(A[j]);
    }
    if (norm > 0.0)
      x1(i) = x0(i) + sc * resi / norm;
    else
      pmtop_error("SparseMatrix::Jacobi2(...) #2");
  }
}

void SparseMatrix::Finalize(int skip_zeros) {
#ifdef NDEBUG
  if (finalize==true) {
    pmtop_error("SparseMatrix::Finalize()");
  }
#endif

  int i, j, nr, nz;
  RowNode *aux;

  RowIndex = new int[size + 1];
  RowIndex[0] = 0;
  for (i = 1; i <= size; i++) {
    nr = 0;
    for (aux = Rows[i - 1]; aux != NULL; aux = aux->Prev)
      if (!skip_zeros || aux->Value != 0.0)
        nr++;
    RowIndex[i] = RowIndex[i - 1] + nr;
  }
  nz = RowIndex[size];
  J = new int[nz];
  A = new double[nz];
  for (j = i = 0; i < size; i++)
    for (aux = Rows[i]; aux != NULL; aux = aux->Prev)
      if (!skip_zeros || aux->Value != 0.0) {
        J[j] = aux->Column;
        A[j] = aux->Value;
        j++;
      }
  for (i = 0; i < size; i++) {
    RowNode *node_p = Rows[i];
    while (node_p != NULL) {
      aux = node_p;
      node_p = node_p->Prev;
      delete aux;
    }
  }
  delete[] Rows;
  Rows = NULL;
  finalize = true;
}

// 该函数用于提取稀疏矩阵的一个子矩阵（子块）。通过指定的行索引和列索引数组，从稀疏矩阵中提取对应的子矩阵，并存储到
// FloatMatrix 类型的对象 subm 中。
void SparseMatrix::GetSubMatrix(const IntArray &rows, const IntArray &cols,
                                FloatMatrix &subm) {
  int i, j, gi, gj, s, t;
  RowNode *aux;

  if (Rows == NULL)
    pmtop_error("SparseMatrix::GetSubMatrix(...) #0");

  for (i = 0; i < rows.giveSize(); i++) {
    if ((gi = rows[i]) < 0)
      gi = -1 - gi, s = -1;
    else
      s = 1;
#ifdef NDEBUG
    if (gi >= size)
      pmtop_error("SparseMatrix::GetSubMatrix(...) #1");
#endif
    for (j = 0; j < cols.giveSize(); j++) {
      if ((gj = cols[j]) < 0)
        gj = -1 - gj, t = -s;
      else
        t = s;
#ifdef NDEBUG
      if (gj >= width)
        pmtop_error("SparseMatrix::GetSubMatrix(...) #2");
#endif
      for (aux = Rows[gi]; 1; aux = aux->Prev)
        if (aux == NULL) {
          subm(i, j) = 0.0;
          break;
        } else if (aux->Column == gj) {
          subm(i, j) = (t < 0) ? (-aux->Value) : (aux->Value);
          break;
        }
    }
  }
}

void SparseMatrix::Set(const int i, const int j, const double A) {
  double a = A;
  int gi, gj, s, t;

  if ((gi = i) < 0)
    gi = -1 - gi, s = -1;
  else
    s = 1;
#ifdef NDEBUG
  if (gi >= size)
    pmtop_error("SparseMatrix::Set (...) #1");
#endif
  if ((gj = j) < 0)
    gj = -1 - gj, t = -s;
  else
    t = s;
#ifdef NDEBUG
  if (gj >= width)
    pmtop_error("SparseMatrix::Set (...) #2");
#endif
  if (t < 0)
    a = -a;
  _Set_(gi, gj, a);
}

void SparseMatrix::Add(const int i, const int j, const double A) {
  int gi, gj, s, t;
  double a = A;

  if ((gi = i) < 0)
    gi = -1 - gi, s = -1;
  else
    s = 1;
#ifdef NDEBUG
  if (gi >= size)
    pmtop_error("SparseMatrix::Add (...) #1");
#endif
  if ((gj = j) < 0)
    gj = -1 - gj, t = -s;
  else
    t = s;
#ifdef NDEBUG
  if (gj >= width)
    pmtop_error("SparseMatrix::Add (...) #2");
#endif
  if (t < 0)
    a = -a;
  _Add_(gi, gj, a);
}

// 该函数用于将一个子矩阵的值设置到稀疏矩阵的指定位置。它根据给定的行和列索引选择子矩阵中的元素，并将这些元素插入到稀疏矩阵的对应位置。通过
// skip_zeros 参数，用户可以选择是否跳过值为零的元素，从而避免不必要的操作。
//
void SparseMatrix::SetSubMatrix(const IntArray &rows, const IntArray &cols,
                                const FloatMatrix &subm, int skip_zeros) {
  int i, j, gi, gj, s, t;
  double a;

  for (i = 0; i < rows.giveSize(); i++) {
    if ((gi = rows[i]) < 0)
      gi = -1 - gi, s = -1;
    else
      s = 1;
#ifdef NDEBUG
    if (gi >= size)
      pmtop_error("SparseMatrix::SetSubMatrix(...) #1");
#endif
    for (j = 0; j < cols.giveSize(); j++) {
      a = subm(i, j);
      if (skip_zeros && a == 0.0)
        continue;
      if ((gj = cols[j]) < 0)
        gj = -1 - gj, t = -s;
      else
        t = s;
#ifdef NDEBUG
      if (gj >= width)
        pmtop_error("SparseMatrix::SetSubMatrix(...) #2");
#endif
      if (t < 0)
        a = -a;
      _Set_(gi, gj, a);
    }
  }
}

void SparseMatrix::SetSubMatrixTranspose(const IntArray &rows,
                                         const IntArray &cols,
                                         const FloatMatrix &subm,
                                         int skip_zeros) {
  int i, j, gi, gj, s, t;
  double a;

  for (i = 0; i < rows.giveSize(); i++) {
    if ((gi = rows[i]) < 0)
      gi = -1 - gi, s = -1;
    else
      s = 1;
#ifdef NDEBUG
    if (gi >= size)
      pmtop_error("SparseMatrix::SetSubMatrixTranspose (...) #1");
#endif
    for (j = 0; j < cols.giveSize(); j++) {
      a = subm(j, i);
      if (skip_zeros && a == 0.0)
        continue;
      if ((gj = cols[j]) < 0)
        gj = -1 - gj, t = -s;
      else
        t = s;
#ifdef NDEBUG
      if (gj >= width)
        pmtop_error("SparseMatrix::SetSubMatrixTranspose (...) #2");
#endif
      if (t < 0)
        a = -a;
      _Set_(gi, gj, a);
    }
  }
}

void SparseMatrix::AddSubMatrix(const IntArray &rows, const IntArray &cols,
                                const FloatMatrix &subm, int skip_zeros) {
  int i, j, gi, gj, s, t;
  double a;

  for (i = 0; i < rows.giveSize(); i++) {
    if ((gi = rows[i]) < 0)
      gi = -1 - gi, s = -1;
    else
      s = 1;
#ifdef NDEBUG
    if (gi >= size)
      pmtop_error("SparseMatrix::AddSubMatrix(...) #1");
#endif
    for (j = 0; j < cols.giveSize(); j++) {
      if ((gj = cols[j]) < 0)
        gj = -1 - gj, t = -s;
      else
        t = s;
#ifdef NDEBUG
      if (gj >= width)
        pmtop_error("SparseMatrix::AddSubMatrix(...) #2");
#endif
      a = subm(i, j);
      if (skip_zeros && a == 0.0) {
        // if the element is zero do not assemble it unless this breaks
        // the symmetric structure
        if (&rows != &cols || subm(j, i) == 0.0)
          continue;
      }
      if (t < 0)
        a = -a;
      _Add_(gi, gj, a);
    }
  }
}
// 该函数用于将一个新行（通过 srow 表示）加到稀疏矩阵的指定行 row
// 中。它根据输入的列索引 cols 和相应的行元素
// srow，将这些元素加到稀疏矩阵中。如果某些元素的值为零，则不执行任何操作。
//
void SparseMatrix::AddRow(const int row, const IntArray &cols,
                          const FloatArray &srow) {
  int j, gi, gj, s, t;
  double a;
#ifdef NDEBUG
  if (Rows == NULL)
    pmtop_error("SparseMatrix::AddRow(...) #0");
#endif

  if ((gi = row) < 0)
    gi = -1 - gi, s = -1;
  else
    s = 1;
#ifdef NDEBUG
  if (gi >= size)
    pmtop_error("SparseMatrix::AddRow(...) #1");
#endif
  for (j = 0; j < cols.giveSize(); j++) {
    if ((gj = cols[j]) < 0)
      gj = -1 - gj, t = -s;
    else
      t = s;
#ifdef NDEBUG
    if (gj >= width)
      pmtop_error("SparseMatrix::AddRow(...) #2");
#endif
    a = srow(j);
    if (a == 0.0)
      continue;
    if (t < 0)
      a = -a;
    _Add_(gi, gj, a);
  }
}

// 该函数用于将稀疏矩阵中指定行的所有元素按比例缩放。给定一个缩放因子
// scale，该函数会将指定行中的所有非零元素乘以 scale。
//
void SparseMatrix::ScaleRow(const int row, const double scale) {
  int i;

  if ((i = row) < 0)
    i = -1 - i;
  if (Rows != NULL) {
    RowNode *aux;

    for (aux = Rows[i]; aux != NULL; aux = aux->Prev)
      aux->Value *= scale;
  } else {
    int j, end = RowIndex[i + 1];

    for (j = RowIndex[i]; j < end; j++)
      A[j] *= scale;
  }
}

// 该函数是 SparseMatrix 类的加法赋值操作符重载。它的作用是将矩阵 B
// 中的所有元素加到当前矩阵（this）中。对于每个非零元素，它通过 _Add_()
// 方法将元素加到当前矩阵的对应位置。
//
SparseMatrix &SparseMatrix::operator+=(SparseMatrix &B) {
  int i;
  RowNode *aux;
#ifdef NDEBUG
  if (Rows == NULL || B.Rows == NULL)
    pmtop_error("SparseMatrix::operator+=(...) #0");

  if (size != B.size || width != B.width)
    pmtop_error("SparseMatrix::operator+=(...) #1");
#endif

  for (i = 0; i < size; i++) {
    for (aux = B.Rows[i]; aux != NULL; aux = aux->Prev) {
      _Add_(i, aux->Column, aux->Value);
    }
  }

  return (*this);
}
// 该函数是 SparseMatrix 类的赋值操作符重载，用于将一个常数值 a
// 赋给当前矩阵。矩阵的所有元素将被设置为该常数 a。
//
SparseMatrix &SparseMatrix::operator=(double a) {
  if (Rows == NULL)
    for (int i = 0, nnz = RowIndex[size]; i < nnz; i++)
      A[i] = a;
  else
    for (int i = 0; i < size; i++)
      for (RowNode *node_p = Rows[i]; node_p != NULL; node_p = node_p->Prev)
        node_p->Value = a;

  return (*this);
}

void SparseMatrix::Print(std::ostream &out, int _width) const {
  int i, j;
#ifdef NDEBUG
  if (A == NULL)
    pmtop_error("SparseMatrix::Print()");
#endif

  out << setiosflags(std::ios::scientific | std::ios::showpos);
  for (i = 0; i < size; i++) {
    out << "[row " << i << "]\n";
    for (j = RowIndex[i]; j < RowIndex[i + 1]; j++) {
      out << "(" << std::setw(3) << J[j] << "," << A[j] << ") ";
      if (!((j + 1 - RowIndex[i]) % _width))
        out << std::endl;
    }
    out << std::endl;
  }
  out << std::endl;
}

void SparseMatrix::PrintMatlab(std::ostream &out) const {
  int i, j;
  std::ios::fmtflags old_fmt = out.setf(std::ios::scientific);
  int old_prec = out.precision(14);

  for (i = 0; i < size; i++)
    for (j = RowIndex[i]; j < RowIndex[i + 1]; j++)
      out << i + 1 << " " << J[j] + 1 << " " << A[j] << std::endl;
  out.precision(old_prec);
  out.setf(old_fmt);
}

void SparseMatrix::PrintMM(std::ostream &out) const {
  int i, j;
  std::ios::fmtflags old_fmt = out.setf(std::ios::scientific);
  int old_prec = out.precision(14);

  out << "%%MatrixMarket matrix coordinate real general" << std::endl
      << "% Generated by AggieFEM" << std::endl;

  out << size << " " << width << " " << NumNonZeroElems() << std::endl;
  for (i = 0; i < size; i++)
    for (j = RowIndex[i]; j < RowIndex[i + 1]; j++)
      out << i + 1 << " " << J[j] + 1 << " " << A[j] << std::endl;
  out.precision(old_prec);
  out.setf(old_fmt);
}

void SparseMatrix::PrintCSR(std::ostream &out) const {
#ifdef NDEBUG
  if (A == NULL)
    pmtop_error("SparseMatrix::PrintCSR()");
#endif

  int i;
  std::ios::fmtflags old_fmt = out.setf(std::ios::scientific);
  int old_prec = out.precision(14);

  out << size << '\n'; // number of rows

  for (i = 0; i <= size; i++)
    out << RowIndex[i] + 1 << '\n';

  for (i = 0; i < RowIndex[size]; i++)
    out << J[i] + 1 << '\n';

  for (i = 0; i < RowIndex[size]; i++)
    out << A[i] << '\n';

  out.precision(old_prec);
  out.setf(old_fmt);
}

void SparseMatrix::PrintCSR2(std::ostream &out) const {
#ifdef NDEBUG
  if (A == NULL)
    pmtop_error("SparseMatrix::PrintCSR2()");
#endif

  int i;
  std::ios::fmtflags old_fmt = out.setf(std::ios::scientific);
  int old_prec = out.precision(14);

  out << size << '\n';  // number of rows
  out << width << '\n'; // number of columns

  for (i = 0; i <= size; i++)
    out << RowIndex[i] << '\n';

  for (i = 0; i < RowIndex[size]; i++)
    out << J[i] << '\n';

  for (i = 0; i < RowIndex[size]; i++)
    out << A[i] << '\n';

  out.precision(old_prec);
  out.setf(old_fmt);
}

double SparseMatrix::IsSymmetric() const {
#ifdef NDEBUG
  if (A == NULL)
    pmtop_error("SparseMatrix::IsSymmetric()");
#endif

  int i, j;
  double a, max;

  max = 0.0;
  for (i = 1; i < size; i++)
    for (j = RowIndex[i]; j < RowIndex[i + 1]; j++)
      if (J[j] < i) {
        a = fabs(A[j] - (*this)(J[j], i));
        if (max < a)
          max = a;
      }

  return max;
}

void SparseMatrix::Symmetrize() {
#ifdef NDEBUG
  if (A == NULL)
    pmtop_error("SparseMatrix::Symmetrize()");
#endif

  int i, j;
  for (i = 1; i < size; i++)
    for (j = RowIndex[i]; j < RowIndex[i + 1]; j++)
      if (J[j] < i) {
        A[j] += (*this)(J[j], i);
        A[j] *= 0.5;
        (*this)(J[j], i) = A[j];
      }
}

int SparseMatrix::NumNonZeroElems() const {
  if (A != NULL) //  matrix is finalized
    return RowIndex[size];
  pmtop_error("SparseMatrix::NumNonZeroElems");
  return -1;
}

int SparseMatrix::CountSmallElems(double tol) {
  int i, counter = 0;

  if (A) {
    int nz = RowIndex[size];
    double *Ap = A;

    for (i = 0; i < nz; i++)
      if (fabs(Ap[i]) < tol)
        counter++;
  } else {
    RowNode *aux;

    for (i = 0; i < size; i++)
      for (aux = Rows[i]; aux != NULL; aux = aux->Prev)
        if (fabs(aux->Value) < tol)
          counter++;
  }

  return counter;
}

#ifdef __MKL_MODULE
void SparseMatrix::covertToMKLCSR() {
#ifdef NDEBUG
  if (A == NULL)
    pmtop_error("SparseMatrix::covertToMKLCSR() #1");
  if (mkl_sparse_matrix != NULL)
    mkl_sparse_destroy(mkl_sparse_matrix);
#endif

  sparse_status_t status;
  status = mkl_sparse_d_create_csr(&mkl_sparse_matrix, SPARSE_INDEX_BASE_ZERO,
                                   size, width, RowIndex, RowIndex + 1, J, A);
  if (status != SPARSE_STATUS_SUCCESS)
    pmtop_error("SparseMatrix::covertToMKLCSR() #2");
  Use_MKL = true;
  mkl_descr.type = SPARSE_MATRIX_TYPE_GENERAL;
  mkl_descr.mode = SPARSE_FILL_MODE_FULL;
  mkl_descr.diag = SPARSE_DIAG_NON_UNIT;
}

FloatMatrix SparseMatrix::Mult(const FloatMatrix &B, double alpha,
                               double beta) {
  FloatMatrix C(size, B.giveNumberOfColumns());
  if (Use_MKL == false)
    pmtop_error("SparseMatrix::Mult(...) #1");
  mkl_sparse_d_mm(SPARSE_OPERATION_NON_TRANSPOSE, alpha, mkl_sparse_matrix,
                  mkl_descr, SPARSE_LAYOUT_COLUMN_MAJOR, B.givePointer(),
                  B.giveNumberOfColumns(), size, beta, C.givePointer(), size);
  return C;
}

FloatMatrix SparseMatrix::RAP(const FloatMatrix &R) {
  FloatMatrix C(R.giveNumberOfColumns(), R.giveNumberOfColumns());
#ifdef NDEBUG
  if (C.giveNumberOfColumns() != size)
    pmtop_error("SparseMatrix::RAP(...) #1");
  if (Use_MKL == false)
    pmtop_error("SparseMatrix::RAP(...) #2");
#endif
  FloatMatrix P = R.transpose();
  mkl_sparse_d_mm(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, mkl_sparse_matrix,
                  mkl_descr, SPARSE_LAYOUT_COLUMN_MAJOR, P.givePointer(),
                  P.giveNumberOfColumns(), size, 0.0, C.givePointer(),
                  C.giveNumberOfColumns());
  return R * C;
}

#endif

// 该函数接受一个稀疏矩阵 S 和一个函数指针 f，然后将矩阵 S
// 中的每个非零元素应用函数 f，即将每个非零元素的值替换为 f(旧值)。
void SparseMatrixFunction(SparseMatrix &S, double (*f)(double)) {
  int n = S.NumNonZeroElems();
  double *s = S.GetData();

  for (int i = 0; i < n; i++)
    s[i] = f(s[i]);
}

SparseMatrix::~SparseMatrix() {
  if (Rows != NULL) {
    for (int i = 0; i < size; i++) {
      RowNode *aux, *node_p = Rows[i];
      while (node_p != NULL) {
        aux = node_p;
        node_p = node_p->Prev;
        delete aux;
      }
    }
    delete[] Rows;
  }
  if (A != NULL) {
    delete[] RowIndex;
    delete[] J;
    delete[] A;
  }
  if (Use_MKL == true) {
    mkl_sparse_destroy(mkl_sparse_matrix);
  }
}

SparseMatrix *Transpose(SparseMatrix &A) {
  int i, j, end;
  int m, n, nnz, *A_i, *A_j, *At_i, *At_j;
  double *A_data, *At_data;

  m = A.Size();  // number of rows of A
  n = A.Width(); // number of columns of A
  nnz = A.NumNonZeroElems();
  A_i = A.GetI();
  A_j = A.GetJ();
  A_data = A.GetData();

  At_i = new int[n + 1];
  At_j = new int[nnz];
  At_data = new double[nnz];

  for (i = 0; i <= n; i++)
    At_i[i] = 0;
  for (i = 0; i < nnz; i++)
    At_i[A_j[i] + 1]++;
  for (i = 1; i < n; i++)
    At_i[i + 1] += At_i[i];

  for (i = j = 0; i < m; i++) {
    end = A_i[i + 1];
    for (; j < end; j++) {
      At_j[At_i[A_j[j]]] = i;
      At_data[At_i[A_j[j]]] = A_data[j];
      At_i[A_j[j]]++;
    }
  }

  for (i = n; i > 0; i--)
    At_i[i] = At_i[i - 1];
  At_i[0] = 0;

  return new SparseMatrix(At_i, At_j, At_data, n, m);
}

SparseMatrix *Mult(SparseMatrix &A, SparseMatrix &B, SparseMatrix *OAB) {

  int nrowsA, ncolsA, nrowsB, ncolsB;
  int *A_i, *A_j, *B_i, *B_j, *C_i, *C_j, *B_marker;
  double *A_data, *B_data, *C_data;
  int ia, ib, ic, ja, jb, num_nonzeros;
  int row_start, counter;
  double a_entry, b_entry;
  SparseMatrix *C;

  nrowsA = A.Size();
  ncolsA = A.Width();
  nrowsB = B.Size();
  ncolsB = B.Width();
#ifdef NDEBUG
  if (ncolsA != nrowsB)
    pmtop_error("Sparse matrix multiplication, Mult (...) #1");
#endif
  A_i = A.GetI();
  A_j = A.GetJ();
  A_data = A.GetData();
  B_i = B.GetI();
  B_j = B.GetJ();
  B_data = B.GetData();

  B_marker = new int[ncolsB];

  for (ib = 0; ib < ncolsB; ib++)
    B_marker[ib] = -1;

  if (OAB == NULL) {
    C_i = new int[nrowsA + 1];

    C_i[0] = num_nonzeros = 0;
    for (ic = 0; ic < nrowsA; ic++) {
      for (ia = A_i[ic]; ia < A_i[ic + 1]; ia++) {
        ja = A_j[ia];
        for (ib = B_i[ja]; ib < B_i[ja + 1]; ib++) {
          jb = B_j[ib];
          if (B_marker[jb] != ic) {
            B_marker[jb] = ic;
            num_nonzeros++;
          }
        }
      }
      C_i[ic + 1] = num_nonzeros;
    }

    C_j = new int[num_nonzeros];
    C_data = new double[num_nonzeros];

    C = new SparseMatrix(C_i, C_j, C_data, nrowsA, ncolsB);

    for (ib = 0; ib < ncolsB; ib++)
      B_marker[ib] = -1;
  } else {
    C = OAB;
#ifdef NDEBUG
    if (nrowsA != C->Size() || ncolsB != C->Width())
      pmtop_error("Sparse matrix multiplication, Mult (...) #2");
#endif
    C_i = C->GetI();
    C_j = C->GetJ();
    C_data = C->GetData();
  }

  counter = 0;
  for (ic = 0; ic < nrowsA; ic++) {
    // row_start = C_i[ic];
    row_start = counter;
    for (ia = A_i[ic]; ia < A_i[ic + 1]; ia++) {
      ja = A_j[ia];
      a_entry = A_data[ia];
      for (ib = B_i[ja]; ib < B_i[ja + 1]; ib++) {
        jb = B_j[ib];
        b_entry = B_data[ib];
        if (B_marker[jb] < row_start) {
          B_marker[jb] = counter;
          if (OAB == NULL)
            C_j[counter] = jb;
          C_data[counter] = a_entry * b_entry;
          counter++;
        } else
          C_data[B_marker[jb]] += a_entry * b_entry;
      }
    }
  }
#ifdef NDEBUG
  if (OAB != NULL && counter != OAB->NumNonZeroElems())
    pmtop_error("Sparse matrix multiplication, Mult (...) #3");
#endif

  delete[] B_marker;

  return C;
}

SparseMatrix *RAP(SparseMatrix &A, SparseMatrix &R, SparseMatrix *ORAP) {
  SparseMatrix *P = Transpose(R);
  SparseMatrix *AP = Mult(A, *P);
  delete P;
  SparseMatrix *_RAP = Mult(R, *AP, ORAP);
  delete AP;
  return _RAP;
}

SparseMatrix *Mult_AtDA(SparseMatrix &A, FloatArray &D, SparseMatrix *OAtDA) {
  int i, At_nnz, *At_j;
  double *At_data;

  SparseMatrix *At = Transpose(A);
  At_nnz = At->NumNonZeroElems();
  At_j = At->GetJ();
  At_data = At->GetData();
  for (i = 0; i < At_nnz; i++)
    At_data[i] *= D(At_j[i]);
  SparseMatrix *AtDA = Mult(*At, A, OAtDA);
  delete At;
  return AtDA;
}

void SparseMatrix::SortColumnIndices() {
#ifdef NDEBUG
 if (finalize==false) {
   std::cerr<<"Sparse matrix does not finalize!\n";
 }
#endif
  if (isSorted)
  {
    return;
  }
  Array<Pair<int,double> > row;
  for (int j = 0, i = 0; i < size; i++) {
    int end = RowIndex[i+1];
    row.SetSize(end - j);
    for (int k = 0; k < row.Size(); k++)
    {
      row[k].one = J[j+k];
      row[k].two = A[j+k];
    }
    row.Sort();
    for (int k = 0; k < row.Size(); k++, j++)
    {
      J[j] = row[k].one;
      A[j] = row[k].two;
    }
  }
  isSorted = true;
}

SparseMatrixComplex::SparseMatrixComplex(const SparseMatrix &Real, std::string type) {

  assert(Real.isFinalized()==true);
  height=Real.Height();
  width=Real.Width();
  RowIndex.resize(width+1);
  J.resize(Real.GetI()[width]);
  A.resize(Real.GetI()[width]);
  for (int i = 0; i < width+1; i++) {
    RowIndex[i] = Real.GetI()[i];
  }
  for (int j=0;j<Real.GetI()[width];j++) {
    J[j] = Real.GetJ()[j];
    if (type == "Real") {
      A.Set(j,Real.GetData()[j],0);
    }
    else if (type == "Imag") {
      A.Set(j,0,Real.GetData()[j]);
    }
    else {
      std::cerr<<"Unknown type "<<type<<"\n";
    }
  }
  covertToMKLCSR();
}

SparseMatrixComplex::SparseMatrixComplex(const SparseMatrix &Real, const SparseMatrix &Imag) {
  assert(Real.isFinalized()==true);
  height=Real.Height();
  width=Real.Width();
#ifdef NDEBUG
  for (int i=0;i<width;i++) {
    assert(Real.GetI()[i]==Imag.GetI()[i]);
  }
#endif
  RowIndex.resize(width+1);
  J.resize(Real.GetI()[width]);
  A.resize(Real.GetI()[width]);
  for (int i = 0; i < width+1; i++) {
    RowIndex[i] = Real.GetI()[i];
  }
  for (int j=0;j<Real.GetI()[width];j++) {
    J[j] = Real.GetJ()[j];
    A.Set(j,Real.GetData()[j],Imag.GetData()[j]);
  }
  covertToMKLCSR();
}


void SparseMatrixComplex::covertToMKLCSR() {
  if (A.empty())
    pmtop_error("SparseMatrix::covertToMKLCSR() #1");
  if (mkl_sparse_matrix != NULL)
    mkl_sparse_destroy(mkl_sparse_matrix);
  sparse_status_t status;
  status = mkl_sparse_z_create_csr(&mkl_sparse_matrix, SPARSE_INDEX_BASE_ZERO,
                                   width, width, RowIndex.givePointer(), RowIndex.givePointer() + 1,
                                   J.givePointer(), A.Data());
  if (status != SPARSE_STATUS_SUCCESS)
    pmtop_error("SparseMatrix::covertToMKLCSR() #2");
  Use_MKL = true;
  mkl_descr.type = SPARSE_MATRIX_TYPE_GENERAL;
  mkl_descr.mode = SPARSE_FILL_MODE_FULL;
  mkl_descr.diag = SPARSE_DIAG_NON_UNIT;
}

void SparseMatrixComplex::EliminateElement(const IntArray & index,bool diagVal ) {
  for (int i=0;i<index.giveSize();i++) {
    ZeroRowColAndKeepDiag(index(i),diagVal);
  }
}




const ComplexScalar &SparseMatrixComplex::operator()(int i, int j) const {
  int k, end;
  static const ComplexScalar zero(0.0,0.0);

#ifdef NDEBUG
  if (i >= height || i < 0 || j < 0)
    pmtop_error("SparseMatrix::operator() const #1");
#endif

  if (A.empty())
    pmtop_error("SparseMatrix::operator() const #2");
  end = RowIndex[i + 1];
  for (k = RowIndex[i]; k < end; k++)
    if (J[k] == j)
      return A[k];

  return zero;
}

const ComplexScalar &SparseMatrixComplex::Elem(int i, int j) const {
  return operator()(i, j);
}

void SparseMatrixComplex::Mult(const ComplexVector &x, ComplexVector &y) const {
  ComplexScalar one(1.0, 0.0);
  AddMult(x, y, one);
}

void SparseMatrixComplex::Mult(const ComplexMatrix &x, ComplexMatrix &y) const {
  ComplexScalar one(1.0, 0.0);
  AddMult(x, y, one);
}

void SparseMatrixComplex::Add(const ComplexScalar &a,SparseMatrixComplex &B, SparseMatrixComplex &C) {
  mkl_sparse_destroy(C.getMKLMatrix());
  mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE,this->getMKLMatrix(),a.givePointer(),B.getMKLMatrix(),&C.getMKLMatrix());
  C.Use_MKL = true;
  MKL_INT *row_start = nullptr;
  MKL_INT *row_end   = nullptr;
  MKL_INT *col_indx  = nullptr;
  MKL_Complex16 *vals = nullptr;
  sparse_index_base_t index;

  mkl_sparse_z_export_csr(
      C.getMKLMatrix(),
      &index,
      &C.height,
      &C.width,
      &row_start,
      &row_end,
      &col_indx,
      &vals
  );
  C.GetI().resize(C.height+1);
  C.GetJ().resize(row_start[C.height]);
  C.GetA().resize(row_start[C.height]);

  std::memcpy(C.GetI().givePointer(), row_start, sizeof(MKL_INT) * (C.height + 1));
  std::memcpy(C.GetJ().givePointer(), col_indx, sizeof(MKL_INT) * row_start[C.height]);
  std::memcpy(C.GetA().Data(), vals, sizeof(MKL_Complex16) * row_start[C.height]);

}

void SparseMatrixComplex::Augment(const ComplexVector &vec, const ComplexScalar &diagVal) {
  // 1. 检查向量维度
    // 假设 ComplexVector 有 Size() 或类似的长度获取方法，这里假设它应该等于 height
    // if (vec.Size() != height) throw std::runtime_error("Vector size mismatch");

    // 2. 如果之前使用了 MKL，结构改变导致句柄失效，必须销毁
    if (Use_MKL && mkl_sparse_matrix != nullptr) {
        mkl_sparse_destroy(mkl_sparse_matrix);
        mkl_sparse_matrix = nullptr;
        Use_MKL = false; // 需要重新调用 covertToMKLCSR 才会启用
    }

    // 3. 计算新的非零元总数 (NNZ)
    int old_nnz = RowIndex[height];
    int new_nnz = old_nnz;

    // 统计向量造成的非零元增加
    // 假设 ComplexScalar 可以与 0 比较，或者你需要自己定义容差
    double tol = 1e-16;
    int vec_nz_count = 0;

    // 注意：这里需要访问 vec 的元素，假设 vec 支持下标访问 vec[i]
    for (int i = 0; i < height; ++i) {
        if (std::abs(vec[i].real()) > tol || std::abs(vec[i].imag()) > tol) {
            vec_nz_count++;
        }
    }

    // 新增的列贡献 vec_nz_count 个非零元
    // 新增的行贡献 vec_nz_count 个非零元 (对称位置)
    // 新增的对角线贡献 1 个非零元 (即使是0通常也保留结构，或者你可以判断)
    new_nnz += vec_nz_count * 2 + 1;

    // 4. 分配新的 CSR 数组
    int new_height = height + 1;
    int new_width = width + 1;

    IntArray new_RowIndex;
    IntArray new_J;
    ComplexVector new_A;

    // 假设 IntArray 和 ComplexVector 有 Resize 方法
    // 如果没有，你需要构造新的临时对象
    new_RowIndex.resize(new_height + 1);
    new_J.resize(new_nnz);
    new_A.resize(new_nnz);

    // 5. 填充数据
    int current_pos = 0;
    new_RowIndex[0] = 0;

    // --- 步骤 A: 处理原有的行 (0 到 height-1) ---
    for (int i = 0; i < height; ++i) {
        // A.1 复制当前行原有的列元素
        int row_start = RowIndex[i];
        int row_end = RowIndex[i + 1];

        for (int k = row_start; k < row_end; ++k) {
            new_J[current_pos] = J[k];
            new_A.Set(current_pos,A[k]);
            current_pos++;
        }

        // A.2 在当前行末尾追加新的一列 (即 vec[i])
        // CSR要求列索引递增，因为新列索引是 height (原最大索引+1)，所以直接追加即可保持有序
        const ComplexScalar& val = vec[i];
        if (std::abs(val.real()) > tol || std::abs(val.imag()) > tol) {
            new_J[current_pos] = width; // 新的列索引就是原来的 width
            new_A.Set(current_pos,val);
            current_pos++;
        }

        // 更新行指针
        new_RowIndex[i + 1] = current_pos;
    }

    // --- 步骤 B: 处理新增的最后一行 (行索引为 height) ---
    // B.1 填入新行的前 n 个元素 (即 vec 的转置部分)
    // 我们按照列索引 0 到 height-1 遍历，这样天然满足 CSR 列排序要求
    for (int j = 0; j < height; ++j) {
        const ComplexScalar& val = vec[j];
        if (std::abs(val.real()) > tol || std::abs(val.imag()) > tol) {
            new_J[current_pos] = j;
            new_A.Set(current_pos,val); // 如果是共轭对称，这里可能需要 std::conj(val)
            current_pos++;
        }
    }

    // B.2 填入右下角的对角线元素
    new_J[current_pos] = width;
    new_A.Set(current_pos,diagVal);
    current_pos++;

    // 结束最后一行
    new_RowIndex[new_height] = current_pos;

    // 6. 替换成员变量
    // 这里利用 IntArray/ComplexVector 的移动赋值或 Swap，减少拷贝开销
    // 假设你的类支持 operator= 或者你需要手动 Swap
    this->RowIndex = new_RowIndex; // 可能会触发拷贝，最好实现 Move Semantics
    this->J = new_J;
    this->A = new_A;

    this->width = new_width;
    this->height = new_height;

}

void SparseMatrixComplex::AddMult(const ComplexVector &x, ComplexVector &y, const ComplexScalar& a) const {
#ifdef NDEBUG
  if ((width != x.Size()))
    pmtop_error("SparseMatrix::AddMult() #1");
#endif

  ComplexScalar zero(0.0, 0.0);

#ifdef __MKL_MODULE
  if (Use_MKL==true)
    mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, a.givePointer(), mkl_sparse_matrix,
                    mkl_descr, x.Data(), zero, y.Data());
  else
    std::cerr<<"SparseMatrix::AddMult() #2"<<std::endl;
#endif
}

void SparseMatrixComplex::AddMult(const ComplexMatrix &x, ComplexMatrix &y, const ComplexScalar& a) const {
#ifdef NDEBUG
  if ((width != x.Size()) )
    pmtop_error("SparseMatrix::AddMult() #1");
#endif

  ComplexScalar zero(0.0, 0.0);

#ifdef __MKL_MODULE
  if (Use_MKL==true)
    mkl_sparse_z_mm(SPARSE_OPERATION_NON_TRANSPOSE,a.givePointer(),mkl_sparse_matrix,mkl_descr,
      SPARSE_LAYOUT_COLUMN_MAJOR,x.Data(),x.Cols(),this->Rows(),zero,y.Data(),this->Rows());
  else
    std::cerr<<"SparseMatrix::AddMult() #2"<<std::endl;
#endif

}


void SparseMatrixComplex::AddMultTranspose(const ComplexVector &x, ComplexVector &y, const ComplexScalar &a) const {
#ifdef NDEBUG
  if ((width != x.Size()) )
    pmtop_error("SparseMatrix::AddMult() #1");
#endif

  ComplexScalar zero(0.0, 0.0);

#ifdef __MKL_MODULE
  if (Use_MKL==true)
    mkl_sparse_z_mv(SPARSE_OPERATION_TRANSPOSE, a.givePointer(), mkl_sparse_matrix,
                    mkl_descr, x.Data(), zero, y.Data());
  else
    std::cerr<<"SparseMatrix::AddMult() #2"<<std::endl;
#endif
}

void SparseMatrixComplex::AddMultConjugateTranspose(const ComplexVector &x, ComplexVector &y, const ComplexScalar &a) const {
#ifdef NDEBUG
  if ((width != x.Size()) )
    pmtop_error("SparseMatrix::AddMult() #1");
#endif

  ComplexScalar zero(0.0, 0.0);

#ifdef __MKL_MODULE
  if (Use_MKL==true)
    mkl_sparse_z_mv(SPARSE_OPERATION_CONJUGATE_TRANSPOSE, a.givePointer(), mkl_sparse_matrix,
                    mkl_descr, x.Data(), zero, y.Data());
  else
    std::cerr<<"SparseMatrix::AddMult() #2"<<std::endl;
#endif
}


void SparseMatrixComplex::SortColumnIndices() {
  if (isSorted)
  {
    return;
  }
  Array<Pair<int,MKL_Complex16> > row;
  for (int j = 0, i = 0; i < Rows(); i++) {
    int end = RowIndex[i+1];
    row.SetSize(end - j);
    for (int k = 0; k < row.Size(); k++)
    {
      row[k].one = J[j+k];
      row[k].two = A[j+k];
    }
    row.Sort();
    for (int k = 0; k < row.Size(); k++, j++)
    {
      J[j] = row[k].one;
      A.Set(j,row[k].two);
    }
  }
  isSorted = true;
}