// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#pragma once
#include "floatarray.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "matrix.hpp"
#include <cmath>
#include <complex.h>
#include <mkl_spblas.h>

#include "ComplexMatrix.h"
#include "ComplexVector.hpp"

#ifdef __MKL_MODULE
#include <mkl.h>
#endif

class RowNode {
public:
  RowNode *Prev;
  int Column;
  double Value;
};

class SparseMatrix : public Matrix {
private:
  /** Arrays for the connectivity information in the CSR storage.
       RowIndex is of size "size+1", J is of size the number of nonzero entries
       in the Sparse matrix (actually stored RowIndex[size]) */
  int *RowIndex, *J, width;

  /// The nonzero entries in the Sparse matrix with size RowIndex[size].
  double *A;

  RowNode **Rows;

  inline double &SearchRow(const int row, const int col);
  inline void _Add_(const int row, const int col, const double a) {
    SearchRow(row, col) += a;
  }
  inline void _Set_(const int row, const int col, const double a) {
    SearchRow(row, col) = a;
  }

  bool externalData;
  bool Use_MKL = false;
  bool finalize = false;
  bool isSorted=false;

#ifdef __MKL_MODULE
  sparse_matrix_t mkl_sparse_matrix;
  matrix_descr mkl_descr;
#endif

public:
  /// Creates sparse matrix.
  SparseMatrix(int nrows, int ncols = 0);

  SparseMatrix(int *i, int *j, double *data, int m, int n,
               bool ExternalData = false)
      : Matrix(m), RowIndex(i), J(j), width(n), A(data),
        externalData(ExternalData) {
    Rows = NULL;
    mkl_sparse_matrix = NULL;
    finalize=externalData;
  }

  bool isFinalized() const { return finalize; }
  /// Return the array RowIndex
  inline int *GetI() const { return RowIndex; }
  /// Return the array J
  inline int *GetJ() const { return J; }
  /// Return element data
  inline double *GetData() const { return A; }
  /// Return the number of columns
  inline int Width() const { return width; }
  /// Returns the number of elements in row i
  int RowSize(int i);

  /// Sort the column indices corresponding to each row.
  void SortColumnIndices();

  /// Returns reference to a_{ij}.  Index i, j = 0 .. size-1
  virtual double &Elem(int i, int j);

  /// Returns constant reference to a_{ij}.  Index i, j = 0 .. size-1
  virtual const double &Elem(int i, int j) const;

  /// Returns reference to A[i][j].  Index i, j = 0 .. size-1
  double &operator()(int i, int j);

  /// Returns reference to A[i][j].  Index i, j = 0 .. size-1
  const double &operator()(int i, int j) const;

  /// Matrix vector multiplication.
  virtual void Mult(const FloatArray &x, FloatArray &y) const;

  virtual void Mult(const double* x, double* y) const;

  /// y += A * x (default)  or  y += a * A * x
  void AddMult(const FloatArray &x, FloatArray &y, const double a = 1.0) const;

  /// Multiply a vector with the transposed matrix. y = At * x
  void MultTranspose(const FloatArray &x, FloatArray &y) const;

  /// y += At * x (default)  or  y += a * At * x
  void AddMultTranspose(const FloatArray &x, FloatArray &y,
                        const double a = 1.0) const;
  void PartMult(const IntArray &rows, const FloatArray &x, FloatArray &y);

  /// Compute y^t A x
  double InnerProduct(const FloatArray &x, const FloatArray &y) const;

  /// Returns a pointer to approximation of the matrix inverse.
  virtual MatrixInverse *Inverse() const;

  /// Eliminates a column from the transpose matrix.
  void EliminateRow(int row, const double sol, FloatArray &rhs);
  void EliminateRow(int row);
  void EliminateCol(int col);
  /// Eliminate all columns 'i' for which cols[i] != 0
  void EliminateCols(IntArray &cols, FloatArray *x = NULL,
                     FloatArray *b = NULL);
  /** Eliminates the column 'rc' to the 'rhs', deletes the row 'rc' and
       replaces the element (rc,rc) with 1.0; assumes that element (i,rc)
       is assembled if and only if the element (rc,i) is assembled.
       If d != 0 then the element (rc,rc) remains the same. */
  void EliminateRowCol(int rc, const double sol, FloatArray &rhs, int d = 0);
  void EliminateRowCol(int rc, int d = 0);
  // Same as above + save the eliminated entries in Ae so that
  // (*this) + Ae is the original matrix
  void EliminateRowCol(int rc, SparseMatrix &Ae, int d = 0);

  /// If a row contains only one diag entry of zero, set it to 1.
  void SetDiagIdentity();
  /// If a row contains only zeros, set its diagonal to 1.
  void EliminateZeroRows();

  // Gauss-Seidel forward and backward iterations over a vector x.
  void Gauss_Seidel_forw(const FloatArray &x, FloatArray &y) const;
  void Gauss_Seidel_back(const FloatArray &x, FloatArray &y) const;

  /// Determine appropriate scaling for Jacobi iteration
  double GetJacobiScaling() const;

  /** One scaled Jacobi iteration for the system A x = b.
       x1 = x0 + sc D^{-1} (b - A x0)  where D is the diag of A. */
  void Jacobi(const FloatArray &b, const FloatArray &x0, FloatArray &x1,
              double sc) const;

  /** x1 = x0 + sc D^{-1} (b - A x0) where \f$ D_{ii} = \sum_j |A_{ij}| \f$. */
  void Jacobi2(const FloatArray &b, const FloatArray &x0, FloatArray &x1,
               double sc = 1.0) const;

  /** Finalize the matrix initialization. The function should be called
       only once, after the matrix has been initialized. It's densenning
       the J and A arrays (by getting rid of -1s in array J). */
  virtual void Finalize(int skip_zeros = 1);

  int Finalized() { return (A != NULL); }

  void GetSubMatrix(const IntArray &rows, const IntArray &cols,
                    FloatMatrix &subm);

  void Set(const int i, const int j, const double a);
  void Add(const int i, const int j, const double a);

  void SetSubMatrix(const IntArray &rows, const IntArray &cols,
                    const FloatMatrix &subm, int skip_zeros = 1);

  void SetSubMatrixTranspose(const IntArray &rows, const IntArray &cols,
                             const FloatMatrix &subm, int skip_zeros = 1);

  void AddSubMatrix(const IntArray &rows, const IntArray &cols,
                    const FloatMatrix &subm, int skip_zeros = 1);

  void AddRow(const int row, const IntArray &cols, const FloatArray &srow);

  void ScaleRow(const int row, const double scale);

  /** Add a sparse matrix to "*this" sparse marix
       Both marices should not be finilized */
  SparseMatrix &operator+=(SparseMatrix &B);

  SparseMatrix &operator=(double a);

  /// Prints matrix to stream out.
  void Print(std::ostream &out = std::cout, int width = 4) const;

  /// Prints matrix in matlab format.
  void PrintMatlab(std::ostream &out = std::cout) const;

  /// Prints matrix in Matrix Market sparse format.
  void PrintMM(std::ostream &out = std::cout) const;

  /// Prints matrix to stream out in hypre_CSRMatrix format.
  void PrintCSR(std::ostream &out) const;

  /// Prints a sparse matrix to stream out in CSR format.
  void PrintCSR2(std::ostream &out) const;

  /// Walks the sparse matrix
  int Walk(int &i, int &j, double &a);

  /// Returns max_{i,j} |(i,j)-(j,i)| for a finalized matrix
  double IsSymmetric() const;

  /// (*this) = 1/2 ((*this) + (*this)^t)
  void Symmetrize();

  /// Returns the number of the nonzero elements in the matrix
  int NumNonZeroElems() const;

  /// Count the number of entries with |a_ij| < tol
  int CountSmallElems(double tol);

  /// Call this if data has been stolen.
  void LoseData() {
    RowIndex = 0;
    J = 0;
    A = 0;
  }
#ifdef __MKL_MODULE
  void covertToMKLCSR();
  bool isUseMKL() { return Use_MKL; };
  const sparse_matrix_t *getMKLMatrix() const { return &mkl_sparse_matrix; };
  FloatMatrix Mult(const FloatMatrix &B, double alpha = 1.0, double beta = 0.0);
  FloatMatrix RAP(const FloatMatrix &R);
#endif
  /// Destroys sparse matrix.
  virtual ~SparseMatrix();
};

/// Applies f() to each element of the matrix (after it is finalized).
void SparseMatrixFunction(SparseMatrix &S, double (*f)(double));

/// Transpose of a sparse matrix. A must be finalized.
SparseMatrix *Transpose(SparseMatrix &A);

/** Matrix product A.B.
    If OAB is not NULL, we assume it has the structure
    of A.B and store the result in OAB.
    If OAB is NULL, we create a new SparseMatrix to store
    the result and return a pointer to it.
    All matrices must be finalized.  */
SparseMatrix *Mult(SparseMatrix &A, SparseMatrix &B, SparseMatrix *OAB = NULL);

/** RAP matrix product. ORAP is like OAB above.
    All matrices must be finalized.  */
SparseMatrix *RAP(SparseMatrix &A, SparseMatrix &R, SparseMatrix *ORAP = NULL);

/** Matrix multiplication A^t D A.
    All matrices must be finalized.  */
SparseMatrix *Mult_AtDA(SparseMatrix &A, FloatArray &D,
                        SparseMatrix *OAtDA = NULL);

// 该函数用于访问稀疏矩阵中指定行(row)和列(col)位置的元素。如果该位置的元素不存在，函数会在需要时创建一个值为
// 0.0 的新元素并返回其引用。
// 该函数支持两种稀疏矩阵的内部存储方式：
// 链表存储：当 Rows 不为 NULL 时，稀疏矩阵的每一行用一个链表表示。
// 压缩行存储（CSR，CompressedSparse Row）：当 Rows 为 NULL 时，矩阵使用 CSR
// 格式存储。

inline double &SparseMatrix::SearchRow(const int row, const int col) {
  if (Rows) {
    RowNode *node_p;

    for (node_p = Rows[row]; 1; node_p = node_p->Prev)
      if (node_p == NULL) {
        node_p = new RowNode;
        node_p->Prev = Rows[row];
        node_p->Column = col;
        node_p->Value = 0.0;
        Rows[row] = node_p;
        break;
      } else if (node_p->Column == col) {
        break;
      }
    return node_p->Value;
  } else {
    int *Ip = RowIndex + row, *Jp = J;
    for (int k = Ip[0], end = Ip[1]; k < end; k++)
      if (Jp[k] == col)
        return A[k];
    std::cerr << "SparseMatrix::SearchRow (...)" << std::endl;
  }
  return A[0];
}


#ifdef __MKL_MODULE
class SparseMatrixComplex {
private:
  int  width;
  int height;
  IntArray RowIndex, J;
  ComplexVector A;

  sparse_matrix_t mkl_sparse_matrix = NULL;
  bool Use_MKL;
  matrix_descr mkl_descr;
  bool isSorted = false;


public:

  SparseMatrixComplex()
      : width(0),
        height(0),
        RowIndex(),
        J(),
        A(),
        mkl_sparse_matrix(NULL),
        Use_MKL(false)            // 默认不开 MKL
  {
    // MKL 描述符初始化为默认值
    mkl_descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    mkl_descr.mode = SPARSE_FILL_MODE_FULL;
    mkl_descr.diag = SPARSE_DIAG_NON_UNIT;
  }
  SparseMatrixComplex(const SparseMatrix &Real,const SparseMatrix &Imag);
  SparseMatrixComplex(const SparseMatrix &Real, std::string type);

  SparseMatrixComplex(const SparseMatrixComplex&) = delete;
  SparseMatrixComplex& operator=(const SparseMatrixComplex&) = delete;

  int Rows() const{ return width; };
  int Cols() const{ return height; };

  // 移动构造
  SparseMatrixComplex(SparseMatrixComplex&& other) noexcept
      : width(other.width), height(other.height),
        RowIndex(std::move(other.RowIndex)),
        J(std::move(other.J)),
        A(std::move(other.A)),
        mkl_sparse_matrix(other.mkl_sparse_matrix),
        Use_MKL(other.Use_MKL),
        mkl_descr(other.mkl_descr)
  {
    other.mkl_sparse_matrix = nullptr;
    other.Use_MKL = false;
    other.width = other.height = 0;
    other.mkl_descr = {};
  }

  // 移动赋值
  SparseMatrixComplex& operator=(SparseMatrixComplex&& other) noexcept {
    if (this != &other) {
      if (Use_MKL && mkl_sparse_matrix != nullptr) {
        mkl_sparse_destroy(mkl_sparse_matrix);
      }

      width = other.width;
      height = other.height;
      RowIndex = std::move(other.RowIndex);
      J = std::move(other.J);
      A = std::move(other.A);
      mkl_sparse_matrix = other.mkl_sparse_matrix;
      Use_MKL = other.Use_MKL;
      mkl_descr = other.mkl_descr;

      other.mkl_sparse_matrix = nullptr;
      other.Use_MKL = false;
      other.width = other.height = 0;
      other.mkl_descr = {};
    }
    return *this;
  }


  /// Returns reference to A[i][j].  Index i, j = 0 .. size-1
  const ComplexScalar &operator()(int i, int j) const;

  const ComplexScalar &Elem(int i, int j) const;

  void Mult(const ComplexVector &x, ComplexVector &y) const;
  void Mult(const ComplexMatrix &x, ComplexMatrix &y) const;

  void Add(const ComplexScalar& a,SparseMatrixComplex&B,SparseMatrixComplex& C);

  sparse_matrix_t& getMKLMatrix() {
    return mkl_sparse_matrix;
  }

  void Augment(const ComplexVector &vec, const ComplexScalar &diagVal);

  /// y += A * x (default)  or  y += a * A * x
  void AddMult(const ComplexVector &x, ComplexVector &y, const ComplexScalar &a) const;

  void AddMult(const ComplexMatrix &x, ComplexMatrix &y, const ComplexScalar &a) const;

  /// y += At * x (default)  or  y += a * At * x
  void AddMultTranspose(const ComplexVector &x, ComplexVector &y,
                        const ComplexScalar& a) const;

  /// y += AH * x (default)  or  y += a * AH * x
  void AddMultConjugateTranspose(const ComplexVector &x, ComplexVector &y,
                        const ComplexScalar& a) const;



  void covertToMKLCSR();

  int NumNonZeroElems() const {
    return RowIndex[width];
  }

  //IntArray GetI(){return RowIndex;};

  IntArray& GetI() { return RowIndex; }
  const IntArray& GetI() const { return RowIndex; }

  IntArray& GetJ() { return J; }
  const IntArray& GetJ() const { return J; }

  ComplexVector& GetA() { return A; }
  const ComplexVector& GetA() const { return A; }

  void EliminateElement(const IntArray& index, bool diagVal= true);

  // 假设 A, RowIndex, J 是可以直接访问的成员或指针
  // 如果不能直接访问，请使用对应的 get 接口，但要避免 Set() 的开销

  void ZeroRowColAndKeepDiag(int j, bool diagval)
  {
    if (j < 0 || j >= height)
      throw std::out_of_range("Index out of range");

    // 获取第 j 行的范围
    int row_start = RowIndex[j];
    int row_end   = RowIndex[j + 1];

    // ---------------------------------------------------------
    // 步骤 1: 遍历第 j 行
    // ---------------------------------------------------------
    for (int idx = row_start; idx < row_end; ++idx) {
      int col = J[idx]; // 这一行连接到的列（即相邻节点）

      // 情况 A: 对角线元素 (Row j, Col j)
      if (col == j) {
        if (!diagval) {
          A.Set(idx, ComplexScalar(0.0, 0.0));
        }
        // 如果 diagval == true，通常建议设为 1.0 而不是保留原值
        // A.Set(idx, ComplexScalar(1.0, 0.0));
        continue;
      }

      // 情况 B: 第 j 行的非对角元素 -> 清零
      A.Set(idx, ComplexScalar(0.0, 0.0));

      // ---------------------------------------------------------
      // 步骤 2: 利用对称性清零对应的列 (Row col, Col j)
      // ---------------------------------------------------------
      // 只有当矩阵结构对称时有效 (FEA 99% 的情况都满足)
      // 我们知道第 col 行一定包含一个指向 j 的元素

      int neighbor_row_start = RowIndex[col];
      int neighbor_row_end   = RowIndex[col + 1];

      // 在第 col 行中寻找列索引为 j 的元素
      // 由于单行元素很少（通常 < 100），这个线性搜索非常快
      for (int k = neighbor_row_start; k < neighbor_row_end; ++k) {
        if (J[k] == j) {
          A.Set(k, ComplexScalar(0.0, 0.0));
          break; // 找到了就停止，不需要继续扫该行
        }
      }
    }
  }


  void PrintMatlab(std::ostream &out) const {
    std::ios::fmtflags old_fmt = out.setf(std::ios::scientific);
    int old_prec = out.precision(14);

    // 遍历行
    for (int i = 0; i < height; ++i) {
      for (int idx = RowIndex[i]; idx < RowIndex[i + 1]; ++idx) {
        int j = J[idx];
        const ComplexScalar &val = A[idx];

        // 只输出非零元素
        if (val.real() != 0.0 || val.imag() != 0.0) {
          out << (i + 1) << " " << (j + 1) << " ";
          // MATLAB 风格复数 a + bi
          out << val.real()
              << (val.imag() >= 0 ? "+" : "-")
              << std::abs(val.imag()) << "i"
              << std::endl;
        }
      }
    }

    out.precision(old_prec);
    out.setf(old_fmt);
  }


  // 输出矩阵信息
  void print() const {
    std::cout << "SparseMatrixComplex: nnz=" << width << "\n";
    std::cout << "RowIndex: ";
    for (int i = 0; i < width+1; i++) std::cout << RowIndex[i] << " ";
    std::cout << "\nJ: ";
    for (int i = 0; i < RowIndex[width]; i++) std::cout << J[i] << " ";
    std::cout << "\nA: ";
    for (int i = 0; i < RowIndex[width]; i++) std::cout << "(" << A[i].real() << "," << A[i].imag() << ") ";
    std::cout << std::endl;
  }

  void SortColumnIndices();

  ~SparseMatrixComplex() {
    if (Use_MKL == true) {
      mkl_sparse_destroy(mkl_sparse_matrix);
      mkl_sparse_matrix = nullptr;
      Use_MKL=false;
    }
  };

};
#endif
