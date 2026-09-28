#pragma once
#include "operator.hpp"

class MatrixInverse;

/// Abstract data type matrix
class Matrix : public Operator {
  friend class MatrixInverse;

public:
  /// Creates matrix of width s.
  Matrix(int s) { size = s; }

  /// Returns reference to a_{ij}.  Index i, j = 0 .. size-1
  virtual double &Elem(int i, int j) = 0;

  /// Returns constant reference to a_{ij}.  Index i, j = 0 .. size-1
  virtual const double &Elem(int i, int j) const = 0;

  /// Returns a pointer to (approximation) of the matrix inverse.
  virtual MatrixInverse *Inverse() const = 0;

  /// Finalizes the matrix initialization.
  virtual void Finalize(int) {}

  /// Prints matrix to stream out.
  virtual void Print(std::ostream &out = std::cout, int width = 4) const;

  /// Destroys matrix.
  virtual ~Matrix() {}
};

/// Abstract data type for matrix inverse
class MatrixInverse : public Operator {
protected:
  const Matrix *a;

public:
  /// Creates approximation of the inverse of square matrix
  MatrixInverse(const Matrix &mat) {
    size = mat.size;
    a = &mat;
  }

  /// Destroys inverse matrix.
  virtual ~MatrixInverse() {}
};