// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#pragma once
#include "error.hpp"
#include "floatarray.h"
#include <iostream>
/// Abstract operator
class Operator {
protected:
  int size;
  int height; ///< Dimension of the output / number of rows in the matrix.
  int width;  ///< Dimension of the input / number of columns in the matrix.


public:
  /// Construct Operator with given size s (default 0)
  Operator(int s = 0) { size = s; }

  /// Returns the size of the input
  inline int Size() const { return size; }

  /// Get the height (size of output) of the Operator. Synonym with NumRows().
  inline int Height() const { return height; }
  /** @brief Get the number of rows (size of output) of the Operator. Synonym
      with Height(). */
  inline int NumRows() const { return height; }

  /// Get the width (size of input) of the Operator. Synonym with NumCols().
  inline int Width() const { return width; }
  /** @brief Get the number of columns (size of input) of the Operator. Synonym
      with Width(). */
  inline int NumCols() const { return width; }



  /// Operator application
  virtual void Mult(const FloatArray &x, FloatArray &y) const = 0;

  /// Action of the transpose operator
  virtual void MultTranspose(const FloatArray &x, FloatArray &y) const {
    pmtop_error("Operator::MultTranspose not implemented.");
  }

  /// Prints operator with input size n and output size m in matlab format.
  void PrintMatlab(std::ostream &out, int n = 0, int m = 0);

  virtual ~Operator() {}
};

/// Operator I: x -> x
class IdentityOperator : public Operator {
public:
  /// Creates I_{nxn}
  IdentityOperator(int n) { size = n; }

  /// Operator application
  virtual void Mult(const FloatArray &x, FloatArray &y) const { y = x; }

  ~IdentityOperator() {}
};

/// The transpose of a given operator (square matrix)
class TransposeOperator : public Operator {
private:
  Operator *A;

public:
  /// Saves the operator
  TransposeOperator(Operator *a) : A(a) { size = A->Size(); }

  /// Operator application
  virtual void Mult(const FloatArray &x, FloatArray &y) const {
    A->MultTranspose(x, y);
  }

  virtual void MultTranspose(const FloatArray &x, FloatArray &y) const {
    A->Mult(x, y);
  }

  ~TransposeOperator() {}
};