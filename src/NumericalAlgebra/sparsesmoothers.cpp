#include "sparsesmoothers.hpp"
#include "floatarray.h"
#include "intarray.h"

/// Create GSSmoother.
GSSmoother::GSSmoother(const SparseMatrix &a) : MatrixInverse(a) {}

/// Matrix vector multiplication with GS Smoother.
void GSSmoother::Mult(const FloatArray &x, FloatArray &y) const {
  y = 0.;

  ((SparseMatrix *)a)->Gauss_Seidel_forw(x, y);
  ((SparseMatrix *)a)->Gauss_Seidel_back(x, y);
}

/// Destroys the GS Smoother.
GSSmoother::~GSSmoother() {}

/// Create the diagonal smoother.
DSmoother::DSmoother(const SparseMatrix &a, double s) : MatrixInverse(a) {
  scale = s;
}

/// Matrix vector multiplication with Diagonal smoother.
void DSmoother::Mult(const FloatArray &x, FloatArray &y) const {
  for (int i = 0; i < x.giveSize(); i++)
    y(i) = scale * x(i) / a->Elem(i, i);
}

/// Destroys the Diagonal smoother.
DSmoother::~DSmoother() {}
