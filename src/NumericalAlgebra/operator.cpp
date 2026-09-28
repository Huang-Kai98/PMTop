#include "operator.hpp"
#include "floatarray.h"
#include <iomanip>

void Operator::PrintMatlab(std::ostream &out, int n, int m) {
  if (n == 0)
    n = size;
  if (m == 0)
    m = size;

  FloatArray x(n), y(m);
  x = 0.0;

  int i, j;
  out << setiosflags(std::ios::scientific | std::ios::showpos);
  for (i = 0; i < n; i++) {
    if (i != 0)
      x(i - 1) = 0.0;
    x(i) = 1.0;
    Mult(x, y);
    for (j = 0; j < m; j++)
      if (y(j))
        out << j + 1 << " " << i + 1 << " " << y(j) << std::endl;
  }
}
