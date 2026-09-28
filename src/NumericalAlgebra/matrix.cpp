#include "matrix.hpp"
#include <iomanip>
#include <iostream>

void Matrix::Print(std::ostream &out, int width) const {
  // output flags = scientific + show sign
  out << setiosflags(std::ios::scientific | std::ios::showpos);
  for (int i = 0; i < size; i++) {
    out << "[row " << i << "]\n";
    for (int j = 0; j < size; j++) {
      out << Elem(i, j) << " ";
      if (!((j + 1) % width))
        out << std::endl;
    }
    out << std::endl;
  }
  out << std::endl;
}
