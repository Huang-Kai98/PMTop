// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

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
