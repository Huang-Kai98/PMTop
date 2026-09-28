// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#include "error.hpp"
void pmtop_error(const char *msg) {
  if (msg)
    std::cerr << msg << std::endl;
  *((int *)NULL) = 0; // force crash by causing segmentation fault
}
void pmtop_assert(const char *msg) {
  if (msg) {
    std::cout << msg << std::endl;
    std::abort();
  }
  *((int *)NULL) = 0; // force crash by causing segmentation fault
}