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