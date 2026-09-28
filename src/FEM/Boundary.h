// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#pragma once
#include <ostream>
#include <vector>

class Boundary {
public:
  std::vector<int> MatrixLocalTion;
  std::vector<std::string> BoundaryLabel;
  std::vector<int> RealValue;
  std::vector<int> ImagValue;

  // 将数据写入输出流
  void saveToStream(std::ostream &out) const {
    out << "Matrix Location: ";
    for (const int &loc : MatrixLocalTion) {
      out << loc << " ";
    }
    out << std::endl;

    out << "Boundary Labels: ";
    for (const std::string &label : BoundaryLabel) {
      out << label << " ";
    }
    out << std::endl;

    out << "Real Values: ";
    for (const int &value : RealValue) {
      out << value << " ";
    }
    out << std::endl;

    out << "Imaginary Values: ";
    for (const int &value : ImagValue) {
      out << value << " ";
    }
    out << std::endl;
  }
};

class BoundaryData {
public:
  int nodeId;
  std::string dofLabel;
  double realPart;
  double imagPart;
};
