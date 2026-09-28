// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#pragma once
#include "Element.h"
class Shell : public Element {};

class Shell181 : public Shell {
public:
  Shell181(){};

  void init() override {
    dof = 6;
    nr = 24;
    jtyp = "Shell181";
    nodes.reserve(4);
    elementStif = new double[24 * 24];
    elementMass = new double[24 * 24];
    elementDamp = new double[24 * 24];
  }

  void setStiff(const int i, const int j, const double value) override {
    elementStif[j * 24 + i] = value;
  }

  void setMass(const int i, const int j, const double value) override {
    elementMass[j * 24 + i] = value;
  }

  void setDamp(const int i, const int j, const double value) override {
    elementDamp[j * 24 + i] = value;
  }

  ~Shell181() {
    delete[] elementStif;
    delete[] elementMass;
    delete[] elementDamp;
  }
};

class Shell63 : public Shell {
public:
  Shell63(){};
  void init() override {
    dof = 6;
    nr = 24;
    jtyp = "Shell63";
    nodes.reserve(4);
    elementStif = new double[24 * 24];
    elementMass = new double[24 * 24];
    elementDamp = new double[24 * 24];
  }

  void setStiff(const int i, const int j, const double value) override {
    elementStif[j * 24 + i] = value;
  }

  void setMass(const int i, const int j, const double value) override {
    elementMass[j * 24 + i] = value;
  }

  void setDamp(const int i, const int j, const double value) override {
    elementDamp[j * 24 + i] = value;
  }

  ~Shell63() {
    delete[] elementStif;
    delete[] elementMass;
    delete[] elementDamp;
  }
};