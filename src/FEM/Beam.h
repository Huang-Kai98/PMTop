#pragma once
#include "Element.h"
class Beam : public Element {};
class Beam188 : public Beam {
public:
  Beam188(){};
  void init() override {
    dof = 6;
    nr = 12;
    jtyp = "Beam188";
    nodes.reserve(2);
    elementStif = new double[12 * 12];
    elementMass = new double[12 * 12];
    elementDamp = new double[12 * 12];
  }

  void setStiff(const int i, const int j, const double value) override {
    elementStif[j * 12 + i] = value;
  }

  void setMass(const int i, const int j, const double value) override {
    elementMass[j * 12 + i] = value;
  }

  void setDamp(const int i, const int j, const double value) override {
    elementDamp[j * 12 + i] = value;
  }

  ~Beam188() {
    delete[] elementStif;
    delete[] elementMass;
    delete[] elementDamp;
  }
};
class COMBIN14 : public Beam {
public:
  COMBIN14(){};
  void init() override {
    dof = 1;
    nr = 2;
    jtyp = "COMBIN14";
    nodes.reserve(2);
    elementStif = new double[2 * 2];
    elementMass = new double[2 * 2];
    elementDamp = new double[2 * 2];
  }

  void setStiff(const int i, const int j, const double value) override {
    elementStif[j * nr + i] = value;
  }

  void setMass(const int i, const int j, const double value) override {
    elementMass[j * nr + i] = value;
  }

  void setDamp(const int i, const int j, const double value) override {
    elementDamp[j * nr + i] = value;
  }

  ~COMBIN14() {
    delete[] elementStif;
    delete[] elementMass;
    delete[] elementDamp;
  }
};