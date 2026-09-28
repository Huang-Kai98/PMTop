// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

//
// Created by huangkai on 25-2-11.
//
#pragma once
#include "Element.h"
class Solid : public Element {};

class Solid185: public Solid {
  public:
    Solid185() {};
    void init() override {
        dof = 3;
        nr = 24;
        jtyp = "Solid185";
        nodes.reserve(8);
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

    ~Solid185() {
        delete[] elementStif;
        delete[] elementMass;
        delete[] elementDamp;
    }

};