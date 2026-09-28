// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

//
// Created by huangkai on 25-1-16.
//

#ifndef OPTIMIZER_H
#define OPTIMIZER_H
#include "GCMMASolver.h"
#include "MMASolver.h"
#include <iostream>
#include <utility>
#include <vector>
#include <cmath>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>

#include "floatarray.h"
#include "floatmatrix.h"

class Problem {
public:
    Problem(int n, int m,
              std::function<void(const FloatArray&, double*, FloatArray&)> objFunc,
              std::function<void(const FloatArray&, double*, FloatArray&, FloatArray&, FloatArray&)> sensFunc,
              FloatArray  x0,
              FloatArray  xmin,
              FloatArray  xmax)
        : n(n), m(m),
          objFunc(std::move(objFunc)),
          sensFunc(std::move(sensFunc)),
          x(x0),
          xold(x),
          xnew(n),
          df(n),
          g(m),
          gnew(m),
          dg(n * m),
          xmin(xmin),
          xmax(xmax)
    {}

    int n, m;
    FloatArray x, xold, xnew;
    FloatArray df, g, gnew, dg;
    FloatArray xmin, xmax;

    /* Parameter description:
     * objFunc: x, f0x, fx
     * sensFunc: x, f0x, fx, df0dx, dfdx
     */
    std::function<void(const FloatArray&, double*, FloatArray&)> objFunc;
    std::function<void(const FloatArray&, double*, FloatArray&, FloatArray&, FloatArray&)> sensFunc;
};

class Optimizer {
public:
    explicit Optimizer(Problem& problem0): problem(problem0) {};
    void SolveGCMMA(int maxIterations, double tolerance);
    void SolveMMA(int maxIterations, double tolerance);
    void Solver(int maxIterations = 8, double tolerance = 0.0002,std::string method="GCMMA");
private:
    Problem& problem;
    static void PrintVector(const FloatArray &vec, const std::string &name) {
        std::cout << name << ":";
        for (const double v : vec) {
            std::cout << " " << v;
        }
        std::cout << std::endl;
    }
    double ComputeChange(const FloatArray &x, FloatArray &xold);
};


#endif //OPTIMIZER_H
