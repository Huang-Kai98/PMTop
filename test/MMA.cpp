
#include "NumericalAlgebra.hpp"
#include "Optimizer.h"
#include <fstream>
#include <iostream>
#include <mkl.h>
#include <mkl_cblas.h>
#include <mkl_service.h>
#include <mkl_spblas.h>
#include <mkl_types.h>
int main() {
    auto objFunc = [](const FloatArray& x, double *f0x, FloatArray& fx) {
        f0x[0] = 0;
        for (int i = 0; i < 3; ++i) {
            f0x[0] += x[i] * x[i];
        }
        fx[0] = std::pow(x[0] - 5, 2) + std::pow(x[1] - 2, 2) + std::pow(x[2] - 1, 2) - 9;
        fx[1] = std::pow(x[0] - 3, 2) + std::pow(x[1] - 4, 2) + std::pow(x[2] - 3, 2) - 9;
    };

    auto sensFunc = [objFunc](const FloatArray& x, double *f0x, FloatArray& fx, FloatArray& df0dx, FloatArray& dfdx) {
        objFunc(x, f0x, fx);  // 显式捕获 objFunc
        for (int i = 0; i < 3; ++i) {
            df0dx[i] = 2 * x[i];
        }
        int k = 0;
        dfdx[k++] = 2 * (x[0] - 5); dfdx[k++] = 2 * (x[0] - 3);
        dfdx[k++] = 2 * (x[1] - 2); dfdx[k++] = 2 * (x[1] - 4);
        dfdx[k++] = 2 * (x[2] - 1); dfdx[k++] = 2 * (x[2] - 3);
    };


    int n = 3;
    int m = 2;
    FloatArray x0 = {4, 3, 2};
    FloatArray xmin(n);
    xmin=0.0;
    FloatArray xmax(n);
    xmax=5.0;
    Problem problem(n,m,objFunc,sensFunc,x0,xmin,xmax);


    Optimizer optimizer(problem);
    optimizer.Solver(10,0.0002,"MMA");



    return 0;
}