// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

//
// Created by huangkai on 25-1-16.
//

#include "Optimizer.h"

void Optimizer::SolveGCMMA(int maxIterations, double tolerance) {
    double f, fnew;
    double ch = 1.0;

    problem.objFunc(problem.x, &f, problem.g);
    PrintVector(problem.x, "Initial x");
    std::cout << "Initial objective: " << f << std::endl;
    PrintVector(problem.g, "Initial constraints");

    GCMMASolver gcmma(problem.n,problem.m);
    for (int iter = 0; ch > tolerance && iter < maxIterations; ++iter) {
        problem.sensFunc(problem.x ,&f, problem.g, problem.df,problem.dg);
        gcmma.OuterUpdate(problem.xnew.givePointer(), problem.x.givePointer(), f, problem.df.givePointer(),
            problem.g.givePointer(), problem.dg.givePointer(), problem.xmin.givePointer(), problem.xmax.givePointer());
        problem.objFunc(problem.xnew, &fnew, problem.gnew);
        bool conserv = gcmma.ConCheck(fnew, problem.gnew.givePointer());
        for (int innerIter = 0; !conserv && innerIter < 15; ++innerIter) {
            gcmma.InnerUpdate(problem.xnew.givePointer(), fnew, problem.gnew.givePointer(), problem.x.givePointer(),
                f, problem.df.givePointer(), problem.g.givePointer(), problem.dg.givePointer(),
                problem.xmin.givePointer(), problem.xmax.givePointer());
            problem.objFunc(problem.xnew, &fnew, problem.gnew);
            conserv = gcmma.ConCheck(fnew, problem.gnew.givePointer());
        }
        problem.x = problem.xnew;
        ch = ComputeChange(problem.x, problem.xold);
        printf("Iteration: %d, Objective: %f, Change: %f\n", iter, f, ch);
        PrintVector(problem.x, "x");
        problem.objFunc(problem.x, &f, problem.g);
        PrintVector(problem.g, "Constraints");
        std::cout << std::endl;
    }
}

void Optimizer::SolveMMA(int maxIterations, double tolerance) {
    double f, fnew;
    double ch = 1.0;

    problem.objFunc(problem.x, &f, problem.g);
    PrintVector(problem.x, "Initial x");
    std::cout << "Initial objective: " << f << std::endl;
    PrintVector(problem.g, "Initial constraints");

    MMASolver mma(problem.n,problem.m);
    for (int iter = 0; ch > tolerance && iter < maxIterations; ++iter) {
        problem.sensFunc(problem.x ,&f, problem.g, problem.df,problem.dg);
        mma.Update(problem.x.givePointer(), problem.df.givePointer(), problem.g.givePointer(),
            problem.dg.givePointer(),problem.xmin.givePointer(), problem.xmax.givePointer());
        ch = ComputeChange(problem.x, problem.xold);
        printf("Iteration: %d, Objective: %f, Change: %f\n", iter, f, ch);
        PrintVector(problem.x, "x");
        problem.objFunc(problem.x, &f, problem.g);
        PrintVector(problem.g, "Constraints");
        std::cout << std::endl;
    }
}

void Optimizer::Solver(int maxIterations, double tolerance, std::string method) {
    std::cout<<"Method: "<< method <<std::endl;
    std::cout<<"MaxIterations: "<<maxIterations<<std::endl;
    std::cout<<"Tolerance: "<<tolerance<<std::endl;
    if (method == "GCMMA") {SolveGCMMA(maxIterations, tolerance);}
    else if (method == "MMA") {SolveMMA(maxIterations, tolerance);}
    else {std::cerr<<"Unkonw Method: "<<method<<std::endl;}
}


double Optimizer::ComputeChange(const FloatArray &x, FloatArray &xold) {
    double ch = 0.0;
    for (int i = 0; i < x.giveSize(); ++i) {
        ch = std::max(ch, std::abs(x[i] - xold[i]));
        xold[i] = x[i];
    }
    return ch;
}
