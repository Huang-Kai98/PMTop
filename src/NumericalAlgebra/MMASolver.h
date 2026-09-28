//
// Created by huangkai on 25-1-16.
//

#ifndef MMASOLVER_H
#define MMASOLVER_H



#include <vector>

class MMASolver {

public:
    MMASolver(int n, int m, double a = 0.0, double c = 1000.0, double d = 0.0);

    void SetAsymptotes(double init, double decrease, double increase);

    void ConstraintModification(bool conMod) {}

    void Update(double *xval, const double *dfdx, const double *gx, const double *dgdx, const double *xmin,
                const double *xmax);

    void Reset() { iter = 0; };

private:
    int n, m, iter;

    const double xmamieps;
    const double epsimin;

    const double raa0;
    const double move, albefa;
    double asyminit, asymdec, asyminc;

    std::vector<double> a, c, d;
    std::vector<double> y;
    double z;

    std::vector<double> lam, mu, s;
    std::vector<double> low, upp, alpha, beta, p0, q0, pij, qij, b, grad, hess;

    std::vector<double> xold1, xold2;

    void GenSub(const double *xval, const double *dfdx, const double *gx, const double *dgdx, const double *xmin,
                const double *xmax);

    void SolveDSA(double *x);
    void SolveDIP(double *x);

    void XYZofLAMBDA(double *x);

    void DualGrad(double *x);
    void DualHess(double *x);
    void DualLineSearch();
    double DualResidual(double *x, double epsi);

    static void Factorize(double *K, int n);
    static void Solve(double *K, double *x, int n);
};



#endif //MMASOLVER_H
