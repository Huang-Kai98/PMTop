// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

//
// Created by huangkai on 25-1-23.
//

#ifndef ANSYSSOLVER_H
#define ANSYSSOLVER_H
#include <floatarray.h>
#include <pardiso.hpp>
#include <Arpack.hpp>
#include <functional>
#include <random>
#include "../FEM/AnsysLoad.h"
#include "../FEM/BilinearForm.h"


class AnsysSolver {
public:
    AnsysSolver(FloatArray &rhoInput, AnsysLoad &ansysInput, std::string test);

    AnsysSolver(FloatArray &rhoInput, AnsysLoad &ansysInput);

    AnsysSolver(FloatArray &rhoInput, AnsysLoad &ansysInput, std::function<double(double x)> StifRho_,
                std::function<double(double x)> MassRho_, std::function<double(double x)> StifRhoD_ = nullptr,
                std::function<double(double x)> MassRhoD_ = nullptr,
                std::function<double(double x)> StifRhoD2_ = nullptr,
                std::function<double(double x)> MassRho2_ = nullptr);

    AnsysSolver(FloatArray &rhoInput, AnsysLoad &ansysInput, std::function<double(double x)> StifRho_,
                std::function<double(double x)> MassRho_, std::function<double(double x)> DampRho_,
                std::function<double(double x)> StifRhoD_ = nullptr,
                std::function<double(double x)> MassRhoD_ = nullptr,
                std::function<double(double x)> DampRhoD_ = nullptr,
                std::function<double(double x)> StifRhoD2_ = nullptr,
                std::function<double(double x)> MassRho2_ = nullptr,
                std::function<double(double x)> DampRhoD2_ = nullptr);

    void EigenProblem(int nev);

    void QudraticEigenProblem(int nev, std::string SolutionOpt);

    FloatArray getEigenVectors();

    FloatArray getEigenValues();

    FloatArray getEigenVector(int i) const;

    FloatMatrix getEigenvector(IntArray &idx) const;

    FloatMatrix getEigenvector(std::vector<int> &idx) const;

    ComplexMatrix getComplexEigenvectors() {
        return arpack_solver->getComplexEigenvectors();
    }

    ComplexVector getComplexEigenvalues() {
        return arpack_solver->getComplexEigenvalues();
    }

    ComplexVector getComplexEigenvalue(IntArray &idx) {
        return arpack_solver->getComplexEigenvalue(idx);
    }

    ComplexMatrix getComplexEigenvector(IntArray &idx) const {
        return arpack_solver->getComplexEigenvector(idx);
    }


    IntArray getInfo(std::string name);

    IntArray getI(std::string name);

    IntArray getJ(std::string name);

    FloatArray getValue(std::string name);

    SparseMatrix *createSparseMatrix(std::string name);

    void TestEigenProblem();

    void NelsonRepeat(ComplexScalar &EigenValue, ComplexMatrix &EigenVector, SparseMatrixComplex &STIF,
                      SparseMatrixComplex &MASS, SparseMatrixComplex &DAMP, int i, ComplexVector &DeriviateEigenValue,
                      ComplexMatrix &DeriviateEigenVector);

    void NelsonSingle(ComplexScalar &EigenValue, ComplexVector &EigenVector, SparseMatrixComplex &STIF,
                      SparseMatrixComplex &MASS, SparseMatrixComplex &DAMP, int i, ComplexScalar &DeriviateEigenValue,
                      ComplexVector &DeriviateEigenVector);

    void NelsonSinglePrepareFactorized(
    ComplexScalar &EigenValue,
    ComplexVector &EigenVector,
    SparseMatrixComplex &STIF,
    SparseMatrixComplex &MASS,
    SparseMatrixComplex &DAMP,
    SparseMatrixComplex &A1,
    IntArray &g,
    PardisoSolver &PA,
    ComplexVector &Cphi,
    ComplexVector &Mphi,
    ComplexVector &RhsCoeff);


    void NelsonSingleBackSubFactorized(
    ComplexScalar &EigenValue,
    ComplexVector &EigenVector,
    SparseMatrixComplex &STIF,
    SparseMatrixComplex &MASS,
    SparseMatrixComplex &DAMP,
    const int index,
    PardisoSolver &PA,
    IntArray &g,
    ComplexVector &Cphi,
    ComplexVector &Mphi,
    ComplexVector &RhsCoeff,
    ComplexScalar &DeriviateEigenValue,
    ComplexVector &DeriviateEigenVector);


    void NelsonAddPrepareFactorized(
        ComplexScalar &EigenValue,
        ComplexVector &EigenVector,
        SparseMatrixComplex &STIF,
        SparseMatrixComplex &MASS,
        SparseMatrixComplex &DAMP,
        SparseMatrixComplex &Aadd,
        PardisoSolver &PA,
        ComplexVector &Mphi,
        ComplexVector &RhsCoeff
    );

    void NelsonAddBackSubFactorized(
        ComplexScalar &EigenValue,
        ComplexVector &EigenVector,
        SparseMatrixComplex &STIF,
        SparseMatrixComplex &MASS,
        SparseMatrixComplex &DAMP,
        const int index,
        PardisoSolver &PA,
        ComplexVector &Mphi,
        ComplexVector &RhsCoeff,
        ComplexScalar &DeriviateEigenValue,
        ComplexVector &DeriviateEigenVector
    );

    void NelsonAdd(ComplexScalar &EigenValue, ComplexVector &EigenVector, SparseMatrixComplex &STIF,
                   SparseMatrixComplex &MASS, SparseMatrixComplex &DAMP, int i, ComplexScalar &DeriviateEigenValue,
                   ComplexVector &DeriviateEigenVector);

    void AdjointNelson(ComplexScalar &EigenValue, ComplexVector &EigenVector, SparseMatrixComplex &STIF,
                       SparseMatrixComplex &MASS, SparseMatrixComplex &DAMP, ComplexScalar Dlembda,
                       ComplexVector &Dsymvector,
                       ComplexScalar &alpha, ComplexVector &V);

    void AdjointNelosnAdd(ComplexScalar &EigenValue, ComplexVector &EigenVector, SparseMatrixComplex &STIF,
                          SparseMatrixComplex &MASS, SparseMatrixComplex &DAMP, ComplexScalar Dlembda,
                          ComplexVector &Dsymvector,
                          ComplexScalar &alpha, ComplexVector &V);

    void AdjointsymmDsymvector(ComplexScalar &EigenValue, ComplexVector &EigenVector, SparseMatrixComplex &STIF,
                               SparseMatrixComplex &MASS, SparseMatrixComplex &DAMP, const PardisoSolver &PA,
                               ComplexVector &Dsymvector,
                               ComplexVector &Dphi, ComplexScalar &DphiL, ComplexScalar &EigenValueSquare,
                               ComplexScalar &Db3Dphib0,
                               const double tol = 1.0e-10, const int maxit = 5);

    double MIR(ComplexVector &EigenVector);

    double MPC(ComplexVector &EigenVector);

    double MAC(ComplexVector &EigenVector_Measured, ComplexVector &EigenVector);

    ComplexVector DMIR(ComplexVector &EigenVector);

    ComplexVector DMPC(ComplexVector &EigenVector);

    ComplexVector DMAC(ComplexVector &EigenVector_Measured, ComplexVector &EigenVector);

    ComplexScalar AdjointElement(ComplexScalar &EigenValue, ComplexVector &EigenVector, const ComplexScalar &alpha,
                                 ComplexVector &V, const int index);

    ComplexScalar AdjointsymmElement(ComplexScalar &EigenValue, ComplexVector &EigenVector, ComplexVector &Dphi,
                                     ComplexScalar &DphiL, const ComplexScalar &EigenValueSquare,
                                     const ComplexScalar &Db3Dphib0, const ComplexScalar &Dlemda,
                                     const int index);

    void Normalize(const ComplexScalar &EigenValue, ComplexVector &EigenVector, const SparseMatrixComplex &DAMP,
                   const SparseMatrixComplex &MASS);

    ComplexMatrix AddNoiseToEigenvectors(const ComplexMatrix& EigenVector,double NoiseLevel = 0.2);

    void Nelson2(ComplexScalar &EigenValue, ComplexMatrix &EigenVector, const SparseMatrix &STIF,
                 const SparseMatrix &MASS, const
                 SparseMatrix &DAMP, const SparseMatrix &STIFD, const SparseMatrix &MASSD, const SparseMatrix &DAMPD,
                 const SparseMatrix
                 &STIFD2, const SparseMatrix &MASSD2, const SparseMatrix &DAMPD2);


    void SetEigenvector(int index, double value);

    void Reorthogonalization();

    void changeEigenvector(IntArray index, FloatMatrix values);

    void Mult(const FloatArray &b, FloatArray &x) { pardiso_solver->Mult(b, x); }

    void Mult(double *b, double *x) { pardiso_solver->Mult(b, x); }

    int getSize() { return ansys.getNumberOfDof(); };
    int getElementSize() { return ansys.getNumberOfElements(); };
    IntArray getElementLocation(int i) { return ansys.getElementLocation(i); };
    FloatMatrix getElementStif(int i, double x) const { return ansys.getElementStif(i, x); };
    FloatMatrix getElementMass(int i, double x) const { return ansys.getElementMass(i, x); };
    FloatMatrix getElementDamp(int i, double x) const { return ansys.getElementDamp(i, x); };


    double getRho(int i) { return rho[i]; };
    SparseMatrix *getStif() { return Stif; };
    SparseMatrix *getMass() { return Mass; };
    SparseMatrix *getDamp() { return Damp; };
    double getStifRho(double x) { return StifRho(x); };
    double getMassRho(double x) { return MassRho(x); };
    double getDampRho(double x) { return DampRho(x); };
    double getStifRhoD(double x) { return StifRhoD(x); };
    double getMassRhoD(double x) { return MassRhoD(x); };
    double getDampRhoD(double x) { return DampRhoD(x); };
    double getStifRhoD2(double x) { return StifRhoD2(x); };
    double getMassRhoD2(double x) { return MassRhoD2(x); };
    double getDampRhoD2(double x) { return DampRhoD2(x); };

    void Adjacency_list_element_node() { ansys.Adjacency_list_element_node(); };
    void Adjacency_list_element() { ansys.Adjacency_list_element(); };
    void Adjacency_list_node() { ansys.Adjacency_list_node(); };

    std::map<int, std::set<int> > Adjacent_cells() { return ansys.adjacent_cells; };
    std::map<int, std::set<int> > Adjacent_nodes() { return ansys.adjacent_cells; };

    KDNode *GetKDNode(std::vector<Node> &points, int dim = 3, int depth = 0) {
        return ansys.buildKDTree(points, dim, depth);
    };
    void printKDTree(const KDNode *node) { ansys.printKDTree(node); };
    std::vector<Node> getAnsysNode() { return ansys.getAnasysNodes(); };
    std::vector<Node> getCentroid() { return ansys.getCentroid(); };

    std::vector<Node> rangeSearch(KDNode *node, const Node &target, double radius, int depth = 0, int dim = 3) {
        return ansys.rangeSearch(node, target, radius, depth, dim);
    };

private:
    FloatArray &rho;
    std::shared_ptr<BilinearForm> a = nullptr;
    std::shared_ptr<PardisoSolver> pardiso_solver;
    std::shared_ptr<ArpackSolver> arpack_solver = nullptr;
    AnsysLoad &ansys;
    SparseMatrix *Stif = nullptr;
    SparseMatrix *Mass = nullptr;
    SparseMatrix *Damp = nullptr;

    int NumberOfDof;
    std::function<double(double x)> StifRho;
    std::function<double(double x)> MassRho;
    std::function<double(double x)> DampRho;
    std::function<double(double x)> StifRhoD;
    std::function<double(double x)> MassRhoD;
    std::function<double(double x)> DampRhoD;
    std::function<double(double x)> StifRhoD2;
    std::function<double(double x)> MassRhoD2;
    std::function<double(double x)> DampRhoD2;
    friend class Top;
    friend class Sensitivity;
    friend class TopProblem;
};


class Sensitivity {
public:
    Sensitivity(std::shared_ptr<AnsysSolver> ansys_solver_) {
        ansys_solver = std::move(ansys_solver_);
    };

    FloatArray AdjointMethod(double Eval, FloatArray &Evec, FloatArray &obj, FloatArray &EigenValueDerivative,
                             FloatArray &EigenVectorDerivative
                             , int MaxIter = 20, double rol = 1e-10) const;

    FloatMatrix DREV(int index, double Eval, FloatArray &EvalD, FloatMatrix &EvecD);

private:
    std::shared_ptr<AnsysSolver> ansys_solver;
};

class Top {
public:
    Top(std::shared_ptr<AnsysSolver> ansys_solver_) {
        ansys_solver = std::move(ansys_solver_);
        x0 = 1.0;
        xmin = 0.01;
        xmax = 1.0;
    }

    void SetTopElement(IntArray set);

    void CreateDensityFilteringMatrix(double radius = std::sqrt(3));

private:
    std::shared_ptr<AnsysSolver> ansys_solver;
    std::unique_ptr<SparseMatrix> DensityFilteringMatrix;
    FloatArray HS;
    FloatArray x0;
    FloatArray xmin;
    FloatArray xmax;
    friend class TopProblem;
};

struct SubEigenProblem {
    FloatArray SubEigenValue;
    FloatMatrix SubEigenVector;

    SubEigenProblem(const FloatArray &SubEigenValue_, const FloatMatrix &SubEigenVector_) {
        SubEigenValue = SubEigenValue_;
        SubEigenVector = SubEigenVector_;
    }
};

#endif //ANSYSSOLVER_H
