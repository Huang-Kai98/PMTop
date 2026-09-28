// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

//
// Created by huangkai on 25-1-23.
//

#include "TOP.h"

#include <chrono>

#include "../FEM/BilinearForm.cpp"
#include "pardiso.hpp"
#include "Arpack.hpp"

#include <iostream>
#include <tic_toc.hpp>

AnsysSolver::AnsysSolver(FloatArray &rhoInput, AnsysLoad &ansysInput, std::string test) : rho(rhoInput),
    ansys(ansysInput) {
    Stif = createSparseMatrix("K");
    Mass = createSparseMatrix("M");
    Damp = createSparseMatrix("C");

    Stif->covertToMKLCSR();
    Mass->covertToMKLCSR();
    Damp->covertToMKLCSR();

    pardiso_solver = std::make_shared<PardisoSolver>();
    pardiso_solver->SetMatrixType(PardisoSolver::REAL_NONSYMMETRIC);
    pardiso_solver->SetOperator(*Stif);
}

AnsysSolver::AnsysSolver(FloatArray &rhoInput, AnsysLoad &ansysInput) : rho(rhoInput), ansys(ansysInput) {
    a = std::make_shared<BilinearForm>(&ansys);
    std::unordered_set<int> ConstrainedDofs;
    ansys.getConstrainedDofs(ConstrainedDofs);
    a->assemble(rho);
    a->EliminateEssentialBC(ConstrainedDofs, 1);
    a->Finalize();
    a->convertToMKL();
    Stif = &a->SpStif();
    Mass = &a->SpMass();
    pardiso_solver = std::make_shared<PardisoSolver>();
    pardiso_solver->SetMatrixType(PardisoSolver::REAL_NONSYMMETRIC);
    pardiso_solver->SetOperator(*Stif);
}

AnsysSolver::AnsysSolver(FloatArray &rhoInput, AnsysLoad &ansysInput, std::function<double(double x)> StifRho_,
                         std::function<double(double x)> MassRho_, std::function<double(double x)> StifRhoD_,
                         std::function<double(double x)> MassRhoD_, std::function<double(double x)> StifRhoD2_,
                         std::function<double(double x)> MassRho2_) : rho(rhoInput),
                                                                      ansys(ansysInput),
                                                                      StifRho(StifRho_),
                                                                      MassRho(MassRho_),
                                                                      StifRhoD(StifRhoD_),
                                                                      MassRhoD(MassRhoD_),
                                                                      StifRhoD2(StifRhoD2_),
                                                                      MassRhoD2(MassRho2_) {
    a = std::make_shared<BilinearForm>(&ansys);
    std::unordered_set<int> ConstrainedDofs;
    ansys.getConstrainedDofs(ConstrainedDofs);
    a->assemble(rho, StifRho, MassRho);
    a->EliminateEssentialBC(ConstrainedDofs, 1);
    a->Finalize();
    a->convertToMKL();
    Stif = &a->SpStif();
    Mass = &a->SpMass();
    pardiso_solver = std::make_shared<PardisoSolver>();
    pardiso_solver->SetMatrixType(PardisoSolver::REAL_NONSYMMETRIC);
    pardiso_solver->SetOperator(*Stif);
}

AnsysSolver::AnsysSolver(FloatArray &rhoInput, AnsysLoad &ansysInput, std::function<double(double x)> StifRho_,
                         std::function<double(double x)> MassRho_, std::function<double(double x)> DampRho_,
                         std::function<double(double x)> StifRhoD_, std::function<double(double x)> MassRhoD_,
                         std::function<double(double x)> DampRhoD_, std::function<double(double x)> StifRhoD2_,
                         std::function<double(double x)> MassRho2_,
                         std::function<double(double x)> DampRhoD2_) : rho(rhoInput),
                                                                       ansys(ansysInput),
                                                                       StifRho(StifRho_),
                                                                       MassRho(MassRho_),
                                                                       StifRhoD(StifRhoD_),
                                                                       MassRhoD(MassRhoD_),
                                                                       StifRhoD2(StifRhoD2_),
                                                                       MassRhoD2(MassRho2_),
                                                                       DampRho(DampRho_),
                                                                       DampRhoD2(DampRhoD2_),
                                                                       DampRhoD(DampRhoD_) {
    a = std::make_shared<BilinearForm>(&ansys);
    std::unordered_set<int> ConstrainedDofs;
    ansys.getConstrainedDofs(ConstrainedDofs);
    a->assemble(rho, StifRho, MassRho, DampRho);
    a->EliminateEssentialBC(ConstrainedDofs, 0);
    a->Finalize();
    a->convertToMKL();
    Stif = &a->SpStif();
    Mass = &a->SpMass();
    Damp = &a->SpDamp();
    pardiso_solver = std::make_shared<PardisoSolver>();
    pardiso_solver->SetMatrixType(PardisoSolver::REAL_NONSYMMETRIC);
    pardiso_solver->SetOperator(*Stif);
}

void AnsysSolver::EigenProblem(int nev) {
    arpack_solver = std::make_shared<ArpackSolver>(Stif->Size(), nev);
    arpack_solver->SetOption(arpack::which::largest_magnitude);
    arpack_solver->Mult(*pardiso_solver, *Mass);
}

void AnsysSolver::QudraticEigenProblem(int nev, std::string SolutionOpt) {
    arpack_solver = std::make_shared<ArpackSolver>(Stif->Size(), nev);
    arpack_solver->SetOption(arpack::which::largest_imaginary);
    arpack_solver->Mult(*pardiso_solver, *Stif, *Mass, *Damp, SolutionOpt);
}


FloatArray AnsysSolver::getEigenVectors() {
    return arpack_solver->getEigenvectors();
};

FloatArray AnsysSolver::getEigenValues() {
    return arpack_solver->getEigenvalues();
}

FloatArray AnsysSolver::getEigenVector(int i) const {
    return arpack_solver->getEigenvector(i);
}

FloatMatrix AnsysSolver::getEigenvector(IntArray &idx) const {
    return arpack_solver->getEigenvector(idx);
}

FloatMatrix AnsysSolver::getEigenvector(std::vector<int> &idx) const {
    return arpack_solver->getEigenvector(idx);
}

IntArray AnsysSolver::getInfo(std::string name) {
    IntArray res;
    std::ifstream infile(name + "_info.txt");
    if (!infile.is_open()) {
        std::cerr << "无法打开文件: " << name + "_info.txt" << std::endl;
        return {};
    }
    int value;
    while (infile >> value) // 连续读取文件中的数字
    {
        res.append(value);
    }
    infile.close();
    return res;
}

IntArray AnsysSolver::getI(std::string name) {
    IntArray res;
    std::ifstream infile(name + "_row.txt");
    if (!infile.is_open()) {
        std::cerr << "无法打开文件: " << name + "_row.txt" << std::endl;
        return {};
    }
    int value;
    while (infile >> value) // 连续读取文件中的数字
    {
        res.append(value - 1);
    }
    infile.close();
    return res;
}

IntArray AnsysSolver::getJ(std::string name) {
    IntArray res;
    std::ifstream infile(name + "_col.txt");
    if (!infile.is_open()) {
        std::cerr << "无法打开文件: " << name + "_col.txt" << std::endl;
        return {};
    }
    int value;
    while (infile >> value) // 连续读取文件中的数字
    {
        res.append(value - 1);
    }
    infile.close();

    return res;
}

FloatArray AnsysSolver::getValue(std::string name) {
    FloatArray res;
    std::ifstream infile(name + "_val.txt");
    if (!infile.is_open()) {
        std::cerr << "无法打开文件: " << name + "_val.txt" << std::endl;
        return {};
    }
    double value;
    while (infile >> value) {
        res.append(value);
    }
    infile.close();
    return res;
}

SparseMatrix *AnsysSolver::createSparseMatrix(std::string name) {
    IntArray II = getI(name);
    IntArray JJ = getJ(name);
    IntArray INFO = getInfo(name);
    FloatArray VAL = getValue(name);

    int *ii = new int[II.giveSize()];
    int *jj = new int[JJ.giveSize()];
    double *val = new double[VAL.giveSize()];
    std::copy(II.begin(), II.end(), ii);
    std::copy(JJ.begin(), JJ.end(), jj);
    std::copy(VAL.begin(), VAL.end(), val);

    SparseMatrix *sparse_matrix = new SparseMatrix(ii, jj, val, INFO[0], INFO[0], true);

    return sparse_matrix;
}

void AnsysSolver::TestEigenProblem() {
    QudraticEigenProblem(10, "Right");

    std::ofstream infile1("StifTest.txt");
    Stif->PrintCSR(infile1);

    infile1.close();
    std::ofstream infile2("MassTest.txt");
    Mass->PrintCSR(infile2);

    infile2.close();
    std::ofstream infile3("DampTest.txt");

    Damp->PrintCSR(infile3);
    infile3.close();

    std::ofstream infile5("ComplexEigenValue.txt");
    infile5 << arpack_solver->getComplexEigenvalues() << std::endl;
    infile5.close();


    std::ofstream infile("ComplexEigenVector.txt");
    if (!infile.is_open()) {
        std::cerr << "无法打开文件: " << "Vector.txt" << std::endl;
    }
    infile << arpack_solver->getComplexEigenvectors() << std::endl;
    infile.close();
}

void AnsysSolver::NelsonRepeat(ComplexScalar &EigenValue, ComplexMatrix &EigenVector, SparseMatrixComplex &STIF,
                               SparseMatrixComplex &MASS,
                               SparseMatrixComplex &DAMP, int index, ComplexVector &DeriviateEigenValue,
                               ComplexMatrix &DeriviateEigenVector) {
    ComplexScalar two(2.0, 0.0);
    ComplexScalar none(-1.0, 0.0);
    ComplexScalar nhalf(-0.5, 0.0);
    auto EigenValue2 = EigenValue * two;
    for (int it = 0; it < 2; ++it) {
        for (int jt = 0; jt < it; ++jt) {
            ComplexVector qj = EigenVector.GetColumn(jt);
            ComplexVector qi = EigenVector.GetColumn(it);
            ComplexVector temp1(STIF.Rows()), temp2(STIF.Rows());
            DAMP.Mult(qi, temp1);
            MASS.Mult(qi, temp2);
            temp2 = temp2 * EigenValue2;
            temp1 = temp2 + temp1;
            auto Rij = qj * temp1;
            qi = qi - qj * Rij;
            EigenVector.SetColumn(it, qi);
        }
        ComplexVector temp1(STIF.Rows()), temp2(STIF.Rows());
        ComplexVector qj = EigenVector.GetColumn(it);
        DAMP.Mult(qj, temp1);
        MASS.Mult(qj, temp2);
        temp2 = temp2 * EigenValue2;
        temp1 = temp2 + temp1;
        auto temp = qj * temp1;
        temp = temp.sqrt();
        EigenVector.SetColumn(it, qj / temp);
    }

    FloatMatrix elementStifD = this->getElementStif(index, this->getStifRhoD(this->getRho(index)));
    FloatMatrix elementMassD = this->getElementMass(index, this->getMassRhoD(this->getRho(index)));
    FloatMatrix elementDampD = this->getElementDamp(index, this->getDampRhoD(this->getRho(index)));
    IntArray matrixLocation = this->getElementLocation(index);
    ComplexMatrix elementStifDC(elementStifD, "Real");
    ComplexMatrix elementMassDC(elementMassD, "Real");
    ComplexMatrix elementDampDC(elementDampD, "Real");
    auto elementEigenVector = EigenVector.GetRow(matrixLocation);
    auto E = elementEigenVector.transpose() * (elementMassDC * EigenValue * EigenValue + elementDampDC * EigenValue +
                                               elementStifDC) * elementEigenVector * none;
    DeriviateEigenValue.resize(E.Rows());
    ComplexMatrix tao(E.Rows(), E.Cols());
    E.ComputeEigenValueAndEigenVector(DeriviateEigenValue, tao);
    EigenVector = EigenVector * tao;

    for (int it = 0; it < 2; ++it) {
        ComplexVector temp1(STIF.Rows()), temp2(STIF.Rows());
        ComplexVector qj = EigenVector.GetColumn(it);
        DAMP.Mult(qj, temp1);
        MASS.Mult(qj, temp2);
        temp2 = temp2 * EigenValue2;
        temp1 = temp2 + temp1;
        auto temp = qj * temp1;
        temp = temp.sqrt();
        EigenVector.SetColumn(it, qj / temp);
    }
    SparseMatrixComplex a1, a2;
    MASS.Add(EigenValue * EigenValue, STIF, a2);
    DAMP.Add(EigenValue, a2, a1);
    ComplexMatrix temp1(EigenVector.Rows(), EigenVector.Cols());
    ComplexMatrix temp2(EigenVector.Rows(), EigenVector.Cols());
    DAMP.Mult(EigenVector, temp1);
    MASS.Mult(EigenVector, temp2);

    temp2 = temp2 * EigenValue2;
    temp1 = temp1 * none - temp2;
    ComplexMatrix DEV(DeriviateEigenValue.Size(), DeriviateEigenValue.Size());
    for (int it = 0; it < 2; ++it) {
        DEV(it, it) = DeriviateEigenValue[it];
    }

    elementEigenVector = EigenVector.GetRow(matrixLocation);

    auto temp3 = (elementMassDC * EigenValue * EigenValue + elementDampDC * EigenValue + elementStifDC) *
                 elementEigenVector * none;

    temp1 = temp1 * DEV;

    temp2.Zero();
    IntArray colindex = {0, 1};
    temp2.SetSubMatrix(matrixLocation, colindex, temp3);

    auto f1 = temp1 + temp2;
    auto g = EigenVector.FindMaxElementLoaction2();

    f1.EliminateElement(g);
    a1.EliminateElement(g);

    PardisoSolver PA;
    PA.SetMatrixType(PardisoSolver::COMPLEX_NONSYMMETRIC);
    PA.SetOperator(a1);

    ComplexMatrix V1(this->ansys.getNumberOfDof(), 2);

    for (int i = 0; i < 2; ++i) {
        ComplexVector DE = f1.GetColumn(i);
        ComplexVector DER(DE.Size());
        PA.Mult(DE, DER);
        for (int it = 0; it < DER.Size(); it++) {
            V1.SetColumn(i, DER);
        }
    }

    ComplexMatrix b(2, 2);


    FloatMatrix elementStifD2 = this->getElementStif(index, this->getStifRhoD2(this->getRho(index)));
    FloatMatrix elementMassD2 = this->getElementMass(index, this->getMassRhoD2(this->getRho(index)));
    FloatMatrix elementDampD2 = this->getElementDamp(index, this->getDampRhoD2(this->getRho(index)));
    ComplexMatrix elementStifDC2(elementStifD2, "Real");
    ComplexMatrix elementMassDC2(elementMassD2, "Real");
    ComplexMatrix elementDampDC2(elementDampD2, "Real");
    auto L1 = elementEigenVector.transpose() * (elementMassDC2 * EigenValue * EigenValue + elementDampDC2 * EigenValue +
                                                elementStifDC2) * elementEigenVector * nhalf - elementEigenVector.
              transpose() * (elementMassDC * EigenValue * two + elementDampDC) * elementEigenVector * DEV;

    MASS.Mult(EigenVector, temp1);
    auto V1T = V1.GetRow(matrixLocation);
    auto L2 = EigenVector.transpose() * temp1 * DEV * DEV * none - elementEigenVector.transpose() * (
                  elementMassDC * EigenValue * EigenValue + elementDampDC * EigenValue + elementStifDC) * V1T;
    auto ET = elementEigenVector.transpose() * (elementMassDC * EigenValue * two + elementDampDC) * elementEigenVector *
              nhalf - EigenVector.transpose() * temp1 * DEV;
    MASS.Mult(V1, temp1);
    DAMP.Mult(V1, temp2);
    auto L3 = EigenVector.transpose() * (temp1 * EigenValue * two + temp2) * DEV * none;
    auto L = L1 + L2 + L3;
    auto VT = (temp1 * EigenValue * two + temp2).transpose() * EigenVector;


    for (int ii = 0; ii < 2; ii++) {
        for (int jj = 0; jj < 2; jj++) {
            if (ii != jj) {
                b.Set(ii, jj, L.At(ii, jj) / (DEV.At(jj, jj) - DEV.At(ii, ii)));
            } else {
                b.Set(ii, jj, ET.At(ii, ii) - VT.At(ii, ii));
            }
        }
    }

    DeriviateEigenVector = V1 + EigenVector * b;
}


void AnsysSolver::NelsonSingle(ComplexScalar &EigenValue, ComplexVector &EigenVector, SparseMatrixComplex &STIF,
                               SparseMatrixComplex &MASS, SparseMatrixComplex &DAMP, int index,
                               ComplexScalar &DeriviateEigenValue,
                               ComplexVector &DeriviateEigenVector) {
    ComplexScalar two(2.0, 0.0);
    ComplexScalar none(-1.0, 0.0);
    ComplexScalar nhalf(-0.5, 0.0);

    auto EigenValue2 = EigenValue * two;


    ComplexVector temp1(STIF.Rows()), temp2(STIF.Rows());

    FloatMatrix elementStifD = this->getElementStif(index, this->getStifRhoD(this->getRho(index)));
    FloatMatrix elementMassD = this->getElementMass(index, this->getMassRhoD(this->getRho(index)));
    FloatMatrix elementDampD = this->getElementDamp(index, this->getDampRhoD(this->getRho(index)));
    IntArray matrixLocation = this->getElementLocation(index);
    ComplexMatrix elementStifDC(elementStifD, "Real");
    ComplexMatrix elementMassDC(elementMassD, "Real");
    ComplexMatrix elementDampDC(elementDampD, "Real");
    auto elementEigenVector = EigenVector.GetRow(matrixLocation);
    DeriviateEigenValue = elementEigenVector * ((elementMassDC * EigenValue * EigenValue + elementDampDC * EigenValue +
                                                 elementStifDC) * elementEigenVector) * none;
    SparseMatrixComplex a1, a2;
    MASS.Add(EigenValue * EigenValue, STIF, a2);
    DAMP.Add(EigenValue, a2, a1);
    DAMP.Mult(EigenVector, temp1);
    MASS.Mult(EigenVector, temp2);
    temp2 = temp2 * EigenValue2;
    temp1 = temp1 * none - temp2;
    elementEigenVector = EigenVector.GetRow(matrixLocation);
    auto temp3 = (elementMassDC * EigenValue * EigenValue + elementDampDC * EigenValue + elementStifDC) *
                 elementEigenVector * none;
    temp1 = temp1 * DeriviateEigenValue;
    temp2.Zero();
    for (int i = 0; i < matrixLocation.giveSize(); i++) {
        temp2.Set(matrixLocation(i), temp3(i));
    }
    auto f1 = temp1 + temp2;

    IntArray g = {EigenVector.FindMinElementLoaction()};

    f1.EliminateElement(g);
    a1.EliminateElement(g);

    PardisoSolver PA;
    PA.SetMatrixType(PardisoSolver::COMPLEX_NONSYMMETRIC);
    PA.SetOperator(a1);

    ComplexVector V1(this->ansys.getNumberOfDof());
    PA.Mult(f1, V1);

    auto L1 = elementEigenVector * ((elementMassDC * EigenValue * two + elementDampDC) * elementEigenVector) *
              DeriviateEigenValue * none;

    MASS.Mult(EigenVector, temp1);
    auto V1T = V1.GetRow(matrixLocation);
    auto L2 = EigenVector * temp1 * DeriviateEigenValue * DeriviateEigenValue * none -
              elementEigenVector * ((elementMassDC * EigenValue * EigenValue + elementDampDC * EigenValue +
                                     elementStifDC) * V1T);
    auto ET = elementEigenVector * ((elementMassDC * EigenValue * two + elementDampDC) * elementEigenVector) * nhalf -
              EigenVector * temp1 * DeriviateEigenValue;
    MASS.Mult(V1, temp1);
    DAMP.Mult(V1, temp2);
    auto L3 = EigenVector * (temp1 * EigenValue * two + temp2) * DeriviateEigenValue * none;
    auto L = L1 + L2 + L3;
    auto VT = (temp1 * EigenValue * two + temp2) * EigenVector;
    auto b = ET - VT;
    DeriviateEigenVector = V1 + EigenVector * b;
}

void AnsysSolver::NelsonSinglePrepareFactorized(
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
    ComplexVector &RhsCoeff)
{
    ComplexScalar two(2.0, 0.0);
    ComplexScalar none(-1.0, 0.0);

    auto EigenValue2 = EigenValue * two;

    SparseMatrixComplex a2;

    MASS.Add(EigenValue * EigenValue, STIF, a2);
    DAMP.Add(EigenValue, a2, A1);

    g = {EigenVector.FindMinElementLoaction()};

    A1.EliminateElement(g);

    PA.SetMatrixType(PardisoSolver::COMPLEX_NONSYMMETRIC);
    PA.SetOperator(A1);

    DAMP.Mult(EigenVector, Cphi);
    MASS.Mult(EigenVector, Mphi);

    RhsCoeff = Cphi * none - Mphi * EigenValue2;
}

void AnsysSolver::NelsonSingleBackSubFactorized(
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
    ComplexVector &DeriviateEigenVector)
{
    ComplexScalar two(2.0, 0.0);
    ComplexScalar none(-1.0, 0.0);
    ComplexScalar nhalf(-0.5, 0.0);

    auto EigenValue2 = EigenValue * two;

    ComplexVector temp1(STIF.Rows());
    ComplexVector temp2(STIF.Rows());

    int rhoId = this->getRho(index);

    FloatMatrix elementStifD =
        this->getElementStif(index, this->getStifRhoD(rhoId));

    FloatMatrix elementMassD =
        this->getElementMass(index, this->getMassRhoD(rhoId));

    FloatMatrix elementDampD =
        this->getElementDamp(index, this->getDampRhoD(rhoId));

    IntArray matrixLocation =
        this->getElementLocation(index);

    ComplexMatrix elementStifDC(elementStifD, "Real");
    ComplexMatrix elementMassDC(elementMassD, "Real");
    ComplexMatrix elementDampDC(elementDampD, "Real");

    auto elementEigenVector =
        EigenVector.GetRow(matrixLocation);

    auto elementDynamicDC =
        elementMassDC * EigenValue * EigenValue
      + elementDampDC * EigenValue
      + elementStifDC;

    DeriviateEigenValue =
        elementEigenVector
      * (elementDynamicDC * elementEigenVector)
      * none;

    auto temp3 =
        elementDynamicDC * elementEigenVector * none;

    temp1 = RhsCoeff * DeriviateEigenValue;

    temp2.Zero();

    for (int i = 0; i < matrixLocation.giveSize(); i++) {
        temp2.Set(matrixLocation(i), temp3(i));
    }

    auto f1 = temp1 + temp2;

    f1.EliminateElement(g);

    ComplexVector V1(this->ansys.getNumberOfDof());

    // 这里只做回代，不再重新分解
    PA.Mult(f1, V1);

    auto L1 =
        elementEigenVector
      * ((elementMassDC * EigenValue * two + elementDampDC)
      * elementEigenVector)
      * DeriviateEigenValue
      * none;

    auto V1T =
        V1.GetRow(matrixLocation);

    auto L2 =
        EigenVector * Mphi
      * DeriviateEigenValue
      * DeriviateEigenValue
      * none
      - elementEigenVector * (elementDynamicDC * V1T);

    auto ET =
        elementEigenVector
      * ((elementMassDC * EigenValue * two + elementDampDC)
      * elementEigenVector)
      * nhalf
      - EigenVector * Mphi * DeriviateEigenValue;

    MASS.Mult(V1, temp1);
    DAMP.Mult(V1, temp2);

    auto L3 =
        EigenVector
      * (temp1 * EigenValue * two + temp2)
      * DeriviateEigenValue
      * none;

    auto L = L1 + L2 + L3;

    auto VT =
        (temp1 * EigenValue * two + temp2)
      * EigenVector;

    auto b = ET - VT;

    DeriviateEigenVector =
        V1 + EigenVector * b;
}


void AnsysSolver::NelsonAdd(ComplexScalar &EigenValue, ComplexVector &EigenVector, SparseMatrixComplex &STIF,
                            SparseMatrixComplex &MASS, SparseMatrixComplex &DAMP, int index,
                            ComplexScalar &DeriviateEigenValue, ComplexVector &DeriviateEigenVector) {
    ComplexScalar two(2.0, 0.0);
    ComplexScalar none(-1.0, 0.0);
    ComplexScalar nhalf(-0.5, 0.0);
    ComplexScalar zero(0.0, 0.0);

    auto EigenValue2 = EigenValue * two;
    ComplexVector temp1(STIF.Rows()), temp2(STIF.Rows());
    FloatMatrix elementStifD = this->getElementStif(index, this->getStifRhoD(this->getRho(index)));
    FloatMatrix elementMassD = this->getElementMass(index, this->getMassRhoD(this->getRho(index)));
    FloatMatrix elementDampD = this->getElementDamp(index, this->getDampRhoD(this->getRho(index)));
    IntArray matrixLocation = this->getElementLocation(index);
    ComplexMatrix elementStifDC(elementStifD, "Real");
    ComplexMatrix elementMassDC(elementMassD, "Real");
    ComplexMatrix elementDampDC(elementDampD, "Real");
    auto elementEigenVector = EigenVector.GetRow(matrixLocation);
    DeriviateEigenValue = elementEigenVector * ((elementMassDC * EigenValue * EigenValue + elementDampDC * EigenValue +
                                                 elementStifDC) * elementEigenVector) * none;
    SparseMatrixComplex a1, a2;
    MASS.Add(EigenValue * EigenValue, STIF, a2);
    DAMP.Add(EigenValue, a2, a1);
    DAMP.Mult(EigenVector, temp1);
    MASS.Mult(EigenVector, temp2);
    temp2 = temp2 * EigenValue2;
    temp1 = temp1 * none - temp2;
    auto a3 = temp1 * none;
    elementEigenVector = EigenVector.GetRow(matrixLocation);
    auto temp3 = (elementMassDC * EigenValue * EigenValue + elementDampDC * EigenValue + elementStifDC) *
                 elementEigenVector * none;
    temp1 = temp1 * DeriviateEigenValue;
    temp2.Zero();
    for (int i = 0; i < matrixLocation.giveSize(); i++) {
        temp2.Set(matrixLocation(i), temp3(i));
    }
    auto f1 = temp1 + temp2;

    ComplexVector F(f1.Size() + 1);
    for (int i = 0; i < f1.Size(); i++) {
        F(i) = f1(i);
    }


    a1.Augment(a3, zero);
    PardisoSolver PA;
    PA.SetMatrixType(PardisoSolver::COMPLEX_NONSYMMETRIC);
    PA.SetOperator(a1);
    ComplexVector Result(f1.Size() + 1);
    ComplexVector V1(f1.Size());
    PA.Mult(F, Result);
    for (int i = 0; i < f1.Size(); i++) {
        V1(i) = Result(i);
    }
    auto L1 = elementEigenVector * ((elementMassDC * EigenValue * two + elementDampDC) * elementEigenVector) *
              DeriviateEigenValue * none;

    MASS.Mult(EigenVector, temp1);
    auto V1T = V1.GetRow(matrixLocation);
    auto L2 = EigenVector * temp1 * DeriviateEigenValue * DeriviateEigenValue * none -
              elementEigenVector * ((elementMassDC * EigenValue * EigenValue + elementDampDC * EigenValue +
                                     elementStifDC) * V1T);
    auto ET = elementEigenVector * ((elementMassDC * EigenValue * two + elementDampDC) * elementEigenVector) * nhalf -
              EigenVector * temp1 * DeriviateEigenValue;
    MASS.Mult(V1, temp1);
    DAMP.Mult(V1, temp2);
    auto L3 = EigenVector * (temp1 * EigenValue * two + temp2) * DeriviateEigenValue * none;
    auto L = L1 + L2 + L3;
    auto VT = (temp1 * EigenValue * two + temp2) * EigenVector;
    auto b = ET - VT;
    DeriviateEigenVector = V1 + EigenVector * b;
}

void AnsysSolver::NelsonAddPrepareFactorized(
    ComplexScalar &EigenValue,
    ComplexVector &EigenVector,
    SparseMatrixComplex &STIF,
    SparseMatrixComplex &MASS,
    SparseMatrixComplex &DAMP,
    SparseMatrixComplex &Aadd,
    PardisoSolver &PA,
    ComplexVector &Mphi,
    ComplexVector &RhsCoeff
)
{
    ComplexScalar two(2.0, 0.0);
    ComplexScalar none(-1.0, 0.0);
    ComplexScalar zero(0.0, 0.0);

    auto EigenValue2 = EigenValue * two;

    ComplexVector Cphi(STIF.Rows());

    SparseMatrixComplex a2;

    // a1 = STIF + lambda * DAMP + lambda^2 * MASS
    MASS.Add(EigenValue * EigenValue, STIF, a2);
    DAMP.Add(EigenValue, a2, Aadd);

    // Cphi = C * phi
    DAMP.Mult(EigenVector, Cphi);

    // Mphi = M * phi
    MASS.Mult(EigenVector, Mphi);

    // 原代码：
    // temp2 = Mphi * EigenValue2;
    // temp1 = Cphi * none - temp2;
    //
    // 因此：
    // RhsCoeff = -Cphi - 2 * lambda * Mphi
    RhsCoeff = Cphi * none - Mphi * EigenValue2;

    // 原代码：
    // auto a3 = temp1 * none;
    //
    // 因此：
    // a3 = Cphi + 2 * lambda * Mphi
    auto a3 = RhsCoeff * none;

    // 原代码：
    // a1.Augment(a3, zero);
    //
    // 这里把增广矩阵构造放到预处理阶段
    Aadd.Augment(a3, zero);

    PA.SetMatrixType(PardisoSolver::COMPLEX_NONSYMMETRIC);

    // 关键：只分解一次
    PA.SetOperator(Aadd);
}

void AnsysSolver::NelsonAddBackSubFactorized(
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
)
{
    ComplexScalar two(2.0, 0.0);
    ComplexScalar none(-1.0, 0.0);
    ComplexScalar nhalf(-0.5, 0.0);

    const int ndof = STIF.Rows();

    ComplexVector temp1(ndof);
    ComplexVector temp2(ndof);

    int rhoId = this->getRho(index);

    FloatMatrix elementStifD =
        this->getElementStif(index, this->getStifRhoD(rhoId));

    FloatMatrix elementMassD =
        this->getElementMass(index, this->getMassRhoD(rhoId));

    FloatMatrix elementDampD =
        this->getElementDamp(index, this->getDampRhoD(rhoId));

    IntArray matrixLocation =
        this->getElementLocation(index);

    ComplexMatrix elementStifDC(elementStifD, "Real");
    ComplexMatrix elementMassDC(elementMassD, "Real");
    ComplexMatrix elementDampDC(elementDampD, "Real");

    auto elementEigenVector =
        EigenVector.GetRow(matrixLocation);

    auto elementDynamicDC =
        elementMassDC * EigenValue * EigenValue
      + elementDampDC * EigenValue
      + elementStifDC;

    // 原代码：
    // DeriviateEigenValue =
    //     elementEigenVector * (elementDynamicDC * elementEigenVector) * none;
    DeriviateEigenValue =
        elementEigenVector
      * (elementDynamicDC * elementEigenVector)
      * none;

    // 原代码：
    // temp3 = elementDynamicDC * elementEigenVector * none;
    auto temp3 =
        elementDynamicDC
      * elementEigenVector
      * none;

    // 原代码：
    // temp1 = temp1 * DeriviateEigenValue;
    //
    // 这里 temp1 对应原来的：
    // temp1 = (-Cphi - 2 lambda Mphi) * DeriviateEigenValue
    temp1 = RhsCoeff * DeriviateEigenValue;

    // 原代码：
    // temp2.Zero();
    // for (...) temp2.Set(matrixLocation(i), temp3(i));
    temp2.Zero();

    for (int i = 0; i < matrixLocation.giveSize(); i++) {
        temp2.Set(matrixLocation(i), temp3(i));
    }

    // 原代码：
    // auto f1 = temp1 + temp2;
    auto f1 = temp1 + temp2;

    // 原代码：
    // ComplexVector F(f1.Size() + 1);
    // for (int i = 0; i < f1.Size(); i++) {
    //     F(i) = f1(i);
    // }
    //
    // 注意最后一个增广自由度右端项默认为 0
    ComplexVector F(f1.Size() + 1);
    F.Zero();

    for (int i = 0; i < f1.Size(); i++) {
        F(i) = f1(i);
    }

    ComplexVector Result(f1.Size() + 1);
    ComplexVector V1(f1.Size());

    // 关键：这里只回代，不再 SetOperator
    PA.Mult(F, Result);

    for (int i = 0; i < f1.Size(); i++) {
        V1(i) = Result(i);
    }

    // ------------------------------------------------------------------
    // 以下部分保持你原 NelsonAdd() 的后处理结构
    // ------------------------------------------------------------------

    auto L1 =
        elementEigenVector
      * ((elementMassDC * EigenValue * two + elementDampDC)
      * elementEigenVector)
      * DeriviateEigenValue
      * none;

    auto V1T =
        V1.GetRow(matrixLocation);

    auto L2 =
        EigenVector
      * Mphi
      * DeriviateEigenValue
      * DeriviateEigenValue
      * none
      - elementEigenVector
      * (elementDynamicDC * V1T);

    auto ET =
        elementEigenVector
      * ((elementMassDC * EigenValue * two + elementDampDC)
      * elementEigenVector)
      * nhalf
      - EigenVector
      * Mphi
      * DeriviateEigenValue;

    MASS.Mult(V1, temp1);
    DAMP.Mult(V1, temp2);

    auto L3 =
        EigenVector
      * (temp1 * EigenValue * two + temp2)
      * DeriviateEigenValue
      * none;

    auto L = L1 + L2 + L3;

    // 如果编译器提示 L 未使用，可以保留这一句
    (void)L;

    auto VT =
        (temp1 * EigenValue * two + temp2)
      * EigenVector;

    auto b = ET - VT;

    DeriviateEigenVector =
        V1 + EigenVector * b;
}

void AnsysSolver::Normalize(const ComplexScalar &EigenValue, ComplexVector &EigenVector,
                            const SparseMatrixComplex &DAMP, const SparseMatrixComplex &MASS) {
    auto EigenValue2 = EigenValue * ComplexScalar{2.0, 0.0};
    ComplexVector temp1(DAMP.Rows()), temp2(DAMP.Rows());
    DAMP.Mult(EigenVector, temp1);
    MASS.Mult(EigenVector, temp2);
    temp2 = temp2 * EigenValue2;
    temp1 = temp2 + temp1;
    auto temp = EigenVector * temp1;
    temp = temp.sqrt();
    EigenVector = EigenVector / temp;
}

ComplexMatrix AnsysSolver::AddNoiseToEigenvectors(const ComplexMatrix &EigenVector,double NoiseLevel) {
    std::mt19937 gen(1024);
    std::normal_distribution<double> dist(0.0, 1.0);
    int current_n = EigenVector.Rows();
    int current_m = EigenVector.Cols();
    ComplexMatrix eigenvector_measured = EigenVector;
    for (int k = 0; k < current_m; ++k) {
        // --- 第一步：获取该阶模态的最大幅值 (amp_max) ---
        double amp_max = 0.0;
        for (int i = 0; i < current_n; ++i) {
            // 假设你的矩阵支持 (row, col) 操作符访问
            // 如果不支持，可能需要用 GetColumn(k)[i]
            ComplexScalar val = EigenVector(i, k);
            double mag = val.norm();
            if (mag > amp_max) {
                amp_max = mag;
            }
        }

        // --- 第二步：计算噪声并叠加 ---
        for (int i = 0; i < current_n; ++i) {
            // 获取理论值
            ComplexScalar phi_exact = EigenVector(i, k);

            // 生成随机噪声 randn
            double noise_matrix_real = dist(gen);
            double noise_matrix_imag = dist(gen);

            // 计算噪声分量
            double n_real = NoiseLevel * amp_max * noise_matrix_real;
            double n_imag = NoiseLevel * amp_max * noise_matrix_imag;

            // 叠加: 实部+噪声, 虚部+噪声
            double measured_real = phi_exact.real() + n_real;
            double measured_imag = phi_exact.imag() + n_imag;

            // 存入结果矩阵
            // 假设 ComplexScalar 构造函数兼容 (double, double)
            eigenvector_measured(i, k) = ComplexScalar(measured_real, measured_imag);
        }

    }
    return eigenvector_measured;
}

void AnsysSolver::AdjointNelson(ComplexScalar &EigenValue, ComplexVector &EigenVector, SparseMatrixComplex &STIF,
                                SparseMatrixComplex &MASS, SparseMatrixComplex &DAMP, ComplexScalar Dlembda,
                                ComplexVector &Dsymvector,
                                ComplexScalar &alpha, ComplexVector &V) {
    ComplexScalar two(2.0, 0.0);
    ComplexScalar none(-1.0, 0.0);

    auto EigenValue2 = EigenValue * two;

    alpha = Dsymvector * EigenVector * none;

    ComplexVector temp1(STIF.Rows()), temp2(STIF.Rows()), temp3(STIF.Rows());
    SparseMatrixComplex a1, a2;
    MASS.Add(EigenValue * EigenValue, STIF, a2);
    DAMP.Add(EigenValue, a2, a1);
    DAMP.Mult(EigenVector, temp1);
    MASS.Mult(EigenVector, temp3);
    temp2 = temp3 * EigenValue2;
    temp1 = temp1 * none - temp2;
    auto f1 = temp1 * (ComplexScalar{2.0, 0.0} * alpha) - Dsymvector * ComplexScalar{2.0, 0.0};

    IntArray g = {f1.FindMaxElementLoaction()};
    f1.EliminateElement(g);
    a1.EliminateElement(g);

    PardisoSolver PA;
    PA.SetMatrixType(PardisoSolver::COMPLEX_NONSYMMETRIC);
    PA.SetOperator(a1);

    ComplexVector V1(this->ansys.getNumberOfDof());
    PA.Mult(f1, V1);

    auto b1 = temp1 * V1 - EigenVector * temp3 * two * alpha - Dlembda * two;
    V = V1 + EigenVector * b1;
}


void AnsysSolver::AdjointNelosnAdd(ComplexScalar &EigenValue, ComplexVector &EigenVector, SparseMatrixComplex &STIF,
                                   SparseMatrixComplex &MASS, SparseMatrixComplex &DAMP, ComplexScalar Dlembda,
                                   ComplexVector &Dsymvector,
                                   ComplexScalar &alpha, ComplexVector &V) {
    ComplexScalar two(2.0, 0.0);
    ComplexScalar none(-1.0, 0.0);

    auto EigenValue2 = EigenValue * two;

    SparseMatrixComplex a1, a2;
    MASS.Add(EigenValue * EigenValue, STIF, a2);
    DAMP.Add(EigenValue, a2, a1);

    ComplexVector temp1(STIF.Rows()), temp2(STIF.Rows());
    DAMP.Mult(EigenVector, temp1);
    MASS.Mult(EigenVector, temp2);
    auto a4 = EigenVector * temp2;
    temp2 = temp2 * EigenValue2;
    temp1 = temp1 * none - temp2;
    auto a3 = temp1 * none;
    a1.Augment(a3, a4);
    ComplexVector b(STIF.Rows() + 1);
    b.Set(STIF.Rows(), Dlembda * two * none);
    for (int i = 0; i < STIF.Rows(); i++) {
        b.Set(i, Dsymvector[i] * two * none);
    }
    PardisoSolver PA;
    PA.SetMatrixType(PardisoSolver::COMPLEX_NONSYMMETRIC);
    PA.SetOperator(a1);

    ComplexVector Result(b.Size());

    PA.Mult(b, Result);
    for (int i = 0; i < this->ansys.getNumberOfDof(); i++) {
        V.Set(i, Result[i]);
    }
    alpha = Result[this->ansys.getNumberOfDof()] / two;
}

void AnsysSolver::AdjointsymmDsymvector(ComplexScalar &EigenValue, ComplexVector &EigenVector,
                                        SparseMatrixComplex &STIF,
                                        SparseMatrixComplex &MASS, SparseMatrixComplex &DAMP, const PardisoSolver &PA,
                                        ComplexVector &Dsymvector,
                                        ComplexVector &Dphi, ComplexScalar &DphiL, ComplexScalar &EigenValueSquare,
                                        ComplexScalar &Db3Dphib0,
                                        const double tol, const int maxit) {
    ComplexScalar two(2.0, 0.0);
    ComplexScalar none(-1.0, 0.0);
    ComplexScalar nhalf(-0.5, 0.0);
    ComplexScalar zero(0.0, 0.0);

    auto EigenValue2 = EigenValue * two;
    ComplexVector Li(STIF.Rows()), temp(STIF.Rows()), t(STIF.Rows());
    DAMP.Mult(EigenVector, Li);
    MASS.Mult(EigenVector, temp);
    auto Mphi = temp;
    temp = temp * EigenValue2;
    Li = Li + temp;


    ComplexVector x(STIF.Rows()), Aq_sparse(STIF.Rows()), Aq(STIF.Rows());
    SparseMatrixComplex A_sparse, a2;
    MASS.Add(EigenValue * EigenValue, STIF, a2);
    DAMP.Add(EigenValue, a2, A_sparse);

    auto r = Dsymvector;


    PA.Mult(r, t);

    auto tao = t.Norm2();
    auto q = t;
    ComplexVector d(STIF.Rows());
    double v_old = 0.0, v_new = 0.0, c = 0.0;
    auto rhoSQMR = r * q;

    for (int i = 0; i < maxit; i++) {
        A_sparse.Mult(q, Aq_sparse);
        auto scalar_alpha = Li * q;
        Aq = Aq_sparse + Li * scalar_alpha;
        auto sigma = q * Aq;
        if (sigma.norm() < 1e-20) {
            std::cerr << "WARNING: sigma is too small" << std::endl;
        }
        auto alpha = rhoSQMR / sigma;
        r = r - Aq * alpha;
        PA.Mult(r, t);
        auto norm_t = t.Norm2();
        v_new = norm_t / tao;
        c = 1 / sqrt(1 + v_new * v_new);
        tao = tao * v_new * c;
        auto c_sq = c * c;
        auto scaler_d = c_sq * v_old * v_old;
        auto scaler_q = alpha * ComplexScalar{c_sq, 0.0};
        d = d * scaler_d + q * scaler_q;
        x = x + d;
        v_old = v_new;
        if (d.Norm2() < tol) {
            break;
        }
        auto rho_new = r * t;
        if (rho_new.norm() < 1e-15) {
            break;
        }
        auto beta = rho_new / rhoSQMR;
        q = t + q * beta;
        rhoSQMR = rho_new;
    }
    Dphi = x;

    Li = Li * nhalf;
    ComplexVector temp1(STIF.Rows()), temp2(STIF.Rows());
    DAMP.Mult(EigenVector, temp1);
    MASS.Mult(EigenVector, temp2);
    temp2 = temp2 * EigenValue2;
    temp1 = temp1 + temp2;
    DphiL = Dphi * Li;

    auto Db3 = DphiL * (EigenVector * Mphi) * 2;

    auto Dphib0 = Dphi * temp1;
    Db3Dphib0 = Db3 - Dphib0;
    EigenValueSquare = EigenValue * EigenValue;
}

double AnsysSolver::MIR(ComplexVector &EigenVector) {
    FloatArray IMAG(EigenVector.Size());
    double mir = 0.0;
#pragma omp parallel for
    for (int i = 0; i < EigenVector.Size(); i++) {
        IMAG(i) = EigenVector(i).imag;
    }
    return IMAG.computeNorm() / EigenVector.Norm2();
}

double AnsysSolver::MPC(ComplexVector &EigenVector) {
    FloatArray IMAG(EigenVector.Size());
    FloatArray REAL(EigenVector.Size());
    double mpc = 0.0;
#pragma omp parallel for
    for (int i = 0; i < EigenVector.Size(); i++) {
        IMAG(i) = EigenVector(i).imag;
        REAL(i) = EigenVector(i).real;
    }
    double SXY, SXX, SYY;
    SXY = REAL * IMAG;
    SXX = REAL * REAL;
    SYY = IMAG * IMAG;
    mpc = (SXX - SYY) * (SXX - SYY) + 4 * SXY * SXY;
    mpc /= (SXX + SYY) * (SXX + SYY);
    return mpc;
}

double AnsysSolver::MAC(ComplexVector &EigenVector_Measured, ComplexVector &EigenVector) {
    auto TNN=EigenVector_Measured.conjDot(EigenVector);
    auto NN=TNN.real()*TNN.real()+TNN.imag()*TNN.imag();
    auto DI=pow(EigenVector_Measured.Norm2(),2);
    auto DJ=pow(EigenVector.Norm2(),2);
    return NN/(DI*DJ);
}


ComplexVector AnsysSolver::DMIR(ComplexVector &EigenVector) {
    ComplexVector RES(EigenVector.Size());
    FloatArray IMAG(EigenVector.Size()), TEMP(EigenVector.Size());
#pragma omp parallel for
    for (int i = 0; i < EigenVector.Size(); i++) {
        IMAG(i) = EigenVector(i).imag;
    }
    double NN = IMAG.computeNorm();
    double DD = EigenVector.Norm2();
    double LEFT = -NN / DD / DD / DD / 2;
    double RIGHT = -1 / NN / DD / 2;
    ComplexScalar left(LEFT);

    RES = EigenVector.Conjugate() * left;
#pragma omp parallel for
    for (int i = 0; i < EigenVector.Size(); i++) {
        RES(i).imag += EigenVector(i).imag * RIGHT;
    }
    return RES;
}

ComplexVector AnsysSolver::DMPC(ComplexVector &EigenVector) {
    ComplexVector RES(EigenVector.Size());
    FloatArray IMAG(EigenVector.Size());
    FloatArray REAL(EigenVector.Size());
    double mpc = 0.0;
#pragma omp parallel for
    for (int i = 0; i < EigenVector.Size(); i++) {
        IMAG(i) = EigenVector(i).imag;
        REAL(i) = EigenVector(i).real;
    }
    double SXY, SXX, SYY;
    SXY = REAL * IMAG;
    SXX = REAL * REAL;
    SYY = IMAG * IMAG;
    mpc = (SXX - SYY) * (SXX - SYY) + 4 * SXY * SXY;
    mpc /= (SXX + SYY) * (SXX + SYY);
    double Delta = SXX - SYY;
    double SUM = SXX + SYY;
    double TempRealA = Delta - mpc * SUM;
    double SXY2 = 2 * SXY;
    double TempImagB = Delta + mpc * SUM;
    double cof = 2 / SUM / SUM;
#pragma omp parallel for
    for (int i = 0; i < EigenVector.Size(); i++) {
        RES.Set(i, cof * TempRealA * EigenVector(i).real + cof * SXY2 * EigenVector(i).imag,
                -cof * SXY2 * EigenVector(i).real + cof * TempImagB * EigenVector(i).imag);
    }
    return RES;
}

ComplexVector AnsysSolver::DMAC(ComplexVector &EigenVector_Measured, ComplexVector &EigenVector) {
    auto TNN=EigenVector_Measured.conjDot(EigenVector);
    auto NN=TNN.real()*TNN.real()+TNN.imag()*TNN.imag();
    auto DI=pow(EigenVector_Measured.Norm2(),2);
    auto DJ=pow(EigenVector.Norm2(),2);
    double mac=NN/(DI*DJ);
    ComplexScalar Cmac(mac,0.0);
    auto ME=EigenVector_Measured.conjDot(EigenVector);
    auto EE=EigenVector.conjDot(EigenVector);
    auto CEigenVector_Measured=EigenVector_Measured.Conjugate();
    auto CEigenVector=EigenVector.Conjugate();
    CEigenVector_Measured=CEigenVector_Measured/ME;
    CEigenVector=CEigenVector/EE;
    return (CEigenVector_Measured-CEigenVector)*Cmac;
}


ComplexScalar AnsysSolver::AdjointElement(ComplexScalar &EigenValue, ComplexVector &EigenVector,
                                          const ComplexScalar &alpha,
                                          ComplexVector &V, const int index) {
    ComplexScalar two(2.0, 0.0);

    auto EigenValue2 = EigenValue * two;
    int rhoId = this->getRho(index);

    FloatMatrix elementStifD = this->getElementStif(index, this->getStifRhoD(rhoId));
    FloatMatrix elementMassD = this->getElementMass(index, this->getMassRhoD(rhoId));
    FloatMatrix elementDampD = this->getElementDamp(index, this->getDampRhoD(rhoId));
    IntArray matrixLocation = this->getElementLocation(index);
    ComplexMatrix elementStifDC(elementStifD, "Real");
    ComplexMatrix elementMassDC(elementMassD, "Real");
    ComplexMatrix elementDampDC(elementDampD, "Real");
    auto elementEigenVector = EigenVector.GetRow(matrixLocation);
    auto elementV = V.GetRow(matrixLocation);

    auto term1 = elementV * ((elementMassDC * EigenValue * EigenValue + elementDampDC * EigenValue + elementStifDC) *
                             elementEigenVector);
    auto term2 = elementEigenVector * ((elementMassDC * EigenValue2 + elementDampDC) * elementEigenVector) * alpha;

    return term1 + term2;
}

ComplexScalar AnsysSolver::AdjointsymmElement(
    ComplexScalar &EigenValue,
    ComplexVector &EigenVector,
    ComplexVector &Dphi,
    ComplexScalar &DphiL,
    const ComplexScalar &EigenValueSquare,
    const ComplexScalar &Db3Dphib0,
    const ComplexScalar &Dlemda,
    const int index)
{
    ComplexScalar two(2.0, 0.0);
    ComplexScalar none(-1.0, 0.0);

    auto EigenValue2 = EigenValue * two;

    int rhoId = this->getRho(index);

    // 注意：这里不要加 const
    IntArray matrixLocation = this->getElementLocation(index);

    ComplexMatrix elementStifDC(
        this->getElementStif(index, this->getStifRhoD(rhoId)), "Real");

    ComplexMatrix elementMassDC(
        this->getElementMass(index, this->getMassRhoD(rhoId)), "Real");

    ComplexMatrix elementDampDC(
        this->getElementDamp(index, this->getDampRhoD(rhoId)), "Real");

    auto elementEigenVector = EigenVector.GetRow(matrixLocation);
    auto elementDphi = Dphi.GetRow(matrixLocation);

    auto Mphi = elementMassDC * elementEigenVector;
    auto Cphi = elementDampDC * elementEigenVector;
    auto Kphi = elementStifDC * elementEigenVector;

    auto b1 = (Mphi * EigenValueSquare + Cphi * EigenValue + Kphi) * none;

    auto b2Vec = Mphi * EigenValue2 + Cphi;

    auto lemdaD = elementEigenVector * b1;
    auto b2 = elementEigenVector * b2Vec;
    auto term1 = elementDphi * b1;
    auto term2 = lemdaD * Db3Dphib0;
    auto term3 = DphiL * b2;

    return (term1 + term2 + term3) * two + lemdaD * Dlemda * two;
}


void AnsysSolver::Nelson2(ComplexScalar &EigenValue, ComplexMatrix &EigenVector, const SparseMatrix &STIF,
                          const SparseMatrix &MASS,
                          const SparseMatrix &DAMP, const SparseMatrix &STIFD, const SparseMatrix &MASSD,
                          const SparseMatrix &DAMPD,
                          const SparseMatrix &STIFD2, const SparseMatrix &MASSD2, const SparseMatrix &DAMPD2) {
    SparseMatrixComplex STIFC(STIF, "Real");
    SparseMatrixComplex MASSC(MASS, "Real");
    SparseMatrixComplex DAMPC(DAMP, "Real");
    SparseMatrixComplex STIFDC(STIFD, "Real");
    SparseMatrixComplex MASSDC(MASSD, "Real");
    SparseMatrixComplex DAMPDC(DAMPD, "Real");
    SparseMatrixComplex STIFD2C(STIFD2, "Real");
    SparseMatrixComplex MASSD2C(MASSD2, "Real");
    SparseMatrixComplex DAMPD2C(DAMPD2, "Real");
    STIFC.covertToMKLCSR();
    MASSC.covertToMKLCSR();
    DAMPC.covertToMKLCSR();
    STIFDC.covertToMKLCSR();
    MASSDC.covertToMKLCSR();
    DAMPDC.covertToMKLCSR();
    STIFD2C.covertToMKLCSR();
    MASSD2C.covertToMKLCSR();
    DAMPD2C.covertToMKLCSR();

    for (int ii = 0; ii < 2; ++ii) {
        ComplexVector temp1(STIF.Size()), temp2(STIF.Size());
        ComplexVector eigenvector(STIF.Size());
        for (int i = 0; i < STIF.Size(); ++i) {
            eigenvector(i) = EigenVector(i, ii);
        }
        DAMPC.Mult(eigenvector, temp1);
        MASSC.Mult(eigenvector, temp2);
        ComplexScalar two(2.0, 0.0);
        auto EigenValue2 = EigenValue * two;
        temp2 = temp2 * EigenValue2;
        temp1 = temp2 + temp1;
        auto temp = eigenvector * temp1;
        temp = temp.sqrt();
        eigenvector = eigenvector / temp;
        for (int i = 0; i < STIF.Size(); ++i) {
            EigenVector(i, ii) = eigenvector(i);
        }
    }
    auto lambda2 = EigenValue * EigenValue;
    SparseMatrixComplex D, dD, d2D;
    ComplexScalar one(1.0, 0.0);
    ComplexScalar none(-1.0, 0.0);
    ComplexScalar two(2.0, 0.0);
    DAMPC.Add(EigenValue, STIFC, D);
    MASSC.Add(lambda2, D, D);
    DAMPDC.Add(EigenValue, STIFDC, dD);
    MASSDC.Add(lambda2, dD, dD);
    DAMPD2C.Add(EigenValue, STIFD2C, d2D);
    MASSD2C.Add(lambda2, d2D, d2D);

    SparseMatrixComplex H, dH;
    MASSC.Add(two * EigenValue, DAMPC, H);
    MASSDC.Add(two * EigenValue, DAMPDC, dH);

    ComplexMatrix temp(EigenVector.Rows(), EigenVector.Cols());
    dD.Mult(EigenVector, temp);
    auto E = (EigenVector.transpose() * temp) * none;
    ComplexVector dar(E.Rows());
    ComplexMatrix grm(E.Rows(), E.Cols());
    E.ComputeEigenValueAndEigenVector(dar, grm);

    ComplexMatrix U = EigenVector * grm;
    ComplexMatrix temp2(EigenVector.Rows(), EigenVector.Cols());
    H.Mult(U, temp2);
    ComplexMatrix DAR(E.Rows(), E.Rows());
    for (int i = 0; i < E.Rows(); ++i) {
        DAR(i, i) = dar(i);
    }
    temp2 = temp2 * DAR;
    auto G = temp2 + temp;
    G = G * none;
    IntArray row, col;
    G.findNonZero(row, col);
    int k1 = row[0];


    std::cout << G << std::endl;
}

void AnsysSolver::SetEigenvector(int index, double value) {
    arpack_solver->SetEigenvalue(index, value);
}

void AnsysSolver::Reorthogonalization() {
    int nev = arpack_solver->GetNev();
    for (int i = 0; i < nev; ++i) {
        FloatArray EigenVector = arpack_solver->getEigenvector(i);
        /*
        int maxIndex=EigenVector.giveIndexMaxElem()-1;
        double maxValue=abs(EigenVector(maxIndex));
        arpack_solver->DevideValue(i,maxValue);
        */

        FloatArray TempVector;
        TempVector.resize(EigenVector.giveSize());
        Mass->Mult(EigenVector, TempVector);
        double value = EigenVector * TempVector;
        arpack_solver->DevideValue(i, std::sqrt(value));
    }
}


void AnsysSolver::changeEigenvector(IntArray index, FloatMatrix values) {
    arpack_solver->changeEigenvector(index, values);
}


FloatArray Sensitivity::AdjointMethod(double Eval, FloatArray &Evec, FloatArray &ObjDerivative,
                                      FloatArray &EigenValueDerivative,
                                      FloatArray &EigenVectorDerivative, int MaxIter, double rol) const {
    FloatArray DC(ansys_solver->getElementSize());
    FloatMatrix L(ansys_solver->getSize(), 1);
    int N = ansys_solver->getSize();
    FloatArray tempL(N);
    auto Mass = ansys_solver->getMass();
    auto Stif = ansys_solver->getStif();
    Mass->Mult(Evec, tempL);
    for (int i = 0; i < N; i++) {
        L(i, 0) = tempL(i);
    }
    //FloatArray DREVb=EigenVectorDerivative;
    FloatArray DREVx(N), r(N), t(N), q(N);
    DREVx = 0.0;
    r = EigenVectorDerivative;
    ansys_solver->Mult(r, t);
    double tao0 = norm(t);
    double delt0 = r.computeSquaredNorm();
    q = t;
    double v0 = 0.0;
    double rou0 = r * q;
    double delt = r.computeSquaredNorm();
    int ll = 0;
    int chongpings = 1;
    FloatArray DREVd(N);
    DREVd = 0.0;
    while ((ll < MaxIter) && ((delt / delt0) > rol)) {
        FloatArray b1(N);
        FloatArray b12(chongpings);
        FloatArray b2(N);
        FloatArray b3(N);
        b1 = 0.0;
        b12 = 0.0;
        b2 = 0.0;
        b3 = 0.0;
        b12.beTProductOf(L, q);
        b2.beProductOf(L, b12);
        Mass->Mult(q, b1);
        Stif->Mult(q, t);
        t -= Eval * b1;
        t += b2;
        double sigema0 = q * t;
        if (std::abs(sigema0) <= 1.0e-25) {
            break;
        }
        double erfa0 = rou0 / sigema0;
        r -= erfa0 * t;
        ansys_solver->Mult(r, t);
        double v1 = norm(t);
        v1 /= tao0;
        double c = 1.0 / std::sqrt((1.0 + v1 * v1));
        tao0 *= v1 * c;
        DREVd.times(c * c * v0 * v0);
        DREVd += c * c * erfa0 * q;
        DREVx += DREVd;
        if (std::abs(rou0) <= 1.0e-20) {
            break;
        }
        double rou1 = r * t;
        double beit = rou1 / rou0;
        q.times(beit);
        q += t;
        rou0 = rou1;
        v0 = v1;
        delt = r.computeSquaredNorm();
        ll += 1;
    }
    //std::cout<<"SQMR Iteration: "<<ll<<"  Rol: "<<delt / delt0<<std::endl;
    FloatArray MPAT;
    MPAT.beTProductOf(L, DREVx);
    for (int index = 0; index < ansys_solver->getElementSize(); index++) {
        if (ansys_solver->ansys.AnasysElements[index + 1]->getTopStatus() == true) {
            IntArray edof = ansys_solver->getElementLocation(index);
            FloatMatrix elementStifD = ansys_solver->getElementStif(
                index, ansys_solver->getStifRhoD(ansys_solver->getRho(index)));
            FloatMatrix elementMassD = ansys_solver->getElementMass(
                index, ansys_solver->getMassRhoD(ansys_solver->getRho(index)));
            FloatArray elementTemp;
            Evec.GetSubVector(edof, elementTemp);
            FloatArray eKDP = elementStifD * elementTemp;
            FloatArray eMDP = elementMassD * elementTemp;
            FloatArray TM = eKDP - Eval * eMDP;
            double lambdaD = TM * elementTemp;
            double PTeMDP = eMDP * elementTemp;
            FloatArray elementAT;
            DREVx.GetSubVector(edof, elementAT);
            double alpha = -1.0 * TM * elementAT;
            alpha -= (-lambdaD + 0.5 * PTeMDP) * MPAT(0);
            DC[index] = ObjDerivative[index] + lambdaD * EigenValueDerivative[0] + alpha;
        } else {
            DC[index] = 0;
        }
    }
    return DC;
}


void Top::SetTopElement(IntArray set) {
    for (int i = 1; i <= set.giveSize(); i++) {
        if (set[i - 1] == 1) {
            ansys_solver->ansys.AnasysElements[i]->setTop(true);
        } else {
            ansys_solver->ansys.AnasysElements[i]->setTop(false);
        }
    }
}


void Top::CreateDensityFilteringMatrix(double radius) {
    DensityFilteringMatrix = std::make_unique<SparseMatrix>(ansys_solver->getElementSize(),
                                                            ansys_solver->getElementSize());
    HS.resize(ansys_solver->getElementSize());
    auto Centroid = ansys_solver->getCentroid();
    auto KDTree = ansys_solver->GetKDNode(Centroid);
    for (const auto &elem: ansys_solver->ansys.AnasysElements) {
        if (elem.second->getTopStatus() == true) {
            auto RNode = ansys_solver->rangeSearch(KDTree, Centroid[elem.first - 1], radius);
            int iH = elem.first - 1;
            for (const auto &rnode: RNode) {
                int jH = rnode.id - 1;
                double dist = Centroid[elem.first - 1].euclideanDistance(rnode);
                double sH = std::max(0.0, radius - dist);
                HS[iH] += sH;
                DensityFilteringMatrix->Add(iH, jH, sH);
            }
        } else {
            int iH = elem.first - 1;
            int jH = elem.first - 1;
            HS[iH] += 1.0;
            DensityFilteringMatrix->Add(iH, jH, 1.0);
        }
    }
    DensityFilteringMatrix->Finalize();
    DensityFilteringMatrix->covertToMKLCSR();
}


FloatMatrix Sensitivity::DREV(int index, double Eval, FloatArray &EvalD, FloatMatrix &EvecD) {
    int NumberOfDof = EvecD.giveNumberOfRows();
    int NumberOfRepeate = EvecD.giveNumberOfColumns();
    FloatMatrix motaiD(NumberOfDof, NumberOfRepeate);
    FloatMatrix L(NumberOfDof, EvecD.giveNumberOfColumns());
    for (int i = 0; i < EvecD.giveNumberOfColumns(); i++) {
        for (int j = 0; j < EvecD.giveNumberOfColumns(); j++) {
            if (j == i) {
                FloatArray EvecDLoacl(NumberOfDof);
                EvecD.GetColumn(j, EvecDLoacl);
                FloatArray Temp(NumberOfDof);
                ansys_solver->Mass->Mult(EvecDLoacl, Temp);
                L.setColumn(Temp, j + 1);
            } else {
                FloatArray EvecDLoacl(NumberOfDof);
                EvecD.GetColumn(j, EvecDLoacl);
                FloatArray b1(NumberOfDof);
                FloatArray b2(NumberOfDof);
                FloatArray b3(NumberOfDof);
                b1 = 0.0;
                b2 = 0.0;
                b3 = 0.0;
                IntArray edof = ansys_solver->getElementLocation(index - 1);
                FloatMatrix elementStifD = ansys_solver->getElementStif(
                    index - 1, ansys_solver->getStifRhoD(ansys_solver->getRho(index - 1)));
                FloatMatrix elementMassD = ansys_solver->getElementMass(
                    index - 1, ansys_solver->getMassRhoD(ansys_solver->getRho(index - 1)));
                FloatArray elementTemp;
                EvecDLoacl.GetSubVector(edof, elementTemp);
                ansys_solver->Mass->Mult(EvecDLoacl, b2);
                FloatArray eKDP = elementStifD * elementTemp;
                FloatArray eMDP = elementMassD * elementTemp;
                b1.SetSubVector(edof, eMDP);
                b3.SetSubVector(edof, eKDP);
                FloatArray b(NumberOfDof);
                b = b3 - Eval * b1 - EvalD(i) * b2;
                L.setColumn(b, j + 1);
            }
        }
        FloatArray EvecDLoacl(NumberOfDof);
        EvecD.GetColumn(i, EvecDLoacl);
        FloatArray b1(NumberOfDof);
        FloatArray b2(NumberOfDof);
        FloatArray b3(NumberOfDof);
        b1 = 0.0;
        b2 = 0.0;
        b3 = 0.0;
        IntArray edof = ansys_solver->getElementLocation(index - 1);
        FloatMatrix elementStifD2 = ansys_solver->getElementStif(
            index - 1, ansys_solver->getStifRhoD2(ansys_solver->getRho(index - 1)));
        FloatMatrix elementMassD2 = ansys_solver->getElementMass(
            index - 1, ansys_solver->getMassRhoD2(ansys_solver->getRho(index - 1)));
        FloatMatrix elementMassD = ansys_solver->getElementMass(
            index - 1, ansys_solver->getMassRhoD(ansys_solver->getRho(index - 1)));
        FloatArray elementTemp;
        EvecDLoacl.GetSubVector(edof, elementTemp);
        FloatArray eKD2P = elementStifD2 * elementTemp;
        FloatArray eMD2P = elementMassD2 * elementTemp;
        FloatArray eMDP = elementMassD * elementTemp;
        b1.SetSubVector(edof, eMD2P);
        b2.SetSubVector(edof, eMDP);
        b3.SetSubVector(edof, eKD2P);
        b1 = -Eval * b1 - 2.0 * EvalD(i) * b2 + b3;
        FloatArray xishu(EvalD.giveSize());
        for (int j = 0; j < EvecD.giveNumberOfColumns(); j++) {
            EvecD.GetColumn(j, EvecDLoacl);
            if (j == i) {
                xishu(j) = 0.5 * EvecDLoacl * b2;
            } else {
                xishu(j) = 0.5 * EvecDLoacl * b1;
            }
        }
        b3 = 0.0;
        FloatMatrix elementStifD = ansys_solver->getElementStif(
            index - 1, ansys_solver->getStifRhoD(ansys_solver->getRho(index - 1)));
        EvecD.GetColumn(i, EvecDLoacl);
        FloatArray eKDP = elementStifD * elementTemp;
        b3.SetSubVector(edof, eKDP);
        FloatArray bTemp(NumberOfDof);
        ansys_solver->Mass->Mult(EvecDLoacl, bTemp);
        FloatArray b = Eval * b2 + EvalD(i) * bTemp - b3;
        b = b - L * xishu;
        FloatArray DREVx(NumberOfDof);
        DREVx = 0.0;
        FloatArray r = b;
        FloatArray t(NumberOfDof);
        ansys_solver->Mult(r, t);
        double tao0 = t.computeNorm();
        double delt0 = r.computeSquaredNorm();
        FloatArray q = t;
        double v0 = 0;
        double rou0 = r * q;
        FloatArray d(NumberOfDof);
        d = 0.0;
        double delt = r.computeSquaredNorm();
        int ll = 0;
        while (ll < 10 && delt / delt0 > 1.0e-8) {
            b1.beTProductOf(L, q);
            b2 = L * b1;
            b1.resize(NumberOfDof);
            ansys_solver->Mass->Mult(q, b1);
            ansys_solver->Stif->Mult(q, t);
            t = t - Eval * b1 + b2;
            double sigema0 = q * t;
            if (std::abs(sigema0) <= 1.0e-20) {
                break;
            }
            double erfa0 = rou0 / sigema0;
            r = r - erfa0 * t;
            ansys_solver->Mult(r, t);
            double v1 = t.computeNorm();
            v1 = v1 / tao0;
            double c = 1.0 / std::sqrt(1.0 + v1 * v1);
            tao0 = tao0 * v1 * c;
            d = c * c * v0 * v0 * d + c * c * erfa0 * q;
            DREVx = DREVx + d;
            if (abs(rou0) <= 1.0e-20) break;
            double rou1 = r * t;
            double beit = rou1 / rou0;
            q = t + beit * q;
            rou0 = rou1;
            v0 = v1;
            delt = r.computeSquaredNorm();
            ll += 1;
        }
        motaiD.setColumn(DREVx, i + 1);
    }
    return motaiD;
}
