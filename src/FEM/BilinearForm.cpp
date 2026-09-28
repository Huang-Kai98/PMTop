// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

//
// Created by huangkai on 25-1-21.
//

#include "BilinearForm.h"

void BilinearForm::assembleStif(FloatArray x, int skip_zeros){
    stif = new SparseMatrix(fes->getNumberOfDof(),fes->getNumberOfDof());
    for (int i=0;i<fes->getNumberOfElements();i++) {
        stif->AddSubMatrix(fes->getElementLocation(i),
            fes->getElementLocation(i),fes->getElementStif(i,x[i]),skip_zeros);
    }
    LOG_TRACE("Stifness Assemble Finish.");
}
void BilinearForm::assembleStif(FloatArray x, std::function<double(double x)> StifRho, int skip_zeros) {
    stif = new SparseMatrix(fes->getNumberOfDof(),fes->getNumberOfDof());
    for (int i=0;i<fes->getNumberOfElements();i++) {
        stif->AddSubMatrix(fes->getElementLocation(i),
            fes->getElementLocation(i),fes->getElementStif(i,StifRho(x[i])),skip_zeros);
    }
    LOG_TRACE("Stifness Assemble Finish.");
}

void BilinearForm::assembleMass(FloatArray x, int skip_zeros){
    mass = new SparseMatrix(fes->getNumberOfDof(),fes->getNumberOfDof());
    for (int i=0;i<fes->getNumberOfElements();i++) {
        mass->AddSubMatrix(fes->getElementLocation(i),
            fes->getElementLocation(i),fes->getElementMass(i,x[i]),skip_zeros);
    }
    LOG_TRACE("Mass Assemble Finish.");
}

void BilinearForm::assembleMass(FloatArray x, std::function<double(double x)> MassRho, int skip_zeros) {
    mass = new SparseMatrix(fes->getNumberOfDof(),fes->getNumberOfDof());
    for (int i=0;i<fes->getNumberOfElements();i++) {
        std::string type = fes->getElementType(i); // 获取一次
        if (type == "Combin14" || type == "COMBIN14") {
            continue;
        }
        mass->AddSubMatrix(fes->getElementLocation(i),
            fes->getElementLocation(i),fes->getElementMass(i,MassRho(x[i])),skip_zeros);
    }
    LOG_TRACE("Mass Assemble Finish.");
}

void BilinearForm::assembleDamp(FloatArray x, std::function<double(double x)> DampRho, int skip_zeros) {
    damp = new SparseMatrix(fes->getNumberOfDof(),fes->getNumberOfDof());
    for (int i=0;i<fes->getNumberOfElements();i++) {
        damp->AddSubMatrix(fes->getElementLocation(i),
            fes->getElementLocation(i), fes->getElementDamp(i,DampRho(x[i])),skip_zeros);
    }
    LOG_TRACE("Damp Assemble Finish.");
}

void BilinearForm::assemble(FloatArray x, std::function<double(double x)> StifRho, std::function<double(double x)> MassRho, std::function<double(double x)> DampRho, int skip_zeros) {
    assembleStif(x,StifRho,skip_zeros);
    assembleMass(x,MassRho,skip_zeros);
    assembleDamp(x,DampRho,skip_zeros);
}



void BilinearForm::assemble(const FloatArray& x, int skip_zeros) {
    assembleStif(x,skip_zeros);
    assembleMass(x,skip_zeros);
}

void BilinearForm::assemble(FloatArray x, std::function<double(double x)> StifRho, std::function<double(double x)> MassRho, int skip_zeros) {
    assembleStif(x,StifRho,skip_zeros);
    assembleMass(x,MassRho,skip_zeros);
}


void BilinearForm::EliminateEssentialBC(std::unordered_set<int> ConstrainedDofs, int d) {

    for (auto & constrainedDof : ConstrainedDofs) {
        if (constrainedDof>=0) {
            stif->EliminateRowCol(constrainedDof,d);
        }
        else {
            stif->EliminateRowCol(-1-constrainedDof,d);
        }
    }
    LOG_TRACE("EliminateEssentialBC Finish.");
}
void BilinearForm::Finalize(int skip_zeros) {
    stif -> Finalize (skip_zeros);
    if (mass!=nullptr)
        mass -> Finalize (skip_zeros);
    if (damp!=nullptr)
        damp -> Finalize (skip_zeros);
    LOG_TRACE("BilinearForm::Finalize Finish.");
}


void BilinearForm::convertToMKL() {
    if (!stif->isUseMKL())
        stif->covertToMKLCSR();
    if (!mass->isUseMKL())
        mass->covertToMKLCSR();
    LOG_TRACE("ConverToMkl Finish.");
}




