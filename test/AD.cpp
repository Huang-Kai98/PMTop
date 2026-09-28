// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.


#include <pardiso.hpp>
#include <tic_toc.hpp>
#include "Arpack.hpp"
#include "Optimizer.h"
#include "BilinearForm.h"
#include "../Algorithm/TOP.h"
#include "sqmr.cpp"

int main() {
    std::string datafile = "aero";
    AnsysLoad anasys(datafile);

    FloatArray x(static_cast<int>(anasys.getNumberOfDof()));
    x=1.0;

    auto StifRho = [](const double x) -> double {
        return x*x*x;
    };

    auto MassRho = [](double x)->double {
        return x;
    };

    auto StifRhoD = [](const double x) -> double {
        return 3*x*x;
    };

    auto MassRhoD = [](double x)->double {
        return 1;
    };

    auto ansys_solver=std::make_shared<AnsysSolver>(x,anasys,StifRho,MassRho,
        StifRhoD,MassRhoD);

    ansys_solver->EigenProblem(10);

    int N=ansys_solver->getSize();
    double Eval=ansys_solver->getEigenValues()[0];
    FloatArray Evec=ansys_solver->getEigenVector(0);
    auto obj = [&N]() -> FloatArray {
        FloatArray x(N);
        x=0.0;
        return x;
    };

    auto EigenValueDerivative=[&Eval,&Evec]() -> FloatArray {
        FloatArray x(1);
        x=-Evec.dotProduct(Evec)/(Eval*Eval);
        return x;
    };

    auto EigenVectorDerivative=[&Eval,&Evec]() -> FloatMatrix {
        FloatArray tx=(2.0/Eval)*Evec;
        FloatMatrix x(tx.giveSize(),1);
        for (int i=0;i<tx.giveSize();i++) {
            x(i,0)=tx(i);
        }
        return x;
    };

    auto objDerivative=std::make_shared<FunctionDerivative>(obj,EigenValueDerivative,EigenVectorDerivative);

    Sensitivity sen(ansys_solver,objDerivative);
    FloatArray DC=sen.AD(0,0);


    anasys.saveToVTK();
    anasys.saveElementValue(DC);

    return 0;
}
