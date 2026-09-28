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
int main() {

    std::string datafile = "BaiCheShen";
    AnsysLoad anasys(datafile);

    BilinearForm a(&anasys);

    FloatArray x(static_cast<int>(anasys.getNumberOfDof()));
    x=1.0;
    a.assemble(x);
    std::unordered_set<int> ConstrainedDofs;
    anasys.getConstrainedDofs(ConstrainedDofs);
    a.EliminateEssentialBC(ConstrainedDofs,1);
    a.Finalize();
    a.convertToMKL();

    const SparseMatrix &Stif = a.SpStif();
    const SparseMatrix &Mass = a.SpMass();

    PardisoSolver pardiso_solver;
    pardiso_solver.SetMatrixType(PardisoSolver::REAL_NONSYMMETRIC);
    pardiso_solver.SetOperator(Stif);

    ArpackSolver M(Stif.Size(),10);
    M.SetOption(arpack::which::largest_magnitude);
    M.Mult(pardiso_solver,Mass);

    FloatArray Eval=M.getEigenvalues();
    FloatArray EVec=M.getEigenvector(1);
    Eval.printYourself();

    anasys.saveToVTK();
    anasys.savePointValue(EVec);
    return 0;
}