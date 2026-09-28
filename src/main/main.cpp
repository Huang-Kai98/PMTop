#include "../FEM/AnsysLoad.h"
#include "../FEM/BilinearForm.h"
#include "Config.h"
#include "Logger.h"

#include <fstream>
#include <iostream>
#include <ctime>
#include <cstdlib>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ComplexMatrix.h"
#include "ComplexScalar.hpp"
#include "ComplexVector.hpp"
#include "pardiso.hpp"

#include "../Algorithm/TOP.h"


static double cpuTimeSeconds() {
    return static_cast<double>(std::clock()) / static_cast<double>(CLOCKS_PER_SEC);
}


static void setSingleThreadMode() {
#ifdef _WIN32
    _putenv_s("OMP_NUM_THREADS", "1");
    _putenv_s("MKL_NUM_THREADS", "1");
    _putenv_s("MKL_DYNAMIC", "FALSE");
    _putenv_s("OMP_DYNAMIC", "FALSE");
    _putenv_s("OPENBLAS_NUM_THREADS", "1");
#else
    setenv("OMP_NUM_THREADS", "1", 1);
    setenv("MKL_NUM_THREADS", "1", 1);
    setenv("MKL_DYNAMIC", "FALSE", 1);
    setenv("OMP_DYNAMIC", "FALSE", 1);
    setenv("OPENBLAS_NUM_THREADS", "1", 1);
#endif

#ifdef _OPENMP
    omp_set_num_threads(1);
    omp_set_dynamic(0);
#endif
}


auto StifRho = [](const double x) -> double {
    double rho;
    if (x < 0.1) {
        rho = x / 100;
    } else {
        rho = x * x * x;
    }
    return rho;
};

auto MassRho = [](double x) -> double {
    return x;
};

auto DampRho = [](double x) -> double {
    return x;
};

auto StifRhoD = [](const double x) -> double {
    double rho;
    if (x < 0.1) {
        rho = 1 / 100;
    } else {
        rho = 3 * x * x;
    }
    return rho;
};

auto MassRhoD = [](double x) -> double {
    return 1;
};

auto DampRhoD = [](double x) -> double {
    return 1;
};

auto StifRhoD2 = [](const double x) -> double {
    double rho;
    if (x < 0.1) {
        rho = 0;
    } else {
        rho = 6 * x;
    }
    return rho;
};

auto MassRhoD2 = [](double x) -> double {
    return 0;
};

auto DampRhoD2 = [](double x) -> double {
    return 0;
};


int main() {
    setSingleThreadMode();

    Logger::getInstance().init(
        "log.txt",
        spdlog::level::trace,
        true,
        true,
        true,
        false
    );

    auto cfg = Config::loadFromFile("Config.json");

    LOG_INFO("Description: {}", cfg.get<std::string>("description"));

    AnsysLoad ansys(cfg);
    ansys.saveToVTK();

    int N = ansys.getNumberOfElements();

    FloatArray x(N);
    x = 1.0;

    AnsysSolver ansys_solver(
        x,
        ansys,
        StifRho,
        MassRho,
        DampRho,
        StifRhoD,
        MassRhoD,
        DampRhoD,
        StifRhoD2,
        MassRhoD2,
        DampRhoD2
    );

    ansys_solver.QudraticEigenProblem(10, "Right");

    SparseMatrixComplex STIF(*ansys_solver.getStif(), "Real");
    SparseMatrixComplex MASS(*ansys_solver.getMass(), "Real");
    SparseMatrixComplex DAMP(*ansys_solver.getDamp(), "Real");

    auto Eigenvectors_measured =
        ansys_solver.AddNoiseToEigenvectors(
            ansys_solver.getComplexEigenvectors()
        );

    int modeIndex = 0;
    IntArray idx1 = {modeIndex};

    auto Eigenvalue = ansys_solver.getComplexEigenvalue(idx1);
    auto Eigenvector = ansys_solver.getComplexEigenvector(idx1);

    ComplexScalar EE = Eigenvalue(modeIndex);
    auto EV = Eigenvector.GetColumn(modeIndex);

    ansys_solver.Normalize(EE, EV, DAMP, MASS);

    ComplexScalar Dlambda(0.0, 0.0);

    ComplexScalar alpha;
    ComplexScalar DphiL;
    ComplexScalar EigenValueSquare;
    ComplexScalar Db3Dphib0;

    ComplexVector V(Eigenvector.Rows());
    ComplexVector Dphi(Eigenvector.Rows());

    std::cout << ansys_solver.MIR(EV) << std::endl;
    std::cout << ansys_solver.MPC(EV) << std::endl;

    auto Eigenvector_measured =
        Eigenvectors_measured.GetColumn(modeIndex);

    std::cout << ansys_solver.MAC(Eigenvector_measured, EV) << std::endl;

    auto Dsymvector = ansys_solver.DMAC(Eigenvector_measured, EV);
    //auto Dsymvector = ansys_solver.DMIR(EV);
    //auto Dsymvector = ansys_solver.DMPC(EV);


    const int Iterations = 1;
    const int elementSize = ansys_solver.getElementSize();

    std::cout << "Number of Dof: "
              << ansys_solver.getSize()
              << std::endl;

    std::cout << "Number of Elements: "
              << elementSize
              << std::endl;

    // =====================================================================
    // 0. 结果数组
    // =====================================================================

    FloatArray DC_NelsonSingle(elementSize);       // 非伴随标准 Nelson
    FloatArray DC_NelsonAdd(elementSize);          // 非伴随加边法

    FloatArray DC_AdjointAdd(elementSize);         // 伴随加边法
    FloatArray DC_AdjointNelson(elementSize);      // 伴随标准 Nelson
    FloatArray DC_Symm(elementSize);               // 对称/改进方法

    FloatArray Diff(elementSize);

    double time_NelsonSingle_factor = 0.0;
    double time_NelsonSingle_solve = 0.0;
    double time_NelsonSingle_post = 0.0;

    SparseMatrixComplex A1_NelsonSingle;
    PardisoSolver PA_NelsonSingle;
    IntArray g_NelsonSingle;

    ComplexVector Cphi(Eigenvector.Rows());
    ComplexVector Mphi(Eigenvector.Rows());
    ComplexVector RhsCoeff(Eigenvector.Rows());

    double t0 = cpuTimeSeconds();

    ansys_solver.NelsonSinglePrepareFactorized(
        EE,
        EV,
        STIF,
        MASS,
        DAMP,
        A1_NelsonSingle,
        g_NelsonSingle,
        PA_NelsonSingle,
        Cphi,
        Mphi,
        RhsCoeff
    );

    time_NelsonSingle_factor += cpuTimeSeconds() - t0;

    for (int i = 0; i < elementSize; i++) {
        ComplexScalar DeriviateEigenValue(0.0, 0.0);
        ComplexVector DeriviateEigenVector(Eigenvector.Rows());

        t0 = cpuTimeSeconds();

        ansys_solver.NelsonSingleBackSubFactorized(
            EE,
            EV,
            STIF,
            MASS,
            DAMP,
            i,
            PA_NelsonSingle,
            g_NelsonSingle,
            Cphi,
            Mphi,
            RhsCoeff,
            DeriviateEigenValue,
            DeriviateEigenVector
        );

        time_NelsonSingle_solve += cpuTimeSeconds() - t0;

        t0 = cpuTimeSeconds();

        DC_NelsonSingle[i] =
            2.0 * (Dsymvector * DeriviateEigenVector).real();

        time_NelsonSingle_post += cpuTimeSeconds() - t0;
    }

    std::cout << "NelsonSingle factorized CPU time: "
              << " factor = " << time_NelsonSingle_factor
              << " solve = " << time_NelsonSingle_solve
              << " post = " << time_NelsonSingle_post
              << " total = "
              << time_NelsonSingle_factor
               + time_NelsonSingle_solve
               + time_NelsonSingle_post
              << std::endl;

    // =====================================================================
    // 2. 非伴随加边法：NelsonAdd，预分解 + 回代版本
    // =====================================================================

    double time_NelsonAdd_factor = 0.0;
    double time_NelsonAdd_backsub = 0.0;
    double time_NelsonAdd_post = 0.0;

    SparseMatrixComplex Aadd_NelsonAdd;
    PardisoSolver PA_NelsonAdd;

    ComplexVector Mphi_Add(Eigenvector.Rows());
    ComplexVector RhsCoeff_Add(Eigenvector.Rows());

    t0 = cpuTimeSeconds();

    ansys_solver.NelsonAddPrepareFactorized(
        EE,
        EV,
        STIF,
        MASS,
        DAMP,
        Aadd_NelsonAdd,
        PA_NelsonAdd,
        Mphi_Add,
        RhsCoeff_Add
    );

    time_NelsonAdd_factor += cpuTimeSeconds() - t0;

    for (int i = 0; i < elementSize; i++) {
        ComplexScalar DeriviateEigenValue(0.0, 0.0);
        ComplexVector DeriviateEigenVector(Eigenvector.Rows());

        t0 = cpuTimeSeconds();

        ansys_solver.NelsonAddBackSubFactorized(
            EE,
            EV,
            STIF,
            MASS,
            DAMP,
            i,
            PA_NelsonAdd,
            Mphi_Add,
            RhsCoeff_Add,
            DeriviateEigenValue,
            DeriviateEigenVector
        );

        time_NelsonAdd_backsub += cpuTimeSeconds() - t0;

        t0 = cpuTimeSeconds();

        DC_NelsonAdd[i] =
            2.0 * (Dsymvector * DeriviateEigenVector).real();

        time_NelsonAdd_post += cpuTimeSeconds() - t0;
    }

    std::cout << "NelsonAdd factorized non-adjoint average CPU time: "
              << " factor: " << time_NelsonAdd_factor / Iterations
              << " backsub: " << time_NelsonAdd_backsub / Iterations
              << " post: " << time_NelsonAdd_post / Iterations
              << " total: "
              << (time_NelsonAdd_factor
                + time_NelsonAdd_backsub
                + time_NelsonAdd_post) / Iterations
              << std::endl;

    // =====================================================================
    // 3. 伴随加边法：AdjointNelosnAdd
    // =====================================================================

    double time1 = 0.0;
    double time2 = 0.0;

    for (int iter = 0; iter < Iterations; iter++) {
        double t0 = cpuTimeSeconds();

        ansys_solver.AdjointNelosnAdd(
            EE,
            EV,
            STIF,
            MASS,
            DAMP,
            Dlambda,
            Dsymvector,
            alpha,
            V
        );

        time1 += cpuTimeSeconds() - t0;

        t0 = cpuTimeSeconds();

        for (int i = 0; i < elementSize; i++) {
            DC_AdjointAdd[i] =
                ansys_solver.AdjointElement(
                    EE,
                    EV,
                    alpha,
                    V,
                    i
                ).real();
        }

        time2 += cpuTimeSeconds() - t0;
    }

    std::cout << "AdjointNelosnAdd average CPU time: "
              << " t1: " << time1 / Iterations
              << " t2: " << time2 / Iterations
              << " total: " << (time1 + time2) / Iterations
              << std::endl;

    // =====================================================================
    // 4. 伴随标准 Nelson：AdjointNelson
    // =====================================================================

    time1 = 0.0;
    time2 = 0.0;

    for (int iter = 0; iter < Iterations; iter++) {
        double t0 = cpuTimeSeconds();

        ansys_solver.AdjointNelson(
            EE,
            EV,
            STIF,
            MASS,
            DAMP,
            Dlambda,
            Dsymvector,
            alpha,
            V
        );

        time1 += cpuTimeSeconds() - t0;

        t0 = cpuTimeSeconds();

        for (int i = 0; i < elementSize; i++) {
            DC_AdjointNelson[i] =
                ansys_solver.AdjointElement(
                    EE,
                    EV,
                    alpha,
                    V,
                    i
                ).real();
        }

        time2 += cpuTimeSeconds() - t0;
    }

    std::cout << "AdjointNelson average CPU time: "
              << " t1: " << time1 / Iterations
              << " t2: " << time2 / Iterations
              << " total: " << (time1 + time2) / Iterations
              << std::endl;

    Diff = DC_AdjointAdd - DC_AdjointNelson;

    std::cout << "Diff AdjointAdd vs AdjointNelson Norm2: "
              << Diff.computeNorm()
              << std::endl;

    // =====================================================================
    // 5. 对称/改进方法：Adjointsymm
    // =====================================================================

    PardisoSolver PA;
    PA.SetMatrixType(PardisoSolver::COMPLEX_NONSYMMETRIC);
    PA.SetOperator(STIF);

    time1 = 0.0;
    time2 = 0.0;

    for (int iter = 0; iter < Iterations; iter++) {
        double t0 = cpuTimeSeconds();

        ansys_solver.AdjointsymmDsymvector(
            EE,
            EV,
            STIF,
            MASS,
            DAMP,
            PA,
            Dsymvector,
            Dphi,
            DphiL,
            EigenValueSquare,
            Db3Dphib0,
            1e-5,
            5
        );

        time1 += cpuTimeSeconds() - t0;

        t0 = cpuTimeSeconds();

        for (int i = 0; i < elementSize; i++) {
            DC_Symm[i] =
                ansys_solver.AdjointsymmElement(
                    EE,
                    EV,
                    Dphi,
                    DphiL,
                    EigenValueSquare,
                    Db3Dphib0,
                    Dlambda,
                    i
                ).real();
        }

        time2 += cpuTimeSeconds() - t0;
    }

    Diff = DC_Symm - DC_AdjointNelson;

    std::cout << "Adjointsymm average CPU time: "
              << " t1: " << time1 / Iterations
              << " t2: " << time2 / Iterations
              << " total: " << (time1 + time2) / Iterations
              << " Norm2_vs_AdjointNelson: " << Diff.computeNorm()
              << std::endl;

    // =====================================================================
    // 6. 代表性单元结果输出
    // =====================================================================

    int maxId = DC_AdjointNelson.Abs().giveIndexMaxElem() - 1;

    std::cout << "Representative element index: "
              << maxId
              << std::endl;

    std::cout << "DC_NelsonSingle: "
              << DC_NelsonSingle[maxId]
              << " DC_NelsonAdd: "
              << DC_NelsonAdd[maxId]
              << " DC_AdjointAdd: "
              << DC_AdjointAdd[maxId]
              << " DC_AdjointNelson: "
              << DC_AdjointNelson[maxId]
              << " DC_Symm: "
              << DC_Symm[maxId]
              << std::endl;

    // =====================================================================
    // 7. 可选保存结果
    // =====================================================================

    // ansys.saveElementValue("NELSON_SINGLE", DC_NelsonSingle);
    // ansys.saveElementValue("NELSON_ADD", DC_NelsonAdd);
    // ansys.saveElementValue("ADJOINT_ADD", DC_AdjointAdd);
    // ansys.saveElementValue("ADJOINT_NELSON", DC_AdjointNelson);
    // ansys.saveElementValue("ADJOINT_SYMM", DC_Symm);

    return 0;
}