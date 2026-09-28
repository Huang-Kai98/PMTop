// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

//
// Created by huangkai on 25-2-15.
//

#ifndef MAIN_H
#define MAIN_H

#include <chrono>
#include <pardiso.hpp>
#include <set>
#include <tic_toc.hpp>
#include "Arpack.hpp"
#include "Optimizer.h"
#include "BilinearForm.h"
#include "../Algorithm/TOP.h"
#include "sqmr.cpp"



auto StifRho = [](const double x) -> double {
    double rho;
    if (x<0.1) {
        rho=x/100;
    }
    else {
        rho = x*x*x;
    }
    return rho;
};

auto MassRho = [](double x)->double {
    return x;
};

auto DampRho = [](double x)->double {
    return x;
};

auto StifRhoD = [](const double x) -> double {
    double rho;
    if (x<0.1) {
        rho=1/100;
    }
    else {
        rho = 3*x*x;
    }
    return rho;
};

auto MassRhoD = [](double x)->double {
    return 1;
};

auto StifRhoD2 = [](const double x) -> double {
    double rho;
    if (x<0.1) {
        rho=0;
    }
    else {
        rho = 6*x;
    }
    return rho;
};

auto MassRhoD2 = [](double x)->double {
    return 0;
};

void ConstraintPoint(FloatArray ConstrainedNodelIndex, int dim, FloatArray& XX) {
    for (int i = 0; i < ConstrainedNodelIndex.giveSize(); i++) {

        for (int j = 0; j < dim; j++) {
            XX(j*dim+j)=1;
        }
    }
}
class TopProblem {
private:
    int NumberOfElement;
    int NumberOfDof;
    double VolumeRatio;
    int NumberOfConstrains;
    FloatArray OldEigenValue;
    FloatMatrix OldEigenVector;
    IntArray OptimizedModalIndex;
    IntArray ConstrainedNodelIndex;
    double BigNumber;
    IntArray Mask;
    std::vector<std::pair<double, std::vector<int>>> Repeat;
    IntArray OMindex;
    IntArray uniqueIndices;
    std::vector<IntArray> repeatedGroups;
    std::shared_ptr<AnsysSolver> ansys_solver;
    std::shared_ptr<AnsysLoad> ansys_load;
    std::unique_ptr<SparseMatrix> DensityFilteringMatrix;
    FloatArray HS;
    FloatArray XX;
    double FrequencyConstrain;




public:
    int n, m;
    FloatArray x, xold, xnew;
    FloatArray df, g, gnew, dg;
    FloatArray xmin, xmax;
    std::vector<FloatArray> EvalIter;
    FloatArray ObjValue;
    std::vector<FloatArray> ConstrainValue;
    KDNode *KDNodeTree;
    TopProblem(FloatArray x0,int NumberOfConstrains_,double VolumeRatio_,std::shared_ptr<AnsysLoad> ansys_load_
        ,FloatArray xmin0,FloatArray xmax0,IntArray& OptimizedModalIndex_,
        double FrequencyConstrain_=10,double BigNmber_=10000000000)
        :ansys_load(ansys_load_)
        ,OptimizedModalIndex(OptimizedModalIndex_)
        ,FrequencyConstrain (FrequencyConstrain_)
        ,NumberOfConstrains(NumberOfConstrains_)
        ,VolumeRatio(VolumeRatio_)
        ,BigNumber(BigNmber_)
        ,x(x0)
        ,xmin(xmin0)
        ,xmax(xmax0)
    {
        NumberOfElement=ansys_load->getNumberOfElements();
        NumberOfDof=ansys_load->getNumberOfDof();
        m=NumberOfConstrains;
        n=ansys_load->NumberOfElement;
        xold=x;
        xnew.resize(n);
        df.resize(n);
        g.resize(m);
        gnew.resize(m);
        dg.resize(n * m);

    }
    void objFunc(FloatArray& x, double& f0x, FloatArray& fx);
    void sensFunc(FloatArray &x, double &f0x,
                  FloatArray &fx, FloatArray &df0dx, FloatArray &dfdx);

    void ConstraintPoint(std::vector<Node>& ConstrainNode,int Dim);
    void ConstraintPoint(FloatArray& XX,int Dim);

    FloatArray LambdaSens(int index, std::shared_ptr<AnsysSolver> &ansys);
    std::map<int,std::vector<SubEigenProblem>> DifferentiableEigenvectorCalculation(double Eval, FloatMatrix Evec,
        std::shared_ptr<AnsysSolver>& ansys_solver);
    IntArray GetOptimizedModalIndex() {
        IntArray optimizedModalIndex;
        FloatArray temp=ansys_solver->getEigenValues();
        for (int i=0;i<temp.giveSize();i++) {
            if (temp(i)>100&&temp(i)<48000) {
                optimizedModalIndex.append(i);
            }
        }
        std::cout<<optimizedModalIndex<<std::endl;
        /*
        IntArray optimizedModalIndex;
        optimizedModalIndex.resize(OMindex.giveSize());
        if (OldEigenVector.isNotEmpty()) {
            for (int i = 0; i < OptimizedModalIndex.giveSize(); i++) {
                int index = OptimizedModalIndex(i);
                FloatArray temp=ansys_solver->getEigenVector(index);
                FloatArray MAC;
                MAC.beProductOf(OldEigenVector, temp);
                FloatArray OldEigenVectorNorm(OldEigenVector.giveNumberOfRows());
                for (int j=0;j<OldEigenVector.giveNumberOfRows();j++) {
                    FloatArray tempOld(OldEigenVector.giveNumberOfColumns());
                    for (int k=0;k<OldEigenVector.giveNumberOfColumns();k++) {
                        tempOld[k]=OldEigenVector(j,k);
                    }
                    OldEigenVectorNorm[j]=tempOld.computeSquaredNorm();
                }
                double tempNorm=temp.computeSquaredNorm();
                for (int j=0;j<OldEigenVector.giveNumberOfRows();j++) {
                    MAC[j]=(MAC[j]*MAC[j])/(OldEigenVectorNorm[j]*tempNorm);
                }

                int newindex= MAC.giveIndexMaxElem()-1;
                optimizedModalIndex(i)=newindex;
            }

        }
        else {
            optimizedModalIndex=OptimizedModalIndex;
        }
        OptimizedModalIndex=optimizedModalIndex;
        std::cout<<optimizedModalIndex<<std::endl;
        */
        return optimizedModalIndex;
    }
    std::pair<IntArray,std::vector<IntArray>>
    ClassifiedEigenValue(const std::vector<std::pair<double, std::vector<int>>>& eigenGroups,
                                 const IntArray& sortedIndices);
    double ComputeConstrains(const FloatArray& x) {
        double G=0.0;
        double V=0.0;
        for (int i=0;i<Mask.giveSize();i++) {
            G+=Mask[i]*x[i];
            if (Mask[i]==1) {
                V+=1.0;
            }
        }
        return G - VolumeRatio * V;
    }
    void TopSet() {
        Mask.resize(NumberOfElement);

        //for (int i=1;i<=NumberOfElement;i++) {
        //    if (ansys_load->elements[i].sec==1) {
        //        Mask[i-1]=1;
        //    }
        //}

        //for (int i=1;i<=NumberOfElement;i++) {
        //    if (ansys_load->elements[i].sec==2) {
        //        Mask[i-1]=1;
        //    }
        //    if (ansys_load->elements[i].sec==3) {
        //        Mask[i-1]=0;
        //    }
        //}
        for (int i=1600;i<=3200;i++) {Mask[i-1]=1;}
        for (int i=1;i<=1600;i++) {Mask[i-1]=0;}


    }

    void SetTopElement() {
        for (int i=1;i<=Mask.giveSize();i++) {
            if (Mask[i-1]==1) {
                ansys_load->AnasysElements[i]->setTop(true);
            }
            else {
                ansys_load->AnasysElements[i]->setTop(false);
            }
        }
    }
    double ComputeChange(const FloatArray &x, FloatArray &xold);
    void SolveMMA(int maxIterations, double tolerance);
    void SolveGCMMA(int maxIterations, double tolerance);

    void CreateDensityFilteringMatrix(double radius=std::sqrt(3));

    void OC();
    static void PrintVector(const FloatArray &vec, const std::string &name) {
        std::cout << name << ":";
        for (const double v : vec) {
            std::cout << " " << v;
        }
        std::cout << std::endl;
    }

};

void TopProblem::ConstraintPoint(std::vector<Node>& ConstrainNode,int Dim) {
    XX.resize(NumberOfDof);
    XX=0.0;

    KDNodeTree = ansys_load->buildKDTree(ansys_load->getAnasysNodes(),3,0);
    for (int i=0;i<ConstrainNode.size();i++) {
        auto RNode= ansys_load->rangeSearch(KDNodeTree,ConstrainNode[i],0.001,0,3);
        for (int j=0;j<Dim/2;j++) {
            for (int k=0;k<RNode.size();k++) {
                XX[(RNode[k].id-1)*Dim+j] =1;
            }
        }
    }

    //IntArray indices={12849,12556,9470,11409,11264,11140,11695,
    //                  12147,11769,9379,10454,9919,9515,9100};
//
 //   for (int j=0;j<Dim/2;j++) {
 //       for (int k=0;k<indices.giveSize();k++) {
 //           XX[(indices[k]-1)*Dim+j] =1;
  //      }
  //  }
}

inline void TopProblem::ConstraintPoint(FloatArray &CN, int Dim) {
    XX.resize(NumberOfDof);
    XX=0.0;


    for (int i=0;i<CN.giveSize();i++) {
        for (int j=0;j<Dim;j++) {
            XX[CN[i]*Dim+j] =1;
        }
    }

}



void TopProblem::CreateDensityFilteringMatrix(double radius) {
    DensityFilteringMatrix = std::make_unique<SparseMatrix>(ansys_load->getNumberOfElements(),
      ansys_load->getNumberOfElements());
    HS.resize(ansys_load->getNumberOfElements());
    auto Centroid=ansys_load->getCentroid();
    auto KDTree=ansys_load->buildKDTree(Centroid,3,0);
    for (const auto& elem:ansys_load->AnasysElements) {
        if (elem.second->getTopStatus()==true) {
            auto RNode= ansys_load->rangeSearch(KDTree,Centroid[elem.first-1],radius,0,3);
            //auto RNode=ansys_solver->rangeSearch(KDTree,Centroid[elem.first-1],radius);
            int iH=elem.first-1;
            for (const auto& rnode:RNode) {
                int jH=rnode.id-1;
                double dist=Centroid[elem.first-1].euclideanDistance(rnode);
                double sH=std::max(0.0,radius-dist);
                HS[iH]+=sH;
                DensityFilteringMatrix->Add(iH,jH,sH);
            }
        }
        else {
            int iH=elem.first-1;
            int jH=elem.first-1;
            HS[iH]+=1.0;
            DensityFilteringMatrix->Add(iH,jH,1.0);
        }
    }
    DensityFilteringMatrix->Finalize();
    DensityFilteringMatrix->covertToMKLCSR();
}

void TopProblem::OC() {
    double f0x;
    FloatArray fx;
    FloatArray f0dx;
    FloatArray dfdx;
    sensFunc(x,f0x,fx,f0dx,dfdx);
    int lowIndex=f0dx.giveIndexMinElem();
    int highIndex=f0dx.giveIndexMaxElem();
    double low=f0dx.at(lowIndex);
    double high=f0dx.at(highIndex);
    double dx_min=0.01;
    double dx_max=0.1;
    double mid=0.0;
    FloatArray dx(NumberOfElement);
    dx=0.5;
    dx=x-dx;
    dx=dx.Abs();
    FloatArray tmp(NumberOfElement);
    tmp=1.0;
    dx=tmp-2*dx;
    dx=(dx_max-dx_min)*dx;
    tmp=dx_min;
    dx+=tmp;

    FloatArray temp01(NumberOfElement);
    temp01=0.001;
    FloatArray temp1(NumberOfElement);
    temp1=1.0;
    FloatArray new_x(NumberOfElement);

    while(((high-low)/(high))>0.00001) {
        tmp=mid;
        tmp-=f0dx;
        tmp=tmp.Sign();
        tmp=tmp.HadamardProduct(dx);
        tmp+=x;
        tmp.beMaxOf(temp01,tmp);
        new_x.beMinOf(temp1,tmp);
        if ((new_x.sum()-VolumeRatio*NumberOfElement)>0) {
            high=mid;
            mid=(low+high)/2;
        }
        else {
            break;
        }
    }
    x=new_x;
    std::cout<<f0x<<" "<<x.sum()/NumberOfElement<<std::endl;
}

void TopProblem::SolveMMA(int maxIterations, double tolerance) {
    double f, fnew;
    double ch = 1.0;

    objFunc(x, f, g);
    //PrintVector(x, "Initial x");
    std::cout << "Initial objective: " << f << std::endl;
    PrintVector(g, "Initial constraints");

    MMASolver mma(n,m);
    for (int iter = 0; ch > tolerance && iter < maxIterations; ++iter) {
        sensFunc(x ,f, g, df,dg);
        mma.Update(x.givePointer(), df.givePointer(), g.givePointer(),
            dg.givePointer(),xmin.givePointer(), xmax.givePointer());
        ch = ComputeChange(x, xold);
        printf("Iteration: %d, Objective: %f, Change: %f\n", iter, f, ch);
        objFunc(x, f, g);
        PrintVector(g, "Constraints");
        std::cout << std::endl;
    }
}

inline void TopProblem::SolveGCMMA(int maxIterations, double tolerance) {
    double f, fnew;
    double ch = 1.0;

    objFunc(x, f, g);

    ansys_load->saveToVTK("Original");
    FloatMatrix EigenVector=ansys_solver->getEigenvector(OMindex);
    ansys_load->savePointValueShell(OMindex,EigenVector);

    std::cout << "Initial objective: " << f << std::endl;
    PrintVector(g, "Initial constraints");
    ConstrainValue.push_back(g);

    GCMMASolver gcmma(n,m);
    for (int iter = 0; ch > tolerance && iter < maxIterations; ++iter) {
        EvalIter.push_back(ansys_solver->getEigenValues());
        sensFunc(x ,f, g, df,dg);
        gcmma.OuterUpdate(xnew.givePointer(), x.givePointer(), f, df.givePointer(),
            g.givePointer(), dg.givePointer(), xmin.givePointer(), xmax.givePointer());
        objFunc(xnew, fnew, gnew);
        bool conserv = gcmma.ConCheck(fnew, gnew.givePointer());
        for (int innerIter = 0; !conserv && innerIter < 15; ++innerIter) {
            gcmma.InnerUpdate(xnew.givePointer(), fnew, gnew.givePointer(), x.givePointer(),
                f, df.givePointer(), g.givePointer(), dg.givePointer(),
                xmin.givePointer(), xmax.givePointer());
            objFunc(xnew, fnew, gnew);
            conserv = gcmma.ConCheck(fnew, gnew.givePointer());
        }
        x = xnew;
        for (int iii=0;iii<1600;iii++) {
            x[iii]=1.0;
        }
        ansys_load->saveToVTK("iter"+std::to_string(iter));
        ansys_load->saveElementValue(x);
        ch = ComputeChange(x, xold);
        printf("Iteration: %d, Objective: %f, Change: %f\n", iter, f, ch);
        ObjValue.push_back(f);
        objFunc(x, f, g);
        PrintVector(g, "Constraints");
        ConstrainValue.push_back(g);
        std::cout << std::endl;
    }

    /*
    for (int index=0;index<x.giveSize();index++) {
        if (x[index]<0.1) {
            auto lowElement=ansys_load->elements[index];
            auto lowNode=lowElement.nodes;
            for (const auto&node:lowNode) {
                for (int j=0;j<6;j++) {
                    ansys_solver->SetEigenvector(node*6+j,0.0);
                }
            }
        }
    }

    ansys_solver->Reorthogonalization();
*/
    ansys_load->saveToVTK("Opt");
    ansys_load->saveElementValue(x);
    EigenVector=ansys_solver->getEigenvector(OMindex);
    ansys_load->savePointValueShell(OMindex,EigenVector);


}


void TopProblem::objFunc(FloatArray& x,double& f0x, FloatArray& fx) {
    ansys_solver = std::make_shared<AnsysSolver>(x, *ansys_load, StifRho, MassRho,
                                                 StifRhoD, MassRhoD, StifRhoD2, MassRhoD2);
    const int nev = 15;
    ansys_solver->EigenProblem(nev);
    FloatArray Eval = ansys_solver->getEigenValues();

    NumberOfDof = ansys_load->getNumberOfDof();
    NumberOfElement = ansys_load->getNumberOfElements();

    Repeat = Eval.findDuplicatesWithIndices();

/*
    for (int index=0;index<x.giveSize();index++) {
        if (x[index]<0.1) {
            ansys_solver->SetEigenvector(index,0.0);
        }
    }
    ansys_solver->Reorthogonalization();

    for (int index=0;index<x.giveSize();index++) {
        if (x[index]<0.1) {
            auto lowElement=ansys_load->elements[index];
            auto lowNode=lowElement.nodes;
            for (const auto&node:lowNode) {
                for (int j=0;j<6;j++) {
                    ansys_solver->SetEigenvector(node*6+j,0.0);
                }
            }
        }
    }
    ansys_solver->Reorthogonalization();
    */

    OMindex = GetOptimizedModalIndex();

    auto [uniqueIndices_,repeatedGroups_] = ClassifiedEigenValue(Repeat, OMindex);
    uniqueIndices = uniqueIndices_;
    repeatedGroups = repeatedGroups_;

     FloatArray TempEigenVector=ansys_solver->getEigenVectors();
    OldEigenValue=ansys_solver->getEigenValues();
    OldEigenVector.resize(nev,NumberOfDof);
    for (int i=0;i<nev;i++) {
        for (int j=0;j<NumberOfDof;j++) {
            OldEigenVector(i,j)=TempEigenVector(i*NumberOfDof+j);
        }
    }



    f0x = 0.0;

    for (int i = 0; i < OMindex.giveSize(); i++) {
        FloatArray TEvec = ansys_solver->getEigenVector(OMindex(i));
        FloatArray XTEvec = TEvec.HadamardProduct(XX);
        f0x += XTEvec * TEvec;
    }
    f0x *= BigNumber;
    fx.resize(NumberOfConstrains);
    fx[0] = ComputeConstrains(x);
    fx[1] = FrequencyConstrain - Eval[0];

}

void TopProblem::sensFunc(FloatArray &x, double &f0x, FloatArray &fx, FloatArray &df0dx, FloatArray &dfdx) {
    objFunc(x, f0x, fx);
    FloatArray f0dx(NumberOfElement);
    f0dx = 0.0;
    dfdx.resize(NumberOfElement * NumberOfConstrains);

    Sensitivity sen(ansys_solver);
    for (int i = 0; i < uniqueIndices.giveSize(); i++) {
        double Eval = ansys_solver->getEigenValues()[uniqueIndices(i)];
        auto EVec = ansys_solver->getEigenVector(uniqueIndices(i));
        FloatArray EigenVectorDerivative = 2 * EVec.HadamardProduct(XX);
        FloatArray EigenValueDerivative(NumberOfElement);
        EigenValueDerivative = 0.0;
        FloatArray ObjDerivative(NumberOfElement);
        ObjDerivative = 0.0;
        f0dx += sen.AdjointMethod(Eval, EVec, ObjDerivative, EigenValueDerivative, EigenVectorDerivative,50,1e-15);
    }

    for (auto RepeatedIndices: repeatedGroups) {
        double Eval = ansys_solver->getEigenValues()[RepeatedIndices[0]];
        FloatMatrix Evec(NumberOfDof, RepeatedIndices.giveSize());
        FloatMatrix motaiD(NumberOfDof, RepeatedIndices.giveSize());
        for (int j = 0; j < RepeatedIndices.giveSize(); j++) {
            auto LocalEvec = ansys_solver->getEigenVector(RepeatedIndices(j));
            Evec.setColumn(LocalEvec, j + 1);
        }
        auto Sub_Eigen_Problem = DifferentiableEigenvectorCalculation(Eval, Evec, ansys_solver);
        FloatArray DC(NumberOfElement);
        DC = 0.0;
        for (const auto &elem: ansys_load->elements) {
            if (ansys_solver->ansys.AnasysElements[elem.first]->getTopStatus()==true) {
                auto SEP = Sub_Eigen_Problem[elem.first - 1][0];
                auto EvalD = SEP.SubEigenValue;
                auto EvecD = Evec * SEP.SubEigenVector;
                motaiD = sen.DREV(elem.first, Eval, EvalD, EvecD);
                FloatMatrix EigenVectorDerivative = EvecD.HadamardProduct(XX);
                EigenVectorDerivative.times(2.0);
                FloatArray EigenValueDerivative(NumberOfElement);
                EigenValueDerivative = 0.0;
                FloatArray ObjDerivative(NumberOfElement);
                ObjDerivative = 0.0;
                for (int i = 0; i < SEP.SubEigenValue.giveSize(); i++) {
                    FloatArray TempMotaiD(NumberOfDof);
                    FloatArray TempEigenVectorDerivative(NumberOfDof);
                    motaiD.GetColumn(i, TempMotaiD);
                    EigenVectorDerivative.GetColumn(i, TempEigenVectorDerivative);
                    DC(elem.first - 1) += TempEigenVectorDerivative * TempMotaiD;
                }
            }
            else {
                DC[elem.first - 1]=0.0;
            }
        }
        f0dx += DC;
    }

    FloatArray LambdaD=LambdaSens(0,ansys_solver);


    for (int i = 0; i < NumberOfElement; i++) {
        for (int j = 0; j < NumberOfConstrains; j++) {
            dfdx[i*NumberOfConstrains]=1;
            dfdx[i*NumberOfConstrains + 1]=LambdaD[i];
        }
    }

    df0dx.resize(NumberOfElement);
    for (int i = 0; i < NumberOfElement; i++) {
        df0dx[i]=f0dx[i]/HS[i];
    }
    FloatArray dft(NumberOfElement);
    DensityFilteringMatrix->Mult(df0dx,dft);

    df0dx=dft*BigNumber;

}

double TopProblem::ComputeChange(const FloatArray &x, FloatArray &xold) {
    double ch = 0.0;
    for (int i = 0; i < x.giveSize(); ++i) {
        ch = std::max(ch, std::abs(x[i] - xold[i]));
        xold[i] = x[i];
    }
    return ch;
}

FloatArray TopProblem::LambdaSens(int index, std::shared_ptr<AnsysSolver>& ansys) {
    FloatArray Result(NumberOfElement);
    Result=0.0;
    int uniqueIndex=uniqueIndices.findSorted(index);
    int repeateGroupindex=0;
    int repeateLocalindex=0;
    for (int i=0;i<repeatedGroups.size();i++) {
        repeateLocalindex=repeatedGroups[i].findSorted(index);
        if (repeateLocalindex!=0) {
            repeateGroupindex=i+1;
        }
    }
    if (uniqueIndex!=0) {
        double Eval = ansys_solver->getEigenValues()[uniqueIndex-1];
        FloatArray Evec=ansys_solver->getEigenVector(uniqueIndex-1);
        for (int index = 0;index<ansys_solver->getElementSize();index++) {
            if (ansys_solver->ansys.AnasysElements[index+1]->getTopStatus()==true) {
                IntArray edof=ansys_solver->getElementLocation(index);
                FloatMatrix elementStifD=ansys_solver->getElementStif(index,ansys_solver->getStifRhoD(ansys_solver->getRho(index)));
                FloatMatrix elementMassD=ansys_solver->getElementMass(index,ansys_solver->getMassRhoD(ansys_solver->getRho(index)));
                FloatArray elementTemp;
                Evec.GetSubVector(edof,elementTemp);
                FloatArray eKDP=elementStifD*elementTemp;
                FloatArray eMDP=elementMassD*elementTemp;
                FloatArray TM = eKDP-Eval*eMDP;
                double lambdaD=TM*elementTemp;
                Result[index]=-1.0*lambdaD;
            }
            else {
                Result[index]=0.0;
            }
        }
    }
    if (repeateLocalindex!=0) {
        double Eval = ansys_solver->getEigenValues()[repeateLocalindex-1];
        FloatMatrix Evec(NumberOfDof, repeatedGroups[repeateLocalindex-1].giveSize());
        for (int j = 0; j < repeatedGroups[repeateLocalindex-1].giveSize(); j++) {
            auto LocalEvec = ansys_solver->getEigenVector(repeatedGroups[repeateLocalindex-1](j));
            Evec.setColumn(LocalEvec, repeatedGroups[repeateLocalindex-1](j) + 1);
        }
        auto Sub_Eigen_Problem = DifferentiableEigenvectorCalculation(Eval, Evec, ansys_solver);
        for (const auto& elem:Sub_Eigen_Problem) {
            Result(elem.first)=-1.0*elem.second[repeateGroupindex-1].SubEigenValue[repeateLocalindex-1];
        }
    }
    return Result;
}


std::map<int,std::vector<SubEigenProblem>> TopProblem::DifferentiableEigenvectorCalculation(double Eval, FloatMatrix Evec,
        std::shared_ptr<AnsysSolver>& ansys_solver) {
    std::map<int, std::vector<SubEigenProblem> > result;
    IntArray ColIndex(Evec.giveNumberOfColumns());
    FloatMatrix zjz_m(Evec.giveNumberOfColumns());
    FloatMatrix zjz_k(Evec.giveNumberOfColumns());
    for (int i = 0; i < Evec.giveNumberOfColumns(); i++) {
        ColIndex[i] = i;
    }


    for (int index = 0; index < ansys_solver->getElementSize(); index++) {
        if (ansys_solver->ansys.AnasysElements[index+1]->getTopStatus()==true) {
            FloatMatrix guodu_mx(Evec.giveNumberOfRows(), Evec.giveNumberOfColumns());
            FloatMatrix guodu_kx(Evec.giveNumberOfRows(), Evec.giveNumberOfColumns());
            IntArray edof = ansys_solver->getElementLocation(index);
            FloatMatrix elementStifD = ansys_solver->getElementStif(
                index, ansys_solver->getStifRhoD(ansys_solver->getRho(index)));
            FloatMatrix elementMassD = ansys_solver->getElementMass(
                index, ansys_solver->getMassRhoD(ansys_solver->getRho(index)));
            FloatMatrix elementTemp(Evec.giveNumberOfColumns(), Evec.giveNumberOfColumns());
            elementTemp.beSubMatrixOf(Evec, edof, ColIndex);
            FloatMatrix temp = elementStifD * elementTemp;
            guodu_kx.setSubMatrix(temp, edof, ColIndex);
            temp = elementMassD * elementTemp;
            guodu_mx.setSubMatrix(temp, edof, ColIndex);
            zjz_m.beTProductOf(Evec, guodu_mx);
            zjz_m.times(Eval);
            zjz_k.beTProductOf(Evec, guodu_kx);
            zjz_k -= zjz_m;
            FloatArray lambda_dao = zjz_k.ComputeEigenValueAndEigenVector();
            SubEigenProblem sub_eigen_problem(lambda_dao, zjz_k);
            result[index].push_back(sub_eigen_problem);
        }
        else {
            zjz_k=0.0;
            FloatArray lambda_dao(Evec.giveNumberOfColumns());
            lambda_dao=0.0;
            SubEigenProblem sub_eigen_problem(lambda_dao, zjz_k);
            result[index].push_back(sub_eigen_problem);
        }
    }

    return result;
}


inline std::pair<IntArray, std::vector<IntArray>>
TopProblem::ClassifiedEigenValue(const std::vector<std::pair<double, std::vector<int>>>& eigenGroups,
                                 const IntArray& sortedIndices) {
    std::unordered_set<int> indexSet(sortedIndices.begin(), sortedIndices.end()); // 存储 OMindex 索引，便于快速查询
    std::unordered_set<int> processedIndices; // 记录已分类的索引
    std::vector<IntArray> repeatedGroups; // 存储重频组
    IntArray uniqueIndices; // 存储单频索引

    // 遍历 eigenGroups 找到重频组
    for (const auto &entry: eigenGroups) {
        const std::vector<int> &indices = entry.second;
        IntArray group;

        // 只保留 sortedIndices 里存在的元素
        for (int num: indices) {
            if (indexSet.count(num)) {
                // O(1) 查询
                group.append(num);
                processedIndices.insert(num); // 记录已处理
            }
        }

        // 只有至少有两个元素时才算重频组
        if (group.giveSize() >= 2) {
            repeatedGroups.push_back(group);
        }

        if (group.giveSize() == 1) {
            uniqueIndices.append(group[0]);
        }
    }
    return {uniqueIndices, repeatedGroups};
}

#endif //MAIN_H
