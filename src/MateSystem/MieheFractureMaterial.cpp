//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.12.30
//+++ Purpose: Calculate the material properties required by Miehe's
//+++          phase field fracture model
//+++           1) viscosity
//+++           2) Gc
//+++           3) L
//+++           4) H
//+++           5) dHdstrain
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/BulkMateSystem.h"
#include "Utils/MathFuns.h"

void BulkMateSystem::MieheFractureMaterial(const int &nDim, const double &t, const double &dt,
                                          const vector<double> &InputParams, const Vector3d &gpCoord,
                                          const vector<double> &gpU, const vector<double> &gpV,
                                          const vector<Vector3d> &gpGradU, const vector<Vector3d> &gpGradV,
                                          vector<double> &gpHist, const vector<double> &gpHistOld) {

    //*****************************************************************************
    //*** just to get rid of warnings, normal users dont need to do this
    //*****************************************************************************
    if(nDim||t||dt||InputParams.size()||
       gpCoord(1)||gpU.size()||gpV.size()||gpGradU.size()||gpGradV.size()||
       gpHist.size()||gpHistOld.size()){}

    if(InputParams.size()<5){
        MessagePrinter::PrintErrorTxt("for Miehe's phase field fracture materials, at least 5 parameters are required, you need to give: lambda, mu, Gc, L, viscosity");
        MessagePrinter::AsFem_Exit();
    }

    int UseHist=0;
    double d;
    double g,dg;// for the degradation function
    const double k=1.0e-6; // for stabilizer

    const double lambda=InputParams[0];
    const double mu=InputParams[1];
    const double Gc=InputParams[2];
    const double L=InputParams[3];
    const double viscosity=InputParams[4];


    _ScalarMaterials["viscosity"]=viscosity;
    _ScalarMaterials["Gc"]=Gc;
    _ScalarMaterials["L"]=L;


    UseHist=0;
    if(InputParams.size()>=6){
        UseHist=static_cast<int>(InputParams[6-1]);
        if(UseHist<0) UseHist=0;
    }


    d=gpU[1];
    g=(1-d)*(1-d);// degradation
    dg=2*(d-1);   // derivative of g

    RankTwoTensor Eps,EpsPos,EpsNeg;
    RankFourTensor ProjPos,ProjNeg,I4Sym;
    RankTwoTensor EigVec;
    double EigVal[3];
    if(nDim==1){
        MessagePrinter::PrintErrorTxt("Miehe's phase field fracture model works only for 2D and 3D case");
        MessagePrinter::AsFem_Exit();
    }
    else if(nDim==2){
        Eps.SetFromGradU(gpU[2],gpU[3]);
    }
    else if(nDim==3){
        Eps.SetFromGradU(gpU[2],gpU[3],gpU[4]);
    }

    ProjPos=Eps.CalcPostiveProjTensor(EigVal,EigVec);
    I4Sym.SetToIdentitySymmetric4();
    ProjNeg=I4Sym-ProjPos;

    // for the positive and negative strain
    EpsPos=ProjPos.DoubleDot(Eps);
    EpsNeg=Eps-EpsPos;

    double trEps,signpos,signneg;
    double psi,psipos,psineg;


    trEps=Eps.Trace();
    psipos=0.5*lambda*BracketPos(trEps)*BracketPos(trEps)+mu*(EpsPos*EpsPos).Trace();
    psineg=0.5*lambda*BracketNeg(trEps)*BracketNeg(trEps)+mu*(EpsNeg*EpsNeg).Trace();
    psi=(g+k)*psipos+psineg;

    _ScalarMaterials["Psi"]=psi;
    _ScalarMaterials["PsiPos"]=psipos;
    _ScalarMaterials["PsiNeg"]=psineg;

    RankTwoTensor StressPos,StressNeg,I;


    StressPos=lambda*BracketPos(trEps)*I+2*mu*EpsPos;
    StressNeg=lambda*BracketNeg(trEps)*I+2*mu*EpsNeg;

    _Rank2Materials["Stress"]=(g+k)*StressPos+StressNeg;
    _Rank2Materials["dStressdD"]=dg*StressPos;

    if(psipos>gpHistOld[0]){
        _ScalarMaterials["Hist"]=psipos;
        _Rank2Materials["dHdstrain"]=StressPos;
    }
    else{
        _ScalarMaterials["Hist"]=gpHistOld[0];
        _Rank2Materials["dHdstrain"].SetToZeros();
    }

    if(UseHist){
        _ScalarMaterials["Hist"]=gpHistOld[0];
        _Rank2Materials["dHdstrain"].SetToZeros();
    }

    signpos=0.0;
    if(BracketPos(trEps)>0) signpos=1.0;

    signneg=0.0;
    if(BracketNeg(trEps)<0) signneg=1.0;

    _Rank4Materials["elasticity_tensor"]=(g+k)*(lambda*signpos*I.CrossDot(I)+2*mu*ProjPos)
                                        +lambda*signneg*I.CrossDot(I)+2*mu*ProjNeg;

}

