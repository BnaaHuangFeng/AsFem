//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.11.14
//+++ Purpose: calculate the material properties for diffusion coupled 
//+++          phase-field fracture solids 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once

#include "MateSystem/MultiphycisMechanicsMaterialBase.h"

/**
 * calculate the constitutive laws for diffusion coupled neohookean hyperlelastic material.
 */
class DiffusionFractureMaterial: public MultiphysicsMechanicsMaterialBase{
public:
    
    /**
     * initialize the material properties
     */
    virtual void InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate) override;
    /**
     * compute the related material properties
     */
    virtual void ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate) override;

private: 
    /**
     * Compute the deformation gradient tensor, it could be small strain \f$\mathbf{\varepsilon}\f$, finite strain deformation gradient tensor \f$\mathbf{F}=\mathbf{U}+\mathbf{I}\f$.
     * @param elmtinfo the current element information data
     * @param elmtsoln the curent element's solution, include the displacement and its gradient
     * @param F the rank-2 strain tensor
     */ 
    virtual void ComputeDeformationGradientTensor(const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, RankTwoTensor &F) override;

    /**
     * Compute the stress \f$\mathbf{\sigma}\f$ and jacobian matrix \f$\mathbb{C}\f$ .
     * @param InputParams the input material parameters from your input file.
     * @param soln the local solution structure.
     * @param Strain the input strain tensor, it is calculated from ComputeStrain function.
     * @param Stress the calculated stress tensor \f$\mathbf{\sigma}\f$ for MechanicsElmt.
     * @param Jacobian the calculated 'elasticity' tensor \f$\mathbb{C}\f$, it is the rank-4 tensor for the general constitutive law.
     */
    virtual void ComputeConstitutiveLaws(const vector<double> &InputParams, const LocalElmtSolution &soln, const RankTwoTensor &F, RankTwoTensor &Strain, RankTwoTensor &Stress, RankFourTensor &Jacobian) override;

private:
    double DegradationFun(const double &x){
        return (1-x)*(1-x)+1.0e-5;
    }
    double DegradationFunDeriv(const double &x){
        return 2*(x-1);
    }

private:
    RankTwoTensor _I,_F,_Strain,_devStress;
    RankTwoTensor _GradU,_MechStrain,_dMechStraindC;
    RankTwoTensor _EigenStrain,_dEigenStraindC;
    RankTwoTensor _EpsPos,_EpsNeg;
    RankTwoTensor _Stress,_StressPos,_StressNeg;
    RankTwoTensor _dStressdC,_dStressdD,_dStressPosdC,_dStressNegdC;
    RankTwoTensor _dEpsPosdC,_dEpsNegdC;
    RankFourTensor _Jac,_ProjPos,_ProjNeg,_I4Sym;
    Vector3d _GradSigmaH;
    double _PsiPos,_PsiNeg,_dSigmaHdC;
    int UseHist;

};