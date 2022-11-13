//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.11.13
//+++ Purpose: Defines the abstract class for plastic materials
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MathUtils/Vector3d.h"
#include "MathUtils/Rank2Tensor.h"
#include "MathUtils/Rank4Tensor.h"

#include "ElmtSystem/LocalElmtData.h"
#include "MateSystem/MaterialsContainer.h"

#include "nlohmann/json.hpp"

/**
 * This abstract class defines the necessary functions for elasto-plastic materials
*/
class PlasticMaterialBase{
protected:
    /**
     * Evaluate the value of the yield function
     * @param t_parameters the input json parameters read from input file
     * @param t_args_vars the input variables, i.e., effective plastic strain, etc.
     * @param t_args_stress the input stress, i.e. cauchy stress, Mandel stress, etc.
     * @param t_args_strain the input strain, i.e. plastic strain, etc.
    */
    virtual double computeYieldFunction(const nlohmann::json &t_parameters,
                                        const vector<double> &t_args_vars,
                                        const vector<Rank2Tensor> &t_args_stress,
                                        const vector<Rank2Tensor> &t_args_strain)=0;
    /**
     * Evaluate the 1st order derivative value of the yield function
     * @param t_parameters the input json parameters read from input file
     * @param t_args_vars the input variables, i.e., effective plastic strain, etc.
     * @param t_args_stress the input stress, i.e. cauchy stress, Mandel stress, etc.
     * @param t_args_strain the input strain, i.e. plastic strain, etc.
    */
    virtual double computeYieldFunctionDeriv(const nlohmann::json &t_parameters,
                                             const vector<double> &t_args_vars,
                                             const vector<Rank2Tensor> &t_args_stress,
                                             const vector<Rank2Tensor> &t_args_strain)=0;
    
    /**
     * Evalute the admissible stresses status during the elasto-plastic deformation, 
     * where one can use the radial return mapping for the calculation.
     * @param t_parameters the input json parameters read from input file
     * @param t_elmtinfo the local element info structure
     * @param t_elmtsoln the local element solution 
     * @param t_mateold the old material constain from previous step
     * @param t_total_strain the total strain or total deformation gradient tensor
     * @param t_mate the current material container
     * @param t_args_stress the output stresses, i.e. cauchy stress, Mandel stress, etc.
     * @param t_args_strain the output strain, i.e. plastic strain, cauchy strain, etc.
     * @param t_args_jacobian the output jacobian and other rank-4 tensor
    */
    virtual void computeAdmissibleStressState(const nlohmann::json &t_parameters,
                                              const LocalElmtInfo &t_elmtinfo,
                                              const LocalElmtSolution &t_elmtsoln,
                                              const MaterialsContainer &t_mateold,
                                              const Rank2Tensor &t_total_strain,
                                              MaterialsContainer &t_mate,
                                              vector<Rank2Tensor> &t_args_stress,
                                              vector<Rank2Tensor> &t_args_strain,
                                              vector<Rank4Tensor> &t_args_jacobian)=0;
};
