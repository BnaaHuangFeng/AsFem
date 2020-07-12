//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.12
//+++ Purpose: define the [nonlinearsolver] block for our input file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <string>

#include "NonlinearSolver/NonlinearSolverType.h"


using namespace std;

class NonlinearSolverBlock{
public:
    NonlinearSolverBlock(){
        _SolverTypeName="newton with line search";
        _SolverType=NonlinearSolverType::NEWTONLS;
        _MaxIters=20;
        _RAbsTol=2.5e-8;
        _RRelTol=1.0e-9;
        _STol=1.0e-16; // |dx|<|x|*stol
        _PCTypeName="lu";
    }

    string              _SolverTypeName;
    NonlinearSolverType _SolverType;
    int                 _MaxIters;
    double _RAbsTol,_RRelTol,_STol;

    string _PCTypeName;

    void Init(){
        _SolverTypeName="newton with line search";
        _SolverType=NonlinearSolverType::NEWTONLS;
        _MaxIters=20;
        _RAbsTol=2.5e-8;
        _RRelTol=1.0e-9;
        _STol=1.0e-16; // |dx|<|x|*stol
        _PCTypeName="lu";
    }
};