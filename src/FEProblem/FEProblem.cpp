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
//+++ Date   : 2022.05.09
//+++ Purpose: the FEProblem class of AsFem, the top level of the
//+++          whole program
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FEProblem/FEProblem.h"

FEProblem::FEProblem(){
    m_timer.resetTimer();
}

void FEProblem::initFEProblem(int args,char *argv[]){
    //***************************************
    // for input file reading
    //***************************************
    m_timer.startTimer();
    m_inputSystem.init(args,argv);

    MessagePrinter::printStars();
    MessagePrinter::printNormalTxt("Start to read the input file");
    m_inputSystem.readInputFile(m_mesh,m_dofhandler,m_elmtsystem,m_fe,
                                m_bcsystem,m_icsystem,
                                m_projsystem,
                                m_nlsolver,
                                m_timestepping,
                                m_output,
                                m_postprocessor,
                                m_jobblock);
    m_timer.endTimer();
    m_timer.printElapseTime("Input file reading is done",false);
    m_mesh.printBulkMeshInfo();

    if(m_inputSystem.isReadOnly()) return;

    //***************************************
    // for dofs init 
    // set dofhandler -> total number of elmts, nodes, dofs data,
    // allocate memory for map of nodal and elemental local dof data id 2 global dof data id. 
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to create dofs map ...");
    m_dofhandler.createBulkDofsMap(m_mesh,m_elmtsystem);
    m_timer.endTimer();
    m_timer.printElapseTime("Dofs map generation is done",false);
    m_dofhandler.printBulkDofsInfo();

    //***************************************
    // for elmt system init
    // assign each elmt element block (m_elemental_elmtblock_id) 
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the Element system ...");
    m_elmtsystem.init(m_mesh);  // Hf: assign element type to every element
    m_timer.endTimer();
    m_timer.printElapseTime("Element system is initialized",false);

    //***************************************
    // for bc system init
    // allocate room for BCSystem -> m_nodes0,m_nodes LocalElmtInfo, m_localR, m_localK
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the BC system ...");
    m_bcsystem.init(m_dofhandler.getMaxDofsPerNode());
    m_timer.endTimer();
    m_timer.printElapseTime("BC system is initialized",false);

    //***************************************
    // for ic system init
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the IC system ...");
    m_icsystem.init(m_dofhandler.getMaxDofsPerNode());
    m_timer.endTimer();
    m_timer.printElapseTime("IC system is initialized",false);

    //***************************************
    // for FE space init
    // allocate room for shapefun (data are still 0), set qpoint (specify data)
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the FE space ...");
    m_fe.init(m_mesh);
    m_timer.endTimer();
    m_timer.printElapseTime("FE space is initialized",false);

    //***************************************
    // for FE system init
    // preallocate BulkFESystem's member, but data are 0 now.
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the FE system ...");
    m_fesystem.init(m_mesh,m_dofhandler);
    m_timer.endTimer();
    m_timer.printElapseTime("FE system is initialized",false);

    //***************************************
    // for Equation system init
    // preallocate m_rhs & m_amatrix (datas are 0)
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the Equation system ...");
    m_equationsystem.init(m_dofhandler);
    MessagePrinter::printNormalTxt("  Start to create Sparsity pattern ...");
    m_equationsystem.createSparsityPattern(m_dofhandler);
    MessagePrinter::printNormalTxt("  Sparsity pattern is ready");
    m_timer.endTimer();
    m_timer.printElapseTime("Equation system is initialized",false);

    //***************************************
    // for Solution system init
    // preallocate solution datas (all dofs) & material datas(all gauss points) (datas are 0)
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the Solution system ...");
    m_solutionsystem.init(m_dofhandler,m_fe);
    m_timer.endTimer();
    m_timer.printElapseTime("Solution system is initialized",false);

    //***************************************
    // for Projection system init
    // preallocate solution datas (a node) & material datas(all nodes) (datas are 0)
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the Projection system ...");
    m_projsystem.init(m_mesh,m_dofhandler);
    m_timer.endTimer();
    m_timer.printElapseTime("Projection system is initialized",false);

    //***************************************
    // for Nonlinear solver system init
    // set SNES, KSP, PC
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the NL solver ...");
    m_nlsolver.init();
    m_timer.endTimer();
    m_timer.printElapseTime("NL solver is initialized",false);

    //***************************************
    // for Postprocess system init
    // preallocate Nodes (in a element) (datas are 0)
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the postprocessor ...");
    m_postprocessor.init();
    m_timer.endTimer();
    m_timer.printElapseTime("Postprocessor is initialized",false);

    //***************************************
    // for fe control init
    // set default ctan, t, dt, timesteppingtype, startstep, finalstep...
    //***************************************
    m_fectrlinfo.init();
    m_fectrlinfo.IsDebug=m_jobblock.m_isdebug;
    m_fectrlinfo.IsDepDebug=m_jobblock.m_isdepdebug;


    //***************************************
    // for print out basic info
    //***************************************
    m_elmtsystem.printElmtSystemInfo();
    m_bcsystem.printBCSystemInfo();
    m_icsystem.printICSystemInfo();
    m_projsystem.printProjectionInfo();
    m_fe.printFEInfo();
    m_nlsolver.printSolverInfo();
    m_output.printInfo();
    m_postprocessor.printInfo();
    if(m_jobblock.m_jobtype==FEJobType::TRANSIENT){
        m_timestepping.printInfo();
    }
    m_jobblock.printJobInfo();
}
//*******************************************
void FEProblem::finalize(){
    m_mesh.releaseMemory();
    m_dofhandler.releaseMemory();
    m_fe.releaseMemory();
    m_elmtsystem.releaseMemory();
    m_bcsystem.releaseMemory();
    m_icsystem.releaseMemory();
    m_equationsystem.releaseMemory();
    m_solutionsystem.releaseMemory();
    m_projsystem.releaseMemory();
    m_nlsolver.releaseMemory();
    m_postprocessor.releaseMemory();
    
}