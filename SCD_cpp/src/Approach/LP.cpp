//
//  LP.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 20/5/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/Solver.hpp"

Solution* LP::solve() {
    unsigned long TiLim = LP_TIME_LIMIT < time_limit_sec ? LP_TIME_LIMIT : time_limit_sec;
    baseCplex->setParam(IloCplex::TiLim, TiLim);
    baseCplex->solve();
    //
    Solution *sol = new Solution(prob);
    if (baseCplex->getStatus() == IloAlgorithm::Infeasible) {
        baseCplex->exportModel(lpPath.c_str());
    }
    sol->objV = baseCplex->getObjValue();
    sol->gap = -1.0;
    sol->cpuT = tt->get_elapsedTimeCPU();
    sol->wallT = tt->get_elapsedTimeWall();
    char note[2048];
    sprintf(note,
            "\"{\'numRows\': %ld, \'numCols\': %ld}\"",
            baseCplex->getNrows(),
            baseCplex->getNcols());
    sol->note = std::string(note);
    return sol;
}
