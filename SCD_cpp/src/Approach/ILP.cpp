//
//  ILP.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 17/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/OtherSolvers.hpp"

Solution* ILP::solve() {
    baseCplex->setParam(IloCplex::TiLim, time_limit_sec);
    baseCplex->solve();
    //
    Solution *sol = new Solution(prob);
    if (baseCplex->getStatus() == IloAlgorithm::Infeasible) {
        baseCplex->exportModel(lpPath.c_str());
    }
    sol->objV = baseCplex->getObjValue();
    sol->gap = baseCplex->getMIPRelativeGap();
    sol->cpuT = tt->get_elapsedTimeCPU();
    sol->wallT = tt->get_elapsedTimeWall();
    char note[2048];
    sprintf(note,
            "\"{\'numRows\': %ld, \'numCols\': %ld}\"",
            baseCplex->getNrows(),
            baseCplex->getNcols());
    sol->note = std::string(note);
    //
    for (int a: prob->A) {
        for (int k: prob->K) {
            sol->y_ary[a][k] = baseCplex->getValue(y_ary[a][k]);
        }
        for (int e: prob->E_a[a]) {
            for (int k: prob->K) {
                sol->z_ary[a][e][k] = baseCplex->getValue(z_ary[a][e][k]);
            }
            for (int i: prob->N_ae[a][e]) {
                for (int j: prob->N_ae[a][e]) {
                    if (bool_x_aeij[a][e][i][j]) {
                        sol->x_ary[a][e][i][j] = baseCplex->getValue(x_ary[a][e][i][j]);
                    } else {
                        sol->x_ary[a][e][i][j] = 0.0;
                    }
                }
                sol->u_ary[a][e][i] = baseCplex->getValue(u_ary[a][e][i]);
            }
        }
    }
    return sol;
}
