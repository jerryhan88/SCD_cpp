//
//  ILP.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 17/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/Base.hpp"

Solution* ILP::solve() {
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
            "\"{numRows: %ld,numCols: %ld}\"",
            baseCplex->getNrows(),
            baseCplex->getNcols());
    sol->note = std::string(note);
    //
    for (int a: prob->A) {
        for (int k: prob->K) {
            sol->y_ak[a][k] = baseCplex->getValue(y_ak[a][k]);
        }
        for (int e: prob->E_a[a]) {
            for (int k: prob->K) {
                sol->z_aek[a][e][k] = baseCplex->getValue(z_aek[a][e][k]);
            }
            for (int i: prob->N_ae[a][e]) {
                for (int j: prob->N_ae[a][e]) {
                    sol->x_aeij[a][e][i][j] = baseCplex->getValue(x_aeij[a][e][i][j]);
                }
                sol->u_aei[a][e][i] = baseCplex->getValue(u_aei[a][e][i]);
            }
        }
    }
    return sol;
}
