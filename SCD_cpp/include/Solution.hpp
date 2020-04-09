//
//  Solution.hpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 18/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef Solution_hpp
#define Solution_hpp

#include <sstream>

#include "Problem.hpp"

class Solution {
public:
    Problem *prob;
    //
    double objV, gap;
    double cpuT, wallT;
    std::string note;
    int lastUpdatedIter;
    //
    double **y_ak, ***z_aek, ****x_aeij, ***u_aei;
    //
    Solution(Problem *prob) {
        this->prob = prob;
        //
        size_t numAgents, numTasks, numNodes;
        numAgents = prob->A.size();
        numTasks = prob->K.size();
        numNodes = prob->cN.size();
        //
        y_ak = new double*[numAgents];
        z_aek = new double**[numAgents];
        x_aeij = new double***[numAgents];
        u_aei = new double**[numAgents];
        for (int a: prob->A) {
            y_ak[a] = new double[numTasks];
            size_t numRR = prob->E_a[a].size();
            z_aek[a] = new double*[numRR];
            x_aeij[a] = new double**[numRR];
            u_aei[a] = new double*[numRR];
            for (int e: prob->E_a[a]) {
                z_aek[a][e] = new double[numTasks];
                x_aeij[a][e] = new double*[numNodes];
                u_aei[a][e] = new double[numNodes];
                for (int i = 0; i < numNodes; i++) {
                    x_aeij[a][e][i] = new double[numNodes];
                }
            }
        }
    }
    ~Solution() {
        for (int a: prob->A) {
            delete [] y_ak[a];
            for (int e: prob->E_a[a]) {
                delete [] z_aek[a][e];
                delete [] u_aei[a][e];
                for (int i = 0; i < prob->cN.size(); i++) {
                    delete [] x_aeij[a][e][i];
                }
                delete [] x_aeij[a][e];
            }
            delete [] z_aek[a];
            delete [] u_aei[a];
            delete [] x_aeij[a];
        }
        delete [] y_ak;
        delete [] z_aek;
        delete [] x_aeij;
        delete [] u_aei;
    }
    //
    void writeSolCSV(std::string solPathCSV);
    void writeSolTXT(std::string solPathTXT);
};

#endif /* Solution_hpp */
