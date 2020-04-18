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
        y_ak = new_dbl_ak(prob);
        z_aek = new_dbl_aek(prob);
        x_aeij = new_dbl_aeij(prob);
        u_aei = new_dbl_aei(prob);
    }
    ~Solution() {
        delete_dbl_ak(prob, y_ak);
        delete_dbl_aek(prob, z_aek);
        delete_dbl_aeij(prob, x_aeij);
        delete_dbl_aei(prob, u_aei);
    }
    //
    void writeSolCSV(std::string solPathCSV);
    void writeSolTXT(std::string solPathTXT);
};

#endif /* Solution_hpp */
