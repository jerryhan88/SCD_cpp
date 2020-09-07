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
    double **y_ary, ***z_ary, ****x_ary, ***u_ary;
    //
    Solution(Problem *prob) {
        this->prob = prob;
        //
        y_ary = new_ak_dbl(prob);
        z_ary = new_aek_dbl(prob);
        x_ary = new_aeij_dbl(prob);
        u_ary = new_aei_dbl(prob);
    }
    ~Solution() {
        del_ak_dbl(prob, y_ary);
        del_aek_dbl(prob, z_ary);
        del_aeij_dbl(prob, x_ary);
        del_aei_dbl(prob, u_ary);
    }
    //
    void writeSolCSV(std::string solPathCSV);
    void writeSolTXT(std::string solPathTXT);
};

#endif /* Solution_hpp */
