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
    //
    double **y_ak, ***z_aek, ****x_aeij, ***u_aei;
    //
    Solution(Problem *prob);
    ~Solution();
    void alloMem4dvs();
    void writeSolCSV(std::string solPathCSV);
    void writeSolTXT(std::string solPathTXT);
};

#endif /* Solution_hpp */
