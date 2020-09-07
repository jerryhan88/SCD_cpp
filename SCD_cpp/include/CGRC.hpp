//
//  CGRC.hpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 7/9/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef CGRC_hpp
#define CGRC_hpp

#include <ilcplex/ilocplex.h>

#include "SolApprBase.hpp"

#include "Problem.hpp"
#include "Solution.hpp"


class CGRC : public SolApprBase {
public:
    std::string lpPath;
    std::string pathPrefix;
    std::string lp_algo;
    //
    IloEnv env;
    IloModel *baseModel;
    IloCplex *baseCplex;
    IloNumVar **y_ary, ***z_ary, ****x_ary, ***u_ary;
    //
    CGRC(Problem *prob, TimeTracker *tt, unsigned long time_limit_sec, int numThreads,
        std::string logPath, std::string lpPath,
        std::string lp_algo): SolApprBase(prob, tt, time_limit_sec, numThreads, logPath) {
        
        this->lpPath = lpPath;
        this->pathPrefix = lpPath.substr(0, lpPath.find(".lp"));
        this->lp_algo = lp_algo;
        

    }
    
};


#endif /* CGRC_hpp */
