//
//  SolBase.hpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 29/3/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef SolApprBase_hpp
#define SolApprBase_hpp

#include "Problem.hpp"
#include "Solution.hpp"
#include "ThreadPool.hpp"

#include "ck_util/util.hpp"         // from util

class SolApprBase {
public:
    Problem *prob;
    TimeTracker *tt;
    unsigned long time_limit_sec;
    std::string logPath;
    ThreadPool& pool = ThreadPool::getInstance(1);
    //
    SolApprBase(Problem *prob, TimeTracker *tt, unsigned long time_limit_sec, int numThreads, std::string logPath) {
        this->prob = prob;
        this->logPath = logPath;
        this->tt = tt;
        this->time_limit_sec = time_limit_sec;
        if (numThreads > 1) {
            pool.resize(numThreads);
        }
    }
    ~SolApprBase() {}
    virtual Solution* solve() {
        throw "Should override solve()";
    }
};
#endif /* SolBase_h */
