//
//  Base.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 17/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "Base.hpp"

Base::Base(Problem *prob, std::string logPath, TimeTracker *tt) {
    this->prob = prob;
    this->logPath = logPath;
    this->tt = tt;
}
