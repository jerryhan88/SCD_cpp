//
//  Base.hpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 17/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef Base_hpp
#define Base_hpp

#include <cfloat>

#include <ilcplex/ilocplex.h>

#include "../Other/Problem.hpp"
#include "../Other/Solution.hpp"
#include "../Other/Etc.hpp"

#define DEFAULT_BUFFER_SIZE 2048

class Base {
public:
    Problem *prob;
    std::string logPath;
    TimeTracker *tt;
    //
    Base(Problem *prob, std::string logPath, TimeTracker *tt);
    virtual Solution* solve() {
        throw "Should override solve()";
    }
};
    
class ILP : public Base {
public:
    std::string lpPath;
    //
    IloEnv env;
    IloCplex *cplex;
    IloModel *cplexModel;
    IloNumVar **y_ak, ***z_aek, ****x_aeij, ***u_aei;
    //
    ILP(Problem *prob, std::string logPath, TimeTracker *tt, std::string lpPath);
    Solution* solve();
protected:
    void build_baseModel();
private:
    void def_dvs();
    void def_objF();
    void def_ETA_cnsts();
    void def_RUT_cnsts();
    void def_COM_cnsts();
    void def_FC_cnsts_aeGiven(int a, int e);
    void def_AT_cnsts_aeGiven(int a, int e);
};
    
#endif /* Base_hpp */
