//
//  Allocator.hpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 25/4/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef Allocator_h
#define Allocator_h

#include <ilcplex/ilocplex.h>

#include "SolApprBase.hpp"

class Allocator {
public:
    Problem *prob;
    //
    IloEnv env;
    IloModel *etaModel; // Evaluation on Task Assingment
    IloCplex *etaCplex;
    IloNumVar **eta_y_ak, ***eta_z_aek;
    //
    IloNum minGap, lastUpdatedST;
    IloBool aborted;
    //
    double ***lrh_l_aek;
    //
    Allocator(Problem *prob, TimeTracker *tt, double ***lrh_l_aek);
    ~Allocator() {
        delete_inv_ak(prob, eta_y_ak);
        delete_inv_aek(prob, eta_z_aek);
        //
        delete etaCplex; delete etaModel;
        env.end();
    }
    void build();
    void update();
    void getSol(double *L1_V, double **lrh_y_ak, double ***lrh_z_aek);
};


#endif /* Allocator_h */
