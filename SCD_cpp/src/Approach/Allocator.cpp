//
//  Allocator.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 25/4/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/Allocator.hpp"

//#define ALLOCATIOR_TIME_LIMIT 0.00001
#define ALLOCATIOR_TIME_LIMIT 30

ILOMIPINFOCALLBACK4(timeLimitCallback,
                    IloNum,   timeLimit,
                    IloNum*,   minGap,
                    IloNum*,   lastUpdatedST,
                    IloBool*,  aborted) {

    if ( !(*aborted) && hasIncumbent()) {
        IloNum gap = getMIPRelativeGap();
        if (*minGap - gap > 0.01) {
            *minGap = gap;
            *lastUpdatedST = getCplexTime();
        }
        //
        IloNum timeUsed = getCplexTime() - *lastUpdatedST;
        if ( timeUsed > timeLimit ) {
            *aborted = IloTrue;
            abort();
        }
    }
}

Allocator::Allocator(Problem *prob, TimeTracker *tt, double ***lrh_l_aek) {
    this->prob = prob;
    //
    eta_y_ak = new_inv_ak(prob, env, 'I');
    eta_z_aek = new_inv_aek(prob, env, 'I');
    this->lrh_l_aek = lrh_l_aek;
    //
    build();
    //
    etaCplex->use(timeLimitCallback(env, ALLOCATIOR_TIME_LIMIT, &minGap, &lastUpdatedST, &aborted));
}

void Allocator::build() {
    etaModel = new IloModel(env);
    //
    IloExpr objF = eta_y_ak[0][0];  // this function is a dummy; later the objective will be updated
    etaModel->add(IloMinimize(env, objF));
    objF.end();
    //
    BaseMM::def_ETA_cnsts(prob, env, eta_y_ak, eta_z_aek, etaModel);
    //
    etaCplex = new IloCplex(*etaModel);
    etaCplex->setOut(env.getNullStream());
}

void Allocator::update() {
    // update coefficient associated with lambda values
    IloExpr objF(env);
    for (int k: prob->K) {
        for (int a : prob->A) {
            for (int e: prob->E_a[a]) {
                objF += prob->r_k[k] * prob->p_ae[a][e] * eta_z_aek[a][e][k];
            }
            objF -= prob->r_k[k] * eta_y_ak[a][k];
        }
    }
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            for (int k: prob->K_ae[a][e]) {
                objF += lrh_l_aek[a][e][k] * eta_y_ak[a][k];
                objF -= lrh_l_aek[a][e][k] * eta_z_aek[a][e][k];
            }
        }
    }
    etaCplex->getObjective().setExpr(objF);
    objF.end();
    //
    minGap = DBL_MAX;
    lastUpdatedST = 0.0;
    aborted = IloFalse;
}

void Allocator::getSol(double *L1_V, double **lrh_y_ak, double ***lrh_z_aek) {
    try {
        *L1_V = -etaCplex->getObjValue();
        for (int a: prob->A) {
            for (int k: prob->K) {
                lrh_y_ak[a][k] = etaCplex->getValue(eta_y_ak[a][k]);
            }
            for (int e: prob->E_a[a]) {
                for (int k: prob->K) {
                    lrh_z_aek[a][e][k] = etaCplex->getValue(eta_z_aek[a][e][k]);
                }
            }
        }
    } catch (IloException& e) {
        throw "the etaModel is not solved";
    }
    if (etaCplex->getStatus() != IloAlgorithm::Optimal) {
        std::cout << "\tthe etaModel is not solved, but the solution is not optimal" << std::endl;
    }
}
