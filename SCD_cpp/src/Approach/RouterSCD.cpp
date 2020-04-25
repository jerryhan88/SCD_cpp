//
//  RouterSCD.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 25/3/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/RouterSCD.hpp"


void update_objF(rmm::RouteMM *rutMM, Problem *prob, int a, int e, double ***lrh_l_aek) {
    IloExpr objF(rutMM->env);
    int k = 0;
    for(std::set<int>::iterator it = prob->K_ae[a][e].begin(); it != prob->K_ae[a][e].end(); ++it) {
        int j = prob->RP_ae[a][e]->n_k[k];
        for (int i: prob->RP_ae[a][e]->N) {
            objF += lrh_l_aek[a][e][*it] * rutMM->x_ij[i][j];
        }
        k++;
    }
    if (rutMM->cc != nullptr) {
        rutMM->cc->bestRelVal = DBL_MAX;
    }
    rutMM->cplex->getObjective().setExpr(objF);
    objF.end();
}

void run_rutAlgo(rmm::RouteMM *rutMM, TimeTracker *tt, unsigned long time_limit_sec) {
    unsigned long timeLimit4RutSolving = time_limit_sec;
    timeLimit4RutSolving -= (unsigned long) tt->get_elapsedTimeWall();
    rutMM->cplex->setParam(IloCplex::TiLim, timeLimit4RutSolving);
    rutMM->run();
    if (rutMM->cplex->getStatus() == IloAlgorithm::Infeasible) {
        throw "the rutMM is infeasible";
    }
}

void getSol_rutMM(rmm::RouteMM *rutMM, double *rut_objV, double **rut_x_ij, double *rut_u_i){
    if (rutMM->cplex->getStatus() == IloAlgorithm::Optimal) {
        *rut_objV = -rutMM->cplex->getObjValue();
        for (int i: rutMM->prob->N) {
            for (int j: rutMM->prob->N) {
                rut_x_ij[i][j] = rutMM->cplex->getValue(rutMM->x_ij[i][j]);
            }
            rut_u_i[i] = rutMM->cplex->getValue(rutMM->u_i[i]);
        }
    } else {
        throw "the rutMM is not solved optimally";
    }
}


void RouterILP::update() {
    update_objF(rutAlgo, prob, a, e, lrh_l_aek);
}

void RouterILP::solve() {
    run_rutAlgo(rutAlgo, tt, time_limit_sec);
}

void RouterILP::getSol() {
    getSol_rutMM(rutAlgo, &rut_objV, rut_x_ij, rut_u_i);
}

void RouterBnC::update() {
    update_objF(rutAlgo, prob, a, e, lrh_l_aek);
    rutAlgo->clear_detectedCuts();
    //
    gh->initIH();
    int k = 0;
    for(std::set<int>::iterator it = prob->K_ae[a][e].begin(); it != prob->K_ae[a][e].end(); ++it) {
        prob->RP_ae[a][e]->r_k[k] = lrh_l_aek[a][e][*it];
        k++;
    }
}

void RouterBnC::solve() {
    gh->run();
    rutAlgo->set_initSol(gh->x_ij, gh->u_i);
    run_rutAlgo(rutAlgo, tt, time_limit_sec);
}

void RouterBnC::getSol() {
    getSol_rutMM(rutAlgo, &rut_objV, rut_x_ij, rut_u_i);
}

void RouterBnCoc::update() {
    update_objF(rutAlgo, prob, a, e, lrh_l_aek);
    rutAlgo->add_detectedCuts2MM();
    //
    gh->initIH();
    int k = 0;
    for(std::set<int>::iterator it = prob->K_ae[a][e].begin(); it != prob->K_ae[a][e].end(); ++it) {
        prob->RP_ae[a][e]->r_k[k] = lrh_l_aek[a][e][*it];
        k++;
    }
}

void RouterGH::update() {
    rutAlgo->initIH();
    //
    int k = 0;
    for(std::set<int>::iterator it = prob->K_ae[a][e].begin(); it != prob->K_ae[a][e].end(); ++it) {
        prob->RP_ae[a][e]->r_k[k] = lrh_l_aek[a][e][*it];
        k++;
    }
}

void RouterGH::solve() {
    rutAlgo->run();
}

void RouterGH::getSol() {
    rut_objV = -rutAlgo->get_objV();
    for (int i: rutAlgo->prob->N) {
        for (int j: rutAlgo->prob->N) {
            rut_x_ij[i][j] = rutAlgo->x_ij[i][j];
        }
        rut_u_i[i] = rutAlgo->u_i[i];
    }
}
