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
            objF += -lrh_l_aek[a][e][*it] * rutMM->x_ij[i][j];
        }
        k++;
    }
    rutMM->cplex->getObjective().setExpr(objF);
    if (rutMM->cplex->getObjective().getSense() != IloObjective::Minimize) {
        rutMM->cplex->getObjective().setSense(IloObjective::Minimize);
    }
    objF.end();
}

void run_rutAlgo(rmm::RouteMM *rutMM, double *rut_objV, double **rut_x_ij, double *rut_u_i) {
    rutMM->run();
    if (rutMM->cplex->getStatus() == IloAlgorithm::Infeasible) {
        throw "the rutMM is infeasible";
    }
    if (rutMM->cplex->getStatus() == IloAlgorithm::Optimal) {
        *rut_objV = rutMM->cplex->getObjValue();
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

void RouterILP::run() {
    rutAlgo->run();
    run_rutAlgo(rutAlgo, &rut_objV, rut_x_ij, rut_u_i);
}


void RouterBnC::update() {
    update_objF(rutAlgo, prob, a, e, lrh_l_aek);
    if (turOnCutPool) {
        rutAlgo->add_detectedCuts2MM();
    }
}

void RouterBnC::run() {
    rutAlgo->run();
    run_rutAlgo(rutAlgo, &rut_objV, rut_x_ij, rut_u_i);
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

void RouterGH::run() {
    rutAlgo->run();
    //
    rut_objV = -rutAlgo->get_objV();
    for (int i: rutAlgo->prob->N) {
        for (int j: rutAlgo->prob->N) {
            rut_x_ij[i][j] = rutAlgo->x_ij[i][j];
        }
        rut_u_i[i] = rutAlgo->u_i[i];
    }
}
