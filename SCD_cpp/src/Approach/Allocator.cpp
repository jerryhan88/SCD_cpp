//
//  Allocator.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 25/4/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/Allocator.hpp"

 
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
//    objF = new IloExpr(env);
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
}

void Allocator::getSol(double *L1_V, double **lrh_y_ak, double ***lrh_z_aek) {
    if (etaCplex->getStatus() == IloAlgorithm::Optimal) {
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
    } else {
        throw "the etaModel is not solved optimally";
    }
    
}
