//
//  PrimalExtractor.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 13/4/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/PrimalExtractor.hpp"


void RouteFixPE::build() {
    pexModel = new IloModel(env);
    //
    IloExpr objF(env);
    for (int k: prob->K) {
        for (int a : prob->A) {
            objF += prob->r_k[k] * pex_y_ak[a][k];
            for (int e: prob->E_a[a]) {
                objF -= prob->r_k[k] * prob->p_ae[a][e] * pex_z_aek[a][e][k];
            }
        }
    }
    pexModel->add(IloMaximize(env, objF));
    objF.end();
    //
    BaseMM::def_ETA_cnsts(prob, env, pex_y_ak, pex_z_aek, pexModel);
    char buf[2048];
    IloExpr linExpr(env);
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            for (int k: prob->K_ae[a][e]) {
                linExpr.clear();
                sprintf(buf, "CC(%d)(%d)(%d)", a, e, k);
                linExpr += pex_y_ak[a][k];
                double sumX = 0.0;
                linExpr -= pex_z_aek[a][e][k];
                pex_COM_cnsts->add(linExpr <= sumX);
                (*pex_COM_cnsts)[pex_COM_cnsts->getSize() - 1].setName(buf);
                pex_COM_cnsts_index[a][e][k] = pex_COM_cnsts->getSize() - 1;
            }
        }
    }
    linExpr.end();
    pexModel->add(*pex_COM_cnsts);
    pexCplex = new IloCplex(*pexModel);
    pexCplex->setOut(env.getNullStream());
    pexCplex->setWarning(env.getNullStream());
}

void RouteFixPE::update() {
    IloExpr linExpr(env);
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            for (int k: prob->K_ae[a][e]) {
                long cnsts_id = pex_COM_cnsts_index[a][e][k];
                double sumX = 0.0;
                for (int j: prob->N_ae[a][e]) {
                    sumX += lrh_x_aeij[a][e][j][prob->n_k[k]];
                }
                (*pex_COM_cnsts)[cnsts_id].setUB(sumX);
            }
        }
    }
}

void RouteFixPE::getSol(Solution *sol) {
    for (int a: prob->A) {
        for (int k: prob->K) {
            sol->y_ak[a][k] = pexCplex->getValue(pex_y_ak[a][k]);
        }
        for (int e: prob->E_a[a]) {
            for (int k: prob->K) {
                sol->z_aek[a][e][k] = pexCplex->getValue(pex_z_aek[a][e][k]);
            }
            for (int i: prob->N_ae[a][e]) {
                for (int j: prob->N_ae[a][e]) {
                    sol->x_aeij[a][e][i][j] = lrh_x_aeij[a][e][i][j];
                }
                sol->u_aei[a][e][i] = lrh_u_aei[a][e][i];
            }
        }
    }
}


void ColGenPE::build() {
    pexModel = new IloModel(env);
    //
    IloExpr objF(env);
    pexModel->add(IloMaximize(env, objF));
    objF.end();
    //
    char buf[2048];
    IloRangeArray cnsts(env);
    IloExpr linExpr(env);
    //
    for (int k : prob->K) {
        linExpr.clear();
        sprintf(buf, "TAS(%d)", k);  // Task Assignment
        for (int a : prob->A) {
            linExpr += pex_y_ak[a][k];
        }
        cnsts.add(linExpr <= 1);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            for (int k: prob->K) {
                if (prob->K_ae[a][e].find(k) == prob->K_ae[a][e].end()) {
                    linExpr.clear();
                    sprintf(buf, "IA(%d)(%d)(%d)", a, e, k);  // Infeasible Assignment
                    linExpr += pex_y_ak[a][k];
                    cnsts.add(linExpr == 0);
                    cnsts[cnsts.getSize() - 1].setName(buf);
                }
            }
        }
    }
    //
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            linExpr.clear();
            sprintf(buf, "ONE(%d)(%d)", a, e);
            pex_ONE_cnsts->add(linExpr <= 1);
            (*pex_ONE_cnsts)[pex_ONE_cnsts->getSize() - 1].setName(buf);
            pex_ONE_cnsts_index[a][e] = pex_ONE_cnsts->getSize() - 1;
        }
    }
    linExpr.end();
    pexModel->add(cnsts);
    pexModel->add(*pex_ONE_cnsts);
    pexCplex = new IloCplex(*pexModel);
    pexCplex->setOut(env.getNullStream());
    pexCplex->setWarning(env.getNullStream());
}

void ColGenPE::update() {
    char buf[2048];
    IloRangeArray cnsts(env);
    IloExpr linExpr(env);
    //
    for (int a : prob->A) {
        std::set<int> newlyAddedW;
        std::set<int> completedTasks_a;
        for (int e: prob->E_a[a]) {
            int col_index = (int) og.size();
            og.push_back(col_index);
            og_a[a].push_back(col_index);
            newlyAddedW.insert(col_index);
            og_ae[a][e].push_back(col_index);
            //
            std::map<int, int> fromToPairs;
            for (int i: prob->N_ae[a][e]) {
                for (int j: prob->N_ae[a][e]) {
                    if (lrh_x_aeij[a][e][i][j] > 0.5) {
                        fromToPairs[i] = j;
                    }
                }
            }
            int n0 = prob->o_ae[a][e];
            std::vector<int> route;
            std::set<int> completedTasks;
            while (n0 != prob->d_ae[a][e]) {
                route.push_back(n0);
                n0 = fromToPairs[n0];
            }
            route.push_back(n0);
            og_rut.push_back(route);
            std::vector<double> arT;
            for (int i: route) {
                arT.push_back(lrh_u_aei[a][e][i]);
                if (prob->dp_k[i] != -1) {
                    completedTasks.insert(prob->dp_k[i]);
                    completedTasks_a.insert(prob->dp_k[i]);
                }
            }
            og_arT.push_back(arT);
            og_tsk.push_back(completedTasks);
            //
            double expProf = 0.0;
            std::vector<int> indicator(prob->K.size(), 0);
            for (int k: completedTasks) {
                expProf += prob->p_ae[a][e] * prob->r_k[k];
                indicator[k] = 1;
            }
            p_w.push_back(expProf);
            e_wk.push_back(indicator);
            //
            pex_th_w->add(IloNumVar(env, 0.0, 1.0, ILOINT));
            sprintf(buf, "o(%d)", col_index);
            (*pex_th_w)[col_index].setName(buf);
            long cnsts_id = pex_ONE_cnsts_index[a][e];
            (*pex_ONE_cnsts)[cnsts_id].setLinearCoef((*pex_th_w)[col_index], 1.0);
        }
        for (int w: newlyAddedW) {
            for (int k: completedTasks_a) {
                if (e_wk[w][k] == 0) {
                    continue;
                }
                linExpr.clear();
                sprintf(buf, "BI(%d)(%d)(%d)", a, w, k);
                linExpr += e_wk[w][k] * (*pex_th_w)[w];
                linExpr -= pex_y_ak[a][k];
                cnsts.add(linExpr <= 0);
                cnsts[cnsts.getSize() - 1].setName(buf);
            }
        }
    }
    linExpr.end();
    pexModel->add(cnsts);
    //
    IloExpr objF(env);
    for (int w: og) {
        objF += p_w[w] * (*pex_th_w)[w];
    }
    pexCplex->getObjective().setExpr(objF);
}

void ColGenPE::getSol(Solution *sol) {
    for (int a: prob->A) {
        std::set<int> assignedTasks_a;
        for (int k: prob->K) {
            sol->y_ak[a][k] = pexCplex->getValue(pex_y_ak[a][k]);
            if (sol->y_ak[a][k] > 0.5) {
                assignedTasks_a.insert(k);
            }
        }
        for (int e: prob->E_a[a]) {
            int w = -1;
            bool chosenColExist = false;
            for (int i = 0; i < og_ae[a][e].size(); i++) {
                w = og_ae[a][e][i];
                if (pexCplex->getValue((*pex_th_w)[w]) > 0.5) {
                    chosenColExist = true;
                    break;
                }
            }
            if (chosenColExist) {
                for (int k: assignedTasks_a) {
                    if (og_tsk[w].find(k) == og_tsk[w].end()) {
                        sol->z_aek[a][e][k] = 1.0;
                    }
                }
                int n0 = og_rut[w][0], n1;
                sol->u_aei[a][e][n0] = og_arT[w][0];
                for (int i = 1; i < og_rut[w].size(); i++) {
                    n1 = og_rut[w][i];
                    sol->x_aeij[a][e][n0][n1] = 1.0;
                    //
                    n0 = n1;
                    sol->u_aei[a][e][n0] = og_arT[w][i];
                }
            } else {
                for (int k: assignedTasks_a) {
                    sol->z_aek[a][e][k] = 1.0;
                }
                int n0 = prob->R_ae[a][e][0], n1;
                sol->u_aei[a][e][n0] = prob->al_i[n0];
                for (int i = 1; i < prob->R_ae[a][e].size(); i++) {
                    n1 = prob->R_ae[a][e][i];
                    sol->x_aeij[a][e][n0][n1] = 1.0;
                    //
                    double erest_arrvTime = sol->u_aei[a][e][n0] + prob->t_ij[n0][n1];
                    double actual_arrvTime = erest_arrvTime > prob->al_i[n1] ? erest_arrvTime : prob->al_i[n1];
                    //
                    n0 = n1;
                    sol->u_aei[a][e][n0] = actual_arrvTime;
                }
            }
        }
    }
}

