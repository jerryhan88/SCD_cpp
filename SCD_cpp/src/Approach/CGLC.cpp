//
//  CGLC.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 13/6/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/CGLC.hpp"


Solution* CGLC::solve() {
    init_LMs();
    //
    int terminationCon = -1;
    while (tt->get_elapsedTimeWall() < time_limit_sec) {
        try {
            solve_dualProblem();
        } catch (const char * str) {
            std::cout << "Exception: " << str << std::endl;
            break;
        }
        dualSols.push_back(L_V);
        //
        solve_rmModel();
        primalSols.push_back(F_V);
        bool zeroDenominator = updateLMs();
        std::cout << "\t" << numIters << "th iter finished; " << TimeTracker::get_curTime();
        if (NUM_ITER_LIMIT == ++numIters) {
            terminationCon = 0;
            break;
        }
        if (zeroDenominator) {
            terminationCon = 2;
            break;
        }
    }
    
    char note[2048];
    sprintf(note,
            "\"{\'lastUpdatedIter\': %d, \'bestSolCpuT\': %f, \'bestSolWallT\': %f, \'terminationCon\': %d, \'numIters\': %d}\"",
            bestSol->lastUpdatedIter, bestSol->cpuT, bestSol->wallT, terminationCon, numIters);
    bestSol->note = std::string(note);
    bestSol->cpuT = tt->get_elapsedTimeCPU();
    bestSol->wallT = tt->get_elapsedTimeWall();
    
    return bestSol;
}


void CGLC::rm_build() {
    rmModel = new IloModel(env);
    //
    init_cols();
    //
    IloExpr objF(env);
    rmModel->add(IloMaximize(env, objF));
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
            linExpr += rm_y_ak[a][k];
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
                    linExpr += rm_y_ak[a][k];
                    cnsts.add(linExpr == 0);
                    cnsts[cnsts.getSize() - 1].setName(buf);
                }
            }
        }
    }
    //
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            int col_index = og_ae[a][e][0];
            sprintf(buf, "ONE(%d)(%d)", a, e);
            rm_ONE_cnsts->add((*rm_th_w)[col_index] == 1);
            (*rm_ONE_cnsts)[rm_ONE_cnsts->getSize() - 1].setName(buf);
            rm_ONE_cnsts_index[a][e] = rm_ONE_cnsts->getSize() - 1;
        }
    }
    linExpr.end();
    rmModel->add(cnsts);
    rmModel->add(*rm_ONE_cnsts);
    rmCplex = new IloCplex(*rmModel);
    rmCplex->setOut(env.getNullStream());
    rmCplex->setWarning(env.getNullStream());
}

void CGLC::init_cols() {
    char buf[2048];
    IloRangeArray cnsts(env);
    IloExpr linExpr(env);
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
            std::vector<int> route(prob->R_ae[a][e]);
            std::set<int> completedTasks;
            og_rut.push_back(route);
            std::vector<double> arT;
            int n0 = route[0];
            arT.push_back(prob->al_i[n0]);
            for (int i = 1; i < route.size(); i++) {
                int n1 = route[i];
                //
                double erest_arrvTime = arT[i - 1] + prob->t_ij[n0][n1];
                double actual_arrvTime = erest_arrvTime > prob->al_i[n1] ? erest_arrvTime : prob->al_i[n1];
                //
                n0 = n1;
                arT.push_back(actual_arrvTime);
            }
            og_arT.push_back(arT);
            og_tsk.push_back(completedTasks);
            //
            double expProf = 0.0;
            std::vector<int> indicator(prob->K.size(), 0);
            p_w.push_back(expProf);
            e_wk.push_back(indicator);
            //
            rm_th_w->add(IloNumVar(env, 0.0, 1.0, ILOINT));
            sprintf(buf, "o(%d)", col_index);
            (*rm_th_w)[col_index].setName(buf);
        }
    }
}

void CGLC::solve_rmModel() {
    rm_update();
    rmCplex->solve();
    if (rmCplex->getStatus() == IloAlgorithm::Infeasible) {
        throw "the rmModel is infeasible";
    }
    assert(rmCplex->getStatus() == IloAlgorithm::Optimal);
    F_V = rmCplex->getObjValue();
    if (F_V > F_V_star) {
        F_V_star = F_V;
        if (bestSol != nullptr) {
            delete bestSol;
        }
        bestSol = new Solution(prob);
        bestSol->objV = F_V_star;
        bestSol->gap = ((L_V > L_V_star ? L_V_star : L_V) - F_V_star) / F_V_star;
        bestSol->cpuT = tt->get_elapsedTimeCPU();
        bestSol->wallT = tt->get_elapsedTimeWall();
        rm_getSol(bestSol);
        bestSol->lastUpdatedIter = numIters;
    }
}


void CGLC::rm_update() {
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
            rm_th_w->add(IloNumVar(env, 0.0, 1.0, ILOINT));
            sprintf(buf, "o(%d)", col_index);
            (*rm_th_w)[col_index].setName(buf);
            long cnsts_id = rm_ONE_cnsts_index[a][e];
            (*rm_ONE_cnsts)[cnsts_id].setLinearCoef((*rm_th_w)[col_index], 1.0);
        }
        for (int w: newlyAddedW) {
            for (int k: completedTasks_a) {
                if (e_wk[w][k] == 0) {
                    continue;
                }
                linExpr.clear();
                sprintf(buf, "BI(%d)(%d)(%d)", a, w, k);
                linExpr += e_wk[w][k] * (*rm_th_w)[w];
                linExpr -= rm_y_ak[a][k];
                cnsts.add(linExpr <= 0);
                cnsts[cnsts.getSize() - 1].setName(buf);
            }
        }
    }
    linExpr.end();
    rmModel->add(cnsts);
    //
    IloExpr objF(env);
    for (int w: og) {
        objF += p_w[w] * (*rm_th_w)[w];
    }
    rmCplex->getObjective().setExpr(objF);
}


void CGLC::rm_getSol(Solution *sol) {
    for (int a: prob->A) {
        std::set<int> assignedTasks_a;
        for (int k: prob->K) {
            sol->y_ak[a][k] = rmCplex->getValue(rm_y_ak[a][k]);
            if (sol->y_ak[a][k] > 0.5) {
                assignedTasks_a.insert(k);
            }
        }
        for (int e: prob->E_a[a]) {
            int w = -1;
            bool chosenColExist = false;
            for (int i = 0; i < og_ae[a][e].size(); i++) {
                w = og_ae[a][e][i];
                if (rmCplex->getValue((*rm_th_w)[w]) > 0.5) {
                    chosenColExist = true;
                    break;
                }
            }
            assert(chosenColExist);
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
            
//            if (chosenColExist) {
//                for (int k: assignedTasks_a) {
//                    if (og_tsk[w].find(k) == og_tsk[w].end()) {
//                        sol->z_aek[a][e][k] = 1.0;
//                    }
//                }
//                int n0 = og_rut[w][0], n1;
//                sol->u_aei[a][e][n0] = og_arT[w][0];
//                for (int i = 1; i < og_rut[w].size(); i++) {
//                    n1 = og_rut[w][i];
//                    sol->x_aeij[a][e][n0][n1] = 1.0;
//                    //
//                    n0 = n1;
//                    sol->u_aei[a][e][n0] = og_arT[w][i];
//                }
//            } else {
//                for (int k: assignedTasks_a) {
//                    sol->z_aek[a][e][k] = 1.0;
//                }
//                int n0 = prob->R_ae[a][e][0], n1;
//                sol->u_aei[a][e][n0] = prob->al_i[n0];
//                for (int i = 1; i < prob->R_ae[a][e].size(); i++) {
//                    n1 = prob->R_ae[a][e][i];
//                    sol->x_aeij[a][e][n0][n1] = 1.0;
//                    //
//                    double erest_arrvTime = sol->u_aei[a][e][n0] + prob->t_ij[n0][n1];
//                    double actual_arrvTime = erest_arrvTime > prob->al_i[n1] ? erest_arrvTime : prob->al_i[n1];
//                    //
//                    n0 = n1;
//                    sol->u_aei[a][e][n0] = actual_arrvTime;
//                }
//            }
        }
    }
}
