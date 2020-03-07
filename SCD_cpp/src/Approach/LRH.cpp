//
//  LRH.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 25/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/Base.hpp"

Solution* LRH::solve() {
    baseCplex->solve();
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            for (int k: prob->K_ae[a][e]) {
                lrh_l_aek[a][e][k] = baseCplex->getDual((*COM_cnsts)[COM_cnsts_index[a][e][k]]);
            }
        }
    }
    logging("Initialization",
            " ", " ", " ", " ",
            " ");
    //
    while (tt->get_elapsedTimeCPU() < TIME_LIMIT_SEC) {
        solve_dualProblem();
        solve_pexModel();
        updateLMs();
        logging("IterSummary",
                " ", " ", " ", " ",
                "Dual Gap: " + std::to_string(bestSol->gap));
        if (bestSol->gap < DUAL_GAP_LIMIT) {
            break;
        }
        if (NUM_ITER_LIMIT == ++numIters) {
            break;
        }
    }
    return bestSol;
}

void LRH::updateLMs() {
    //
    // Update the step isze first!
    //
    if (L_V < L_V_star) {
        // Better solution
        L_V_star = L_V;
        noUpdateCounter_L = 0;
    } else {
        // No improvement
        noUpdateCounter_L += 1;
        if (noUpdateCounter_L == NO_IMPROVEMENT_LIMIT) {
            ut *= STEP_DECREASE_RATE;
            noUpdateCounter_L = 0;
        }
    }
    // Set at (Scale for the movement)
    double numerator = ut * (L_V_star - F_V_star);
    double denominator2 = 0.0;
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            for (int k: prob->K_ae[a][e]) {
                double sumX = 0.0;
                for (int j: prob->N_ae[a][e]) {
                    sumX += lrh_x_aeij[a][e][j][prob->n_k[k]];
                }
                denominator2 += pow(lrh_y_ak[a][k] - lrh_z_aek[a][e][k] - sumX, 2);
            }
        }
    }
    double at = numerator / denominator2;
    //
    // Update multipliers
    //
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            for (int k: prob->K_ae[a][e]) {
                double sumX = 0.0;
                for (int j: prob->N_ae[a][e]) {
                    sumX += lrh_x_aeij[a][e][j][prob->n_k[k]];
                }
                double lm = lrh_l_aek[a][e][k] + at * (lrh_y_ak[a][k] - lrh_z_aek[a][e][k] - sumX);
                lrh_l_aek[a][e][k] = lm < 0.0 ? 0.0 : lm;
            }
        }
    }
    logging("updateLMs",
            " ", " ", " ", " ",
            "at: " + std::to_string(at) +"  ut: " + std::to_string(at));
}

void LRH::solve_dualProblem() {
    L1_V = 0.0; L2_V = 0.0;
    //
    solve_etaModel();
    logging("solve_etaModel",
            " ", " ", " ", " ",
            "L1V: " + std::to_string(L1_V));
    solve_rutModels();
    logging("solve_rutModels",
            " ", " ", " ", " ",
            "L2V: " + std::to_string(L2_V));
    L_V = L1_V + L2_V;
    if (numIters == 0) {
        logging("solveDuals",
                "inf", std::to_string(L_V), " ", " ",
                " ");
    } else {
        logging("solve_dualProblem",
                std::to_string(L_V_star), std::to_string(L_V),
                " ", " ",
                " ");
    }
}

void LRH::build_etaModel() {
    eta_y_ak = Base::gen_y_ak(prob, env, 'I');
    eta_z_aek = Base::gen_z_aek(prob, env, 'I');
    //
    etaModel = new IloModel(env);
    //
    IloExpr objF = eta_y_ak[0][0];  // this function is a dummy; later the objective will be updated
    etaModel->add(IloMinimize(env, objF));
    objF.end();
    //
    Base::def_ETA_cnsts(prob, env, eta_y_ak, eta_z_aek, etaModel);
    //
    etaCplex = new IloCplex(*etaModel);
    etaCplex->setOut(env.getNullStream());
}

void LRH::update_etaModel() {
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
}

void LRH::solve_etaModel() {
    update_etaModel();
    etaCplex->solve();
    if (etaCplex->getStatus() == IloAlgorithm::Infeasible) {
        etaCplex->exportModel(lpPath.c_str());
        throw "the etaModel is infeasible";
    }
    if (etaCplex->getStatus() == IloAlgorithm::Optimal) {
        L1_V = -etaCplex->getObjValue();
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

void LRH::build_rutModels() {
    std::vector<std::string> cut_names {"hSE", "hCA", "hRS", "hIP"};
    IloCplex::CutManagement cutManagerType = IloCplex::UseCutForce;
    IloBool isLocalCutAdd = IloFalse;
    //
    rmm::RouteMM *rutMM;
    rmm::BnC *bncMM;
    for (int a: prob->A) {
        std::vector<rmm::RouteMM*> a_rutMMs;
        std::vector<rmm::BnC*> a_BnCs;
        std::vector<double**> a_x_ij;
        std::vector<double*> a_u_i;
        for (int e: prob->E_a[a]) {
            if (_ROUTER == "ILP") {
                rutMM = new rmm::ILP(prob->RP_ae[a][e], "", true);
                a_rutMMs.push_back(rutMM);
            } else if (_ROUTER == "BnC") {
                bncMM = new rmm::BnC(prob->RP_ae[a][e], "", cut_names, true, cutManagerType, isLocalCutAdd, tt);
                rutMM = bncMM;
                a_rutMMs.push_back(rutMM);
                a_BnCs.push_back(bncMM);
            } else {
                throw "a undefied router is requested";
            }
            rutMM->cplex->setWarning(rutMM->env.getNullStream());
            if (turnOnInitSol) {
                double **x_ij = new double*[prob->RP_ae[a][e]->N.size()];
                double *u_i = new double[prob->RP_ae[a][e]->N.size()];
                for (int i: prob->RP_ae[a][e]->N) {
                    x_ij[i] = new double[prob->RP_ae[a][e]->N.size()];
                    for (int j: prob->RP_ae[a][e]->N) {
                        x_ij[i][j] = 0.0;
                    }
                    u_i[i] = 0;
                }
                a_x_ij.push_back(x_ij);
                a_u_i.push_back(u_i);
            }
        }
        rutMMs.push_back(a_rutMMs);
        if (_ROUTER == "BnC") {
            rutBnCs.push_back(a_BnCs);
        }
        if (turnOnInitSol) {
            rut_x_ij.push_back(a_x_ij);
            rut_u_i.push_back(a_u_i);
        }
    }
}

void LRH::update_rutMM(int a, int e) {
    rmm::RouteMM *rutMM = rutMMs[a][e];
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
    //
    if (turnOnInitSol) {
        rutMM->set_initSol(rut_x_ij[a][e], rut_u_i[a][e]);
    }
}

void LRH::solve_rutModels() {
    rmm::RouteMM *rutMM;
    rmm::BnC *bncMM;
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            update_rutMM(a, e);
            if (_ROUTER == "BnC" && turOnCutPool) {
                bncMM = rutBnCs[a][e];
                bncMM->add_detectedCuts2MM();
            }
            rutMM = rutMMs[a][e];
            rutMM->cplex->solve();
            if (rutMM->cplex->getStatus() == IloAlgorithm::Infeasible) {
                rutMM->cplex->exportModel(lpPath.c_str());
                throw "the rutMM is infeasible";
            }
            if (rutMM->cplex->getStatus() == IloAlgorithm::Optimal) {
                L2_V += -rutMM->cplex->getObjValue();
                double x, u;
                for (int rutID0: prob->RP_ae[a][e]->N) {
                    int scdID0 = prob->rutID_scdID_ae[a][e][rutID0];
                    for (int rutID1: prob->RP_ae[a][e]->N) {
                        int scdID1 = prob->rutID_scdID_ae[a][e][rutID1];
                        x = rutMM->cplex->getValue(rutMM->x_ij[rutID0][rutID1]);
                        lrh_x_aeij[a][e][scdID0][scdID1] = x;
                        if (turnOnInitSol) {
                            rut_x_ij[a][e][rutID0][rutID1] = x;
                        }
                    }
                    u = rutMM->cplex->getValue(rutMM->u_i[rutID0]);
                    lrh_u_aei[a][e][scdID0] = u;
                    if (turnOnInitSol) {
                        rut_u_i[a][e][rutID0] = u;
                    }
                }
            } else {
                throw "the rutMM is not solved optimally";
            }
        }
    }
}

void LRH::build_pexModel() {
    pex_y_ak = Base::gen_y_ak(prob, env, 'I');
    pex_z_aek = Base::gen_z_aek(prob, env, 'I');
    //
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
    Base::def_ETA_cnsts(prob, env, pex_y_ak, pex_z_aek, pexModel);
    char buf[DEFAULT_BUFFER_SIZE];
    IloExpr linExpr(env);
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            for (int k: prob->K_ae[a][e]) {
                linExpr.clear();
                sprintf(buf, "CC(%d)(%d)(%d)", a, e, k);
                linExpr += pex_y_ak[a][k];
                double sumX = 0.0;
                for (int j: prob->N_ae[a][e]) {
                    sumX += lrh_x_aeij[a][e][j][prob->n_k[k]];
                }
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

void LRH::update_pexModel() {
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

void LRH::solve_pexModel() {
    update_pexModel();
    pexCplex->solve();
    if (pexCplex->getStatus() == IloAlgorithm::Infeasible) {
        pexCplex->exportModel(lpPath.c_str());
        throw "the pexCplex is infeasible";
    }
    if (pexCplex->getStatus() == IloAlgorithm::Optimal) {
        F_V = pexCplex->getObjValue();
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
            for (int a: prob->A) {
                for (int k: prob->K) {
                    bestSol->y_ak[a][k] = pexCplex->getValue(pex_y_ak[a][k]);
                }
                for (int e: prob->E_a[a]) {
                    for (int k: prob->K) {
                        bestSol->z_aek[a][e][k] = pexCplex->getValue(pex_z_aek[a][e][k]);
                    }
                    for (int i: prob->N_ae[a][e]) {
                        for (int j: prob->N_ae[a][e]) {
                            bestSol->x_aeij[a][e][i][j] = lrh_x_aeij[a][e][i][j];
                        }
                        bestSol->u_aei[a][e][i] = lrh_u_aei[a][e][i];
                    }
                }
            }
        }
    } else {
        throw "the pexCplex is not solved optimally";
    }
    logging("solve_pexModel",
            " ", " ",
            std::to_string(F_V_star), std::to_string(F_V),
            " ");
}

void LRH::logging(std::string indicator,
                  std::string _L_V_star, std::string _L_V, std::string _F_V_star, std::string _F_V,
                  std::string note) {
    if (logPath != "") {
        std::string _log(std::to_string(numIters) +
                    "," + std::to_string(tt->get_elapsedTimeWall()) +
                     "," + std::to_string(tt->get_elapsedTimeCPU()));
        std::vector<std::string> tempVec {indicator, _L_V_star, _L_V, _F_V_star, _F_V, note};
        for (std::string str: tempVec) {
            _log += "," + str;
        }
        char log[_log.size() + 1];
        std::strcpy(log, _log.c_str());
        appendRow(logPath, log);
    }
}
