//
//  LRH.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 25/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/Solver.hpp"


std::mutex lrh_mtx_vec;

Solution* LRH::solve() {
    baseCplex->solve();
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            for (int k: prob->K_ae[a][e]) {
                lrh_l_aek[a][e][k] = baseCplex->getDual((*COM_cnsts)[COM_cnsts_index[a][e][k]]);
            }
        }
    }
    if (logPath != "") {
        logging("Initialization",
        " ", " ", " ", " ",
        " ");
    }
    std::cout << "\t" << "The LP model has been solved " << tt->get_curTime();
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
        solve_pexModel();
        primalSols.push_back(F_V);
        bool zeroDenominator = updateLMs();
        if (logPath != "") {
            logging("IterSummary",
            " ", " ", " ", " ",
            "Dual Gap: " + std::to_string(bestSol->gap));
        }
        std::cout << "\t" << numIters << "th iter finished; " << tt->get_curTime();
        if (NUM_ITER_LIMIT == ++numIters) {
            terminationCon = 0;
            break;
        }
        if (bestSol->gap < DUAL_GAP_LIMIT) {
            terminationCon = 1;
            break;
        }
        if (zeroDenominator) {
            terminationCon = 2;
            break;
        }
    }
    //
    char note[2048];
    sprintf(note,
            "\"{\'lastUpdatedIter\': %d, \'bestSolCpuT\': %f, \'bestSolWallT\': %f, \'terminationCon\': %d, \'numIters\': %d}\"",
            bestSol->lastUpdatedIter, bestSol->cpuT, bestSol->wallT, terminationCon, numIters);
    bestSol->note = std::string(note);
    bestSol->cpuT = tt->get_elapsedTimeCPU();
    bestSol->wallT = tt->get_elapsedTimeWall();
    //
    char row[2048];
    std::fstream fout_csv;
    //
    std::string solQual_csv = pathPrefix + "_solQual.csv";
    fout_csv.open(solQual_csv, std::ios::out);
    fout_csv << "numIter,dualSol,primalSol" << "\n";
    for (int i = 0; i < numIters; i++) {
        sprintf(row, "%d,%f,%f", i, dualSols[i], primalSols[i]);
        fout_csv << row << "\n";
    }
    fout_csv.close();
    //
    std::string subProbSolWallT_csv = pathPrefix + "_subProbSolWallT.csv";
    fout_csv.open(subProbSolWallT_csv, std::ios::out);
    fout_csv << "numIter,a,e,elapsedWallT" << "\n";
    int inum, a, e;
    double elapsedWallT;
    for (iiidTup ele: subProbSolwallTs) {
        std::tie(inum, a, e, elapsedWallT) = ele;
        sprintf(row, "%d,%d,%d,%f", inum, a, e, elapsedWallT);
        fout_csv << row << "\n";
    }
    fout_csv.close();
    
    return bestSol;
}

bool LRH::updateLMs() {
    //
    // Update the step isze first!
    //
    if (L_V < L_V_star) {
        // Better solution
        L_V_star = L_V;
        noUpdateCounter_L = 0;
        //
        bestSol->gap = (L_V_star - F_V_star) / F_V_star;
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
    //
    // Update multipliers
    //
    
    
//    char note[2048];
//    std::cout << "\n";
    double at;
    if (denominator2 != 0) {
        at = numerator / denominator2;
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                for (int k: prob->K_ae[a][e]) {
                    double sumX = 0.0;
                    for (int j: prob->N_ae[a][e]) {
                        sumX += lrh_x_aeij[a][e][j][prob->n_k[k]];
                    }
                    double lm = lrh_l_aek[a][e][k] + at * (lrh_y_ak[a][k] - lrh_z_aek[a][e][k] - sumX);
                    lrh_l_aek[a][e][k] = lm < 0.0 ? 0.0 : lm;
//                    sprintf(note, "%d; (%d,%d,%d): %e", numIters, a, e, k, lrh_l_aek[a][e][k]);
//                    std::cout << note << "\n";
                }
            }
        }
    } else {
        at = 0.0;
    }
    if (logPath != "") {
        logging("updateLMs",
        " ", " ", " ", " ",
        "at: " + std::to_string(at) +"  ut: " + std::to_string(ut));
    }
    if (at == 0.0) {
        return true;
    } else {
        return false;
    }
}

void LRH::solve_dualProblem() {
    L1_V = 0.0; L2_V = 0.0;
    //
    solve_etaModel();
    if (logPath != "") {
        logging("solve_etaModel",
        " ", " ", " ", " ",
        "L1V: " + std::to_string(L1_V));
    }
    solve_rutModels();
    if (logPath != "") {
        logging("solve_rutModels",
        " ", " ", " ", " ",
        "L2V: " + std::to_string(L2_V));
    }
    L_V = L1_V + L2_V;
    if (logPath != "") {
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
}

void LRH::build_etaModel() {
    eta_y_ak = gen_y_ak(prob, env, 'I');
    eta_z_aek = gen_z_aek(prob, env, 'I');
    //
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
    RouterSCD *rut;
    unsigned long rut_time_limit_sec = time_limit_sec;
    for (int a: prob->A) {
        std::vector<RouterSCD*> a_routers;
        std::vector<double**> a_x_ij;
        std::vector<double*> a_u_i;
        for (int e: prob->E_a[a]) {
            std::string rutLogPath = "";
            std::string rutLpLogPath = "";
            if (logPath != "") {
                rutLogPath += logPath.substr(0, logPath.find_last_of("."));
                rutLogPath += "_" + std::to_string(a) + std::to_string(e);
                rutLogPath += ".log";
                //
                rutLpLogPath += logPath.substr(0, logPath.find_last_of("."));
                rutLpLogPath += "_" + std::to_string(a) + std::to_string(e);
                rutLpLogPath += ".lp";
            }
            if (_ROUTER == "ILP") {
                rut = new RouterILP(prob, tt, rut_time_limit_sec,
                                    a, e, lrh_l_aek,
                                    rutLogPath, rutLpLogPath);
            } else if (_ROUTER.rfind("BnC", 0) == 0) {
                std::string tmp("BnC");
                std::string option = _ROUTER.substr(tmp.size(), _ROUTER.size());
                bool turOnCutPool = option == "oc";
                rut = new RouterBnC(prob, tt, rut_time_limit_sec,
                                    a, e, lrh_l_aek,
                                    rutLogPath, rutLpLogPath, turOnCutPool);
            } else if (_ROUTER == "GH") {
                rut = new RouterGH(prob, tt, rut_time_limit_sec,
                                   a, e, lrh_l_aek,
                                   rutLogPath);
            } else {
                throw "a undefied router is requested";
            }
            a_routers.push_back(rut);;
        }
        routers.push_back(a_routers);
    }
}

void LRH::solve_rutModels() {
    std::vector<std::shared_future<void>> threadJobs;

    
    auto solve_a_rutModel = [](int numIters, RouterSCD* rut, std::vector<iiidTup> *vec, std::string &pathPrefix) {
        double start = rut->tt->get_elapsedTimeWall();
        rut->update();
//        if (rut->a == 0 && rut->e == 0) {
//            std::cout << 12;
//        }
//        char fpath[2048];
//        sprintf(fpath, "%s_a%02d-e%02d-ni%02d.lp", pathPrefix.c_str(), rut->a, rut->e, numIters);
//        rut->rutCplex->exportModel(fpath);
        
        rut->run();
        double elapsedTime = rut->tt->get_elapsedTimeWall() - start;
        std::unique_lock<std::mutex> lock(lrh_mtx_vec);
        vec->push_back(std::make_tuple(numIters, rut->a, rut->e, elapsedTime));
    };
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            std::shared_future<void> returnValue = (std::shared_future<void>) pool.push(solve_a_rutModel, numIters, routers[a][e], &subProbSolwallTs, pathPrefix);
            threadJobs.push_back(returnValue);
        }
    }
    std::for_each(threadJobs.begin(), threadJobs.end(), [](std::shared_future<void> x){
        x.get();});
    for (int a: prob->A) {
        for (int e: prob->E_a[a]) {
            L2_V += -routers[a][e]->rut_objV;
            for (int rutID0: prob->RP_ae[a][e]->N) {
                int scdID0 = prob->rutID_scdID_ae[a][e][rutID0];
                for (int rutID1: prob->RP_ae[a][e]->N) {
                    int scdID1 = prob->rutID_scdID_ae[a][e][rutID1];
                    lrh_x_aeij[a][e][scdID0][scdID1] = routers[a][e]->rut_x_ij[rutID0][rutID1];
                }
                lrh_u_aei[a][e][scdID0] = routers[a][e]->rut_u_i[rutID0];
            }
        }
    }
}

void LRH::solve_pexModel() {
    pex->solve(lrh_x_aeij, lrh_u_aei);
    F_V = pex->F_V;
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
                bestSol->y_ak[a][k] = pex->y_ak[a][k];
            }
            for (int e: prob->E_a[a]) {
                for (int k: prob->K) {
                    bestSol->z_aek[a][e][k] = pex->z_aek[a][e][k];
                }
                for (int i: prob->N_ae[a][e]) {
                    for (int j: prob->N_ae[a][e]) {
                        bestSol->x_aeij[a][e][i][j] = pex->x_aeij[a][e][i][j];
                    }
                    bestSol->u_aei[a][e][i] = pex->u_aei[a][e][i];
                }
            }
        }
        bestSol->lastUpdatedIter = numIters;
    }
    if (logPath != "") {
        logging("solve_pexModel",
        " ", " ",
        std::to_string(F_V_star), std::to_string(F_V),
        " ");
    }
}

void LRH::logging(std::string indicator,
                  std::string _L_V_star, std::string _L_V, std::string _F_V_star, std::string _F_V,
                  std::string note) {
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
