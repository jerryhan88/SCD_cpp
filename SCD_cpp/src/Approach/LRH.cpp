//
//  LRH.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 25/2/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/LRH.hpp"


std::mutex lrh_mtx_vec;

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

Solution* LRH::solve() {
    std::cout << "\t" << "The LP model has been build " << TimeTracker::get_curTime();
    logging("Begin_Initialization", " ");
    init_LMs();
    logging("End_Initialization", " ");
    std::cout << "\t" << "The LP model has been solved " << TimeTracker::get_curTime();
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
            logging("IterSummary", "Dual Gap:" + std::to_string(bestSol->gap));
        }
        std::cout << "\t" << numIters << "th iter finished; " << TimeTracker::get_curTime();
        if (NUM_ITER_LIMIT == ++numIters) {
            terminationCon = 0;
            break;
        }
//        if (bestSol->gap < DUAL_GAP_LIMIT) {
//            terminationCon = 1;
//            break;
//        }
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

void LRH::init_LMs() {
    logging("Begin_baseSolve", " ");
    unsigned long TiLim = LP_TIME_LIMIT < time_limit_sec ? LP_TIME_LIMIT : time_limit_sec;
    baseCplex->setParam(IloCplex::TiLim, TiLim);
    baseCplex->solve();
    logging("End_baseSolve", " ");
    //
    logging("Begin_baseGetSol", " ");
    for (int a : prob->A) {
        for (int e: prob->E_a[a]) {
            for (int k: prob->K_ae[a][e]) {
                lrh_l_aek[a][e][k] = baseCplex->getDual((*COM_cnsts)[COM_cnsts_index[a][e][k]]);
            }
        }
    }
    logging("End_baseGetSol", " ");
}

bool LRH::updateLMs() {
    logging("Begin_updateLMs", " ");
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
                }
            }
        }
    } else {
        at = 0.0;
    }
    logging("End_updateLMs", "at: " + std::to_string(at) +"  ut: " + std::to_string(ut));
    if (at == 0.0) {
        return true;
    } else {
        return false;
    }
}

void LRH::solve_dualProblem() {
    logging("Begin_solve_dualProblem", " ");
    L1_V = 0.0; L2_V = 0.0;
    //
    solve_etaModel();
    solve_rutModels();
    L_V = L1_V + L2_V;
    if (numIters == 0) {
        logging("End_solve_dualProblem", "L_V*:inf;L_V:" + std::to_string(L_V));
    } else {
        logging("End_solve_dualProblem", "L_V*:" + std::to_string(L_V_star) + ";L_V:" + std::to_string(L_V));
    }
}

void LRH::build_allocator() {
    alr = new Allocator(prob, tt, lrh_l_aek);
}

void LRH::solve_etaModel() {
    logging("Begin_solve_etaModel", " ");
    //
    logging("Begin_solve_etaModel_updateModel", " ");
    alr->update();
    
    logging("End_solve_etaModel_updateModel", " ");
    //
    logging("Begin_solve_etaModel_solve", " ");
    
    unsigned long timeLimit4AssSolving = time_limit_sec;
    timeLimit4AssSolving -= (unsigned long) tt->get_elapsedTimeWall();
    alr->etaCplex->setParam(IloCplex::TiLim, timeLimit4AssSolving);
    try {
        alr->etaCplex->solve();
    } catch (IloException& e) {
       std::cerr << "Concert exception caught: " << e << std::endl;
    }

    if (alr->etaCplex->getStatus() == IloAlgorithm::Infeasible) {
        alr->etaCplex->exportModel(lpPath.c_str());
        throw "the etaModel is infeasible";
    }
    logging("End_solve_etaModel_solve", " ");
    //
    logging("Begin_solve_etaModel_getSol", " ");
    alr->getSol(&L1_V, lrh_y_ak, lrh_z_aek);
    logging("End_solve_etaModel_getSol", " ");
    //
    logging("End_solve_etaModel", "L1V: " + std::to_string(L1_V));
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
//                rutLogPath += logPath.substr(0, logPath.find_last_of("."));
//                rutLogPath += "_" + std::to_string(a) + std::to_string(e);
//                rutLogPath += ".log";
//                //
//                rutLpLogPath += logPath.substr(0, logPath.find_last_of("."));
//                rutLpLogPath += "_" + std::to_string(a) + std::to_string(e);
//                rutLpLogPath += ".lp";
            }
            if (_ROUTER == "ILP") {
                rut = new RouterILP(prob, tt, rut_time_limit_sec,
                                    a, e, lrh_l_aek,
                                    rutLogPath, rutLpLogPath);
            } else if (_ROUTER.rfind("BnC", 0) == 0) {
                std::string tmp("BnC");
                std::string option = _ROUTER.substr(tmp.size(), _ROUTER.size());
                bool turOnCutPool = option == "oc";
                if (!turOnCutPool) {
                    rut = new RouterBnC(prob, tt, rut_time_limit_sec,
                    a, e, lrh_l_aek,
                    rutLogPath, rutLpLogPath);
                } else {
                    rut = new RouterBnCoc(prob, tt, rut_time_limit_sec,
                    a, e, lrh_l_aek,
                    rutLogPath, rutLpLogPath);
                }
                
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
    logging("Begin_solve_rutModels", " ");
    auto rut_update = [](RouterSCD* rut) {
        rut->update();
//        char fpath[2048];
//        sprintf(fpath, "%s_a%02d-e%02d-ni%02d.lp", pathPrefix.c_str(), rut->a, rut->e, numIters);
//        rut->rutCplex->exportModel(fpath);
    };
    auto rut_solve = [](RouterSCD* rut, int numIters, std::vector<iiidTup> *vec, std::string &pathPrefix) {
        double start = rut->tt->get_elapsedTimeWall();
        rut->solve();
        double elapsedTime = rut->tt->get_elapsedTimeWall() - start;
        std::unique_lock<std::mutex> lock(lrh_mtx_vec);
        vec->push_back(std::make_tuple(numIters, rut->a, rut->e, elapsedTime));
    };
    auto rut_getSol = [](RouterSCD* rut) {
        rut->getSol();
    };
    //
    if (routerPool.getThreadCount() > 1) {
        std::vector<std::shared_future<void>> threadJobs;
        logging("Begin_solve_rutModels_update", " ");
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                std::shared_future<void> returnValue = (std::shared_future<void>) routerPool.push(rut_update, routers[a][e]);
                threadJobs.push_back(returnValue);
            }
        }
        std::for_each(threadJobs.begin(), threadJobs.end(), [](std::shared_future<void> x){ x.get(); });
        logging("End_solve_rutModels_update", " ");
        threadJobs.clear();
        logging("Begin_solve_rutModels_solve", " ");
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                std::shared_future<void> returnValue = (std::shared_future<void>) routerPool.push(rut_solve, routers[a][e], numIters, &subProbSolwallTs, pathPrefix);
                threadJobs.push_back(returnValue);
            }
        }
        std::for_each(threadJobs.begin(), threadJobs.end(), [](std::shared_future<void> x){ x.get(); });
        logging("End_solve_rutModels_solve", " ");
        threadJobs.clear();
        logging("Begin_solve_rutModels_getSol", " ");
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                std::shared_future<void> returnValue = (std::shared_future<void>) routerPool.push(rut_getSol, routers[a][e]);
                threadJobs.push_back(returnValue);
            }
        }
        std::for_each(threadJobs.begin(), threadJobs.end(), [](std::shared_future<void> x){ x.get(); });
        logging("End_solve_rutModels_getSol", " ");
    } else {
        logging("Begin_solve_rutModels_update", " ");
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                rut_update(routers[a][e]);
            }
        }
        logging("End_solve_rutModels_update", " ");
        logging("Begin_solve_rutModels_solve", " ");
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                rut_solve(routers[a][e], numIters, &subProbSolwallTs, pathPrefix);
            }
        }
        logging("End_solve_rutModels_solve", " ");
        logging("Begin_solve_rutModels_getSol", " ");
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                rut_getSol(routers[a][e]);
            }
        }
        logging("End_solve_rutModels_getSol", " ");
    }
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
    logging("End_solve_rutModels", "L2V: " + std::to_string(L2_V));
}

void LRH::build_extractor() {
    pex = new PrimalExtractor(prob, lrh_x_aeij, lrh_u_aei);
}

void LRH::solve_pexModel() {
    logging("Begin_solve_pexModel", " ");
    logging("Begin_solve_pexModel_update", " ");
    pex->update();
    logging("End_solve_pexModel_update", " ");
    logging("Begin_solve_pexModel_solve", " ");
    pex->pexCplex->solve();
    if (pex->pexCplex->getStatus() == IloAlgorithm::Infeasible) {
        throw "the pexCplex is infeasible";
    }
    assert(pex->pexCplex->getStatus() == IloAlgorithm::Optimal);
    logging("End_solve_pexModel_solve", " ");
    F_V = pex->pexCplex->getObjValue();
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
        logging("Begin_solve_pexModel_getSol", " ");
        pex->getSol(bestSol);
        logging("End_solve_pexModel_getSol", " ");
        bestSol->lastUpdatedIter = numIters;
    }
    logging("End_solve_pexModel", "F_V*:" + std::to_string(F_V_star) + ";F_V:" + std::to_string(F_V));
}

void LRH::logging(std::string indicator, std::string note) {
    if (logPath != "") {
        std::string _log(std::to_string(numIters) +
                    "," + std::to_string(tt->get_elapsedTimeWall()) +
                     "," + std::to_string(tt->get_elapsedTimeCPU()));
        _log += "," + indicator + "," + note;
        char log[_log.size() + 1];
        std::strcpy(log, _log.c_str());
        appendRow(logPath, log);
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


void PrimalExtractor::build() {
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

void PrimalExtractor::update() {
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

void PrimalExtractor::getSol(Solution *sol) {
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
