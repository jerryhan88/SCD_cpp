//
//  LRH.hpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 14/6/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef LRH_h
#define LRH_h

#include <vector>

#include <ilcplex/ilocplex.h>

#include "SolApprBase.hpp"
#include "RouterSCD.hpp"

#include "Problem.hpp"
#include "Solution.hpp"


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

class PrimalExtractor {
public:
    Problem *prob;
    //
    IloEnv env;
    IloModel *pexModel;
    IloCplex *pexCplex;
    //
    double ****lrh_x_aeij, ***lrh_u_aei;
    //
    IloNumVar **pex_y_ak, ***pex_z_aek;
    IloRangeArray *pex_COM_cnsts;
    long ***pex_COM_cnsts_index;
    //
    PrimalExtractor(Problem *prob, double ****lrh_x_aeij, double ***lrh_u_aei) {
        this->prob = prob;
        this->lrh_x_aeij = lrh_x_aeij;
        this->lrh_u_aei = lrh_u_aei;
        //
        pex_y_ak = new_inv_ak(prob, env, 'I');
        pex_z_aek = new_inv_aek(prob, env, 'I');
        //
        pex_COM_cnsts = new IloRangeArray(env);
        pex_COM_cnsts_index = new long**[prob->A.size()];
        for (int a : prob->A) {
            pex_COM_cnsts_index[a] = new long*[prob->E_a[a].size()];
            for (int e: prob->E_a[a]) {
                pex_COM_cnsts_index[a][e] = new long[prob->K.size()];
            }
        }
        build();
        
    }
    ~PrimalExtractor() {
        delete_inv_ak(prob, pex_y_ak);
        delete_inv_aek(prob, pex_z_aek);
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                delete [] pex_COM_cnsts_index[a][e];
            }
            delete [] pex_COM_cnsts_index[a];
        }
        delete [] pex_COM_cnsts_index;
        //
        delete pexCplex; delete pexModel;
        delete pex_COM_cnsts;
        env.end();
    }
    void build();
    void update();
    void getSol(Solution *sol);
};

class LRH : public BaseMM {
public:
    std::string _ROUTER, _EXTRACTOR;
    double STEP_DECREASE_RATE, DUAL_GAP_LIMIT;
    unsigned int NO_IMPROVEMENT_LIMIT, NUM_ITER_LIMIT;
    //
    double L1_V, L2_V, L_V, L_V_star, F_V, F_V_star;
    double ut;
    int numIters, noUpdateCounter_L;
    double **lrh_y_ak, ***lrh_z_aek, ****lrh_x_aeij, ***lrh_u_aei, ***lrh_l_aek;
    //
    Allocator *alr;
    std::vector<std::vector<RouterSCD*>> routers;
    PrimalExtractor *pex;
    //
    Solution *bestSol;
    std::vector<iiidTup> subProbSolwallTs;
    std::vector<double> dualSols, primalSols;
    //
    LRH(Problem *prob, TimeTracker *tt,
        unsigned long time_limit_sec, int numThreads,
        std::string logPath, std::string lpPath, std::string lp_algo,
        std::string _router,
        double dual_gap_limit, unsigned int num_iter_limit, unsigned int no_improvement_limit) : BaseMM(prob, tt, time_limit_sec, numThreads, logPath, lpPath, lp_algo, 'C') {
        _ROUTER = _router;
        DUAL_GAP_LIMIT = dual_gap_limit;
        NUM_ITER_LIMIT = num_iter_limit;
        NO_IMPROVEMENT_LIMIT = no_improvement_limit;
        L_V_star = DBL_MAX; F_V_star = -DBL_MAX;
        STEP_DECREASE_RATE = 0.5;
        ut = 2.0;
        noUpdateCounter_L = 0;
        numIters = 0;
        //
        lrh_y_ak = new_dbl_ak(prob);
        lrh_z_aek = new_dbl_aek(prob);
        lrh_x_aeij = new_dbl_aeij(prob);
        lrh_u_aei = new_dbl_aei(prob);
        lrh_l_aek = new_dbl_aek(prob);
        //
        build_allocator();
        build_rutModels();
        build_extractor();
        //
        bestSol = nullptr;
        //
        if (logPath != "") {
            std::string _header("Iteration,wallT,cpuT,Function,Note");
            char header[_header.size() + 1];
            std::strcpy(header, _header.c_str());
            createCSV(logPath, header);
        }
    }
    ~LRH() {
        delete_dbl_ak(prob, lrh_y_ak);
        delete_dbl_aek(prob, lrh_z_aek);
        delete_dbl_aeij(prob, lrh_x_aeij);
        delete_dbl_aei(prob, lrh_u_aei);
        delete_dbl_aek(prob, lrh_l_aek);
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                delete routers[a][e];
            }
        }
        //
        delete alr;
        delete pex;
        env.end();
        routers.clear();
    }
    //
    Solution* solve();
    //
    void init_LMs();
    void solve_dualProblem();
    bool updateLMs();

private:
    void build_allocator();
    void solve_etaModel();
    //
    void build_rutModels();
    void solve_rutModels();
    //
    void build_extractor();
    void solve_pexModel();
    //
    void logging(std::string indicator, std::string note);
};


#endif /* LRH_h */
