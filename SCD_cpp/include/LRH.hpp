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

class LRH;

class Allocator {
public:
    Problem *prob;
    //
    IloEnv env;
    IloModel *etaModel; // Evaluation on Task Assingment
    IloCplex *etaCplex;
    IloNumVar **eta_y_ary, ***eta_z_ary;
    //
    IloNum minGap, lastUpdatedST;
    IloBool aborted;
    //
    LRH *lrh;
    //
    Allocator(Problem *prob, TimeTracker *tt, LRH *lrh);
    
    
    ~Allocator() {
        del_ak_inv(prob, eta_y_ary);
        del_aek_inv(prob, eta_z_ary);
        //
        delete etaCplex; delete etaModel;
        env.end();
    }
    void build();
    void update();
    void getSol(double *L1_V, double **lrh_y_ary, double ***lrh_z_ary);
};

class PrimalExtractor {
public:
    Problem *prob;
    //
    IloEnv env;
    IloModel *pexModel;
    IloCplex *pexCplex;
    //
    LRH *lrh;
    //
    IloNumVar **pex_y_ary, ***pex_z_ary;
    IloRangeArray *pex_COM_cnsts;
    long ***pex_COM_cnsts_index;
    //
    PrimalExtractor(Problem *prob, LRH *lrh);
    ~PrimalExtractor() {
        del_ak_inv(prob, pex_y_ary);
        del_aek_inv(prob, pex_z_ary);
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
    void def_ETA_cnsts();
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
    double **lrh_y_ary, ***lrh_z_ary, ****lrh_x_ary, ***lrh_u_ary, ***lrh_l_ary;
    std::map<iiTup, double> lrh_y_map;
    std::map<iiiTup, double> lrh_z_map;
    std::map<iiiiTup, double> lrh_x_map;
    std::map<iiiTup, double> lrh_u_map;
    std::map<iiiTup, double> lrh_l_map;
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
        gen_y_dbl();
        gen_z_dbl();
        gen_x_dbl();
        gen_u_dbl();
        gen_l_dbl();
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
        del_y_dbl();
        del_z_dbl();
        del_x_dbl();
        del_u_dbl();
        del_l_dbl();
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
    //
    void set_y_dbl(int a, int k, double v);
    void set_z_dbl(int a, int e, int k, double v);
    //
    double get_y_dbl(int a, int k);
    double get_z_dbl(int a, int e, int k);
    double get_x_dbl(int a, int e, int i, int j);
    double get_u_dbl(int a, int e, int i);
    double get_l_dbl(int a, int e, int k);
private:
    void gen_y_dbl();
    void del_y_dbl();
    //
    void gen_z_dbl();
    void del_z_dbl();
    //
    void gen_x_dbl();
    void del_x_dbl();
    void set_x_dbl(int a, int e, int i, int j, double v);
    //
    void gen_u_dbl();
    void del_u_dbl();
    void set_u_dbl(int a, int e, int i, double v);
    //
    void gen_l_dbl();
    void del_l_dbl();
    void set_l_dbl(int a, int e, int k, double v);
    //
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
