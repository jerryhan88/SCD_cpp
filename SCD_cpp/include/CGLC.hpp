//
//  CGLC.hpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 14/6/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef CGLC_h
#define CGLC_h

#include <vector>

#include <ilcplex/ilocplex.h>

#include "LRH.hpp"

#include "Problem.hpp"
#include "Solution.hpp"


class CGLC : public LRH {
public:
    IloModel *rmModel;
    IloCplex *rmCplex;
    //
    IloNumVar **rm_y_ary;
    IloNumVarArray *rm_th_w;
    IloRangeArray *rm_ONE_cnsts;
    long **rm_ONE_cnsts_index;
    //
    std::vector<int> og;
    std::vector<std::vector<int>> og_a;
    std::vector<std::vector<std::vector<int>>> og_ae;
    std::vector<double> p_w;
    std::vector<std::vector<int>> e_wk;
    //
    std::vector<std::vector<int>> og_rut;
    std::vector<std::vector<double>> og_arT;
    std::vector<std::set<int>> og_tsk;
    //
    CGLC(Problem *prob, TimeTracker *tt,
        unsigned long time_limit_sec, int numThreads,
        std::string logPath, std::string lpPath, std::string lp_algo,
        std::string _router,
        double dual_gap_limit, unsigned int num_iter_limit, unsigned int no_improvement_limit) : LRH(prob, tt, time_limit_sec, numThreads, logPath, lpPath, lp_algo, _router, dual_gap_limit, num_iter_limit, no_improvement_limit) {
        
        rm_y_ary = new_ak_inv(prob, env, 'I');
        rm_th_w = new IloNumVarArray(env);
        rm_ONE_cnsts = new IloRangeArray(env);
        rm_ONE_cnsts_index = new long*[prob->A.size()];
        for (int a: prob->A) {
            rm_ONE_cnsts_index[a] = new long[prob->E_a[a].size()];
            std::vector<int> a_og;
            std::vector<std::vector<int>> ae_og;
            for (int e = 0; e < prob->E_a[a].size(); e++) {
                std::vector<int> e_og;
                ae_og.push_back(e_og);
            }
            og_a.push_back(a_og);
            og_ae.push_back(ae_og);
        }
        rm_build();
    }
    void rm_build();
    Solution* solve();
    ~CGLC() {
        del_ak_inv(prob, rm_y_ary);
        delete rm_th_w;
        delete rm_ONE_cnsts;
        for (int a: prob->A) {
            delete [] rm_ONE_cnsts_index[a];
        }
        delete [] rm_ONE_cnsts_index;
        //
        delete rmCplex;
        delete rmModel;
    }
    
private:
    void init_cols();
    void solve_rmModel();
    void rm_update();
    void rm_getSol(Solution *sol);
};


#endif /* CGLC_h */
