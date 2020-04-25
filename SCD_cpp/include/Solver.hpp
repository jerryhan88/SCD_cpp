//
//  Solver.hpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 13/4/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef Solver_h
#define Solver_h


#include <tuple>
#include <set>
#include <vector>
//
#include "SolApprBase.hpp"
#include "Problem.hpp"
#include "Solution.hpp"
#include "Allocator.hpp"
#include "RouterSCD.hpp"
#include "PrimalExtractor.hpp"
//
#include "ck_route/Sequencer.h"     // from BnC_CPLEX
#include "ck_route/CutBase.hpp"     // from BnC_CPLEX
#include "ck_route/Router.hpp"     // from BnC_CPLEX

typedef std::tuple<int, int, int, double> iiidTup;
typedef std::tuple<int, double> idTup;

std::vector <int> getIntersection(std::vector < std::vector <int> > &sets);

class Agent {
public:
    Problem *prob;
    int a;
    double volumeLimit, weightLimit;
    std::vector<double> u_e;
    //
    std::set<int> KnM, KnP;
    std::vector<std::set<int>> KnM_e, KnP_e, Hn_e;
    std::vector<std::vector<int>> Sn_e;
    std::vector<int> index4Search_e;
    double currentReward, currentVolume, currentWeight;
    //
    int best_tid;
    double best_expReward;
    std::vector<std::vector<int>> best_tid_eBestSeq;
    //
    Agent(Problem *prob, int a) {
        this->prob = prob;
        this->a = a;
        this->volumeLimit = prob->v_a[a];
        this->weightLimit = prob->w_a[a];
        //
        this->currentReward = 0.0;
        this->currentVolume = 0.0;
        this->currentWeight = 0.0;
        //
        std::vector < std::vector <int> > sets;
        for (int e : prob->E_a[a]) {
            std::vector <int> set(prob->K_ae[a][e].begin(), prob->K_ae[a][e].end());
            sets.push_back(set);
        }
        std::vector <int> feasibleTasks4allRR = getIntersection(sets);
        for (int e : prob->E_a[a]) {
            std::set<int> e_KnM, e_KnP, e_Hn;
            std::vector<int> e_Sn;
            for (int k: feasibleTasks4allRR) {
                KnM.insert(k);
                e_KnM.insert(k);
            }
            KnM_e.push_back(e_KnM);
            KnP_e.push_back(e_KnP);
            Hn_e.push_back(e_Hn);
            for (int i: prob->R_ae[a][e]) {
                e_Sn.push_back(i);
            }
            Sn_e.push_back(e_Sn);
            index4Search_e.push_back(1);
            u_e.push_back(prob->u_ae[a][e]);
        }
    }
    ~Agent() {}
    //
    void find_bestTask();
    void updateBestTask();
    void updateFeasibleTasks();
};

class PureGH : public SolApprBase {
public:
    std::set<int> unSelectedTasks;
    std::vector<Agent*> agents;
    //
    PureGH(Problem *prob, TimeTracker *tt, unsigned long time_limit_sec, int numThreads, std::string logPath): SolApprBase(prob, tt, time_limit_sec, numThreads, logPath) {
        for (int k: prob->K) {
            unSelectedTasks.insert(k);
        }
        for (int a: prob->A) {
            agents.push_back(new Agent(prob, a));
        }
    }
    ~PureGH() {
    }
    //
    Solution* solve();
    bool confirmTaskSelection(Agent *agt);
};

    
class ILP : public BaseMM {
public:
    ILP(Problem *prob, TimeTracker *tt,
        unsigned long time_limit_sec, int numThreads,
        std::string logPath, std::string lpPath) : BaseMM(prob, tt, time_limit_sec, numThreads, logPath, lpPath, 'I') {
        this->lpPath = lpPath;
    }
    ~ILP() {
    }
    Solution* solve();
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
        std::string logPath, std::string lpPath,
        std::string _router, std::string _extractor,
        double dual_gap_limit, unsigned int num_iter_limit, unsigned int no_improvement_limit) : BaseMM(prob, tt, time_limit_sec, numThreads, logPath, lpPath, 'C') {
        _ROUTER = _router;
        _EXTRACTOR = _extractor;
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
        if (_EXTRACTOR == "RF") {
            pex = new RouteFixPE(prob, lrh_x_aeij, lrh_u_aei);
        } else {
            assert (_EXTRACTOR == "CG");
            pex = new ColGenPE(prob, lrh_x_aeij, lrh_u_aei);
        }
        alr = new Allocator(prob, lrh_l_aek);
        build_rutModels();
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
private:
    void solve_dualProblem();
    bool updateLMs();
    //
    void solve_etaModel();
    //
    void build_rutModels();
    void solve_rutModels();
    //
    void solve_pexModel();
    //
    void logging(std::string indicator, std::string note);
};

#endif /* Solver_h */
