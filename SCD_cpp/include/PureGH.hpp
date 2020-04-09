//
//  PureGH.hpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 29/3/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef PureGH_hpp
#define PureGH_hpp

#include <set>
#include <vector>
//
#include "SolApprBase.hpp"
#include "Problem.hpp"
#include "Solution.hpp"
//
#include "ck_route/Sequencer.h"     // from BnC_CPLEX


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



#endif /* PureGH_hpp */
