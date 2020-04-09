//
//  PureGH.cpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 29/3/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/PureGH.hpp"

std::vector <int> getIntersection(std::vector < std::vector <int> > &sets) {
    std::vector <int> result;  // To store the reaultant set
    int smallSetInd = 0;  // Initialize index of smallest set
    unsigned long minSize = sets[0].size(); // Initialize size of smallest set

    // sort all the sets, and also find the smallest set
    for (int i = 1 ; i < sets.size() ; i++)
    {
        // sort this set
        std::sort(sets[i].begin(), sets[i].end());

        // update minSize, if needed
        if (minSize > sets[i].size())
        {
            minSize = sets[i].size();
            smallSetInd = i;
        }
    }

    std::map<int,int> elementsMap;

    // Add all the elements of smallest set to a map, if already present,
    // update the frequency
    for (int i = 0; i < sets[smallSetInd].size(); i++)
    {
        if (elementsMap.find( sets[smallSetInd][i] ) == elementsMap.end())
            elementsMap[ sets[smallSetInd][i] ] = 1;
        else
            elementsMap[ sets[smallSetInd][i] ]++;
    }

    // iterate through the map elements to see if they are present in
    // remaining sets
    std::map<int,int>::iterator it;
    for (it = elementsMap.begin(); it != elementsMap.end(); ++it)
    {
        int elem = it->first;
        long freq = it->second;

        bool bFound = true;

        // Iterate through all sets
        for (int j = 0 ; j < sets.size() ; j++)
        {
            // If this set is not the smallest set, then do binary search in it
            if (j != smallSetInd)
            {
                // If the element is found in this set, then find its frequency
                if (binary_search( sets[j].begin(), sets[j].end(), elem ))
                {
                   long lInd = lower_bound(sets[j].begin(), sets[j].end(), elem) - sets[j].begin();
                   long rInd = upper_bound(sets[j].begin(), sets[j].end(), elem) - sets[j].begin();

                   // Update the minimum frequency, if needed
                   if ((rInd - lInd) < freq)
                       freq = rInd - lInd;
                }
                // If the element is not present in any set, then no need
                // to proceed for this element.
                else
                {
                    bFound = false;
                    break;
                }
            }
        }

        // If element was found in all sets, then add it to result 'freq' times
        if (bFound)
        {
            for (int k = 0; k < freq; k++)
                result.push_back(elem);
        }
    }
    return result;
}


Solution* PureGH::solve() {
    auto taskSelFunc = [](Agent* agt){ agt->find_bestTask(); };
    //
    while (!unSelectedTasks.empty()) {
        size_t sizeBefore = unSelectedTasks.size();
        //
        std::vector<std::shared_future<void>> threadJobs;
        for (Agent* agt: agents) {
            std::shared_future<void> returnValue = (std::shared_future<void>) pool.push(taskSelFunc, agt);
            threadJobs.push_back(returnValue);
        }
        for (std::vector<std::shared_future<void>>::iterator it = threadJobs.begin(); it != threadJobs.end(); it++) {
            (*it).get();
        }
        //
        for (Agent* agt: agents) {
            bool successTaskSelection = confirmTaskSelection(agt);
            if (successTaskSelection) {
                agt->updateBestTask();
            } else {
                agt->updateFeasibleTasks();
            }
        }
        if (sizeBefore == unSelectedTasks.size()) {
            break;
        }
    }
    //
    Solution *sol = new Solution(prob);
    sol->cpuT = tt->get_elapsedTimeCPU();
    sol->wallT = tt->get_elapsedTimeWall();
    for (int a : prob->A) {
        for (int k : prob->K) {
            sol->y_ak[a][k] = 0.0;
        }
        for (int e: prob->E_a[a]) {
            for (int k : prob->K) {
                sol->z_aek[a][e][k] = 0.0;
            }
            for (int i: prob->N_ae[a][e]) {
                for (int j: prob->N_ae[a][e]) {
                    sol->x_aeij[a][e][i][j] = 0.0;
                }
                sol->u_aei[a][e][i] = 0.0;
            }
        }
    }
    //
    int n0, n1;
    sol->objV = 0.0;
    Agent *agt;
    for (int a : prob->A) {
        agt = agents[a];
        sol->objV += agt->currentReward;
        for (int k: prob->K) {
            if (agt->KnP.find(k) != agt->KnP.end()) {
                sol->y_ak[a][k] = 1.0;
                for (int e: prob->E_a[a]) {
                    if (agt->KnP_e[e].find(k) != agt->KnP_e[e].end()) {
                        sol->z_aek[a][e][k] = 1.0;
                    }
                }
            }
        }
        for (int e: prob->E_a[a]) {
            n0 = agt->Sn_e[e][0];
            sol->u_aei[a][e][n0] = prob->al_i[n0];
            for (int s = 1; s < agt->Sn_e[e].size(); s++) {
                n1 = agt->Sn_e[e][s];
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
    sol->gap = -1.0;
    return sol;
}

bool PureGH::confirmTaskSelection(Agent *agt) {
    bool successTaskSelection = false;
    if (unSelectedTasks.find(agt->best_tid) != unSelectedTasks.end()) {
        unSelectedTasks.erase(agt->best_tid);
        successTaskSelection = true;
    }
    return successTaskSelection;
}

void Agent::updateBestTask() {
    KnM.erase(best_tid);
    KnP.insert(best_tid);
    //
    currentReward += best_expReward;
    currentVolume += prob->v_k[best_tid];
    currentWeight += prob->w_k[best_tid];
    for (int e: prob->E_a[a]) {
        if (KnM_e[e].find(best_tid) != KnM_e[e].end()) {
            if (best_tid_eBestSeq[e].size() != 0) {
                KnM_e[e].erase(best_tid);
                KnP_e[e].insert(best_tid);
                Sn_e[e].clear();
                Sn_e[e].insert(Sn_e[e].end(), best_tid_eBestSeq[e].begin(), best_tid_eBestSeq[e].end());
                Hn_e[e].insert(prob->h_k[best_tid]);
                std::vector<int>::iterator it = std::find(Sn_e[e].begin(), Sn_e[e].end(), prob->n_k[best_tid]);
                index4Search_e[e] = (int) std::distance(Sn_e[e].begin(), it) + 1;
            }
        }
    }
}

void Agent::updateFeasibleTasks() {
    KnM.erase(best_tid);
}

void Agent::find_bestTask() {
    double kr, kv, kw, p, expReward;
    bool isFeasible;
    //
    best_tid = -1;
    best_expReward = -DBL_MAX;
    best_tid_eBestSeq.clear();
    for (int k: KnM) {
        kr = prob->r_k[k];
        kv = prob->v_k[k];
        kw = prob->w_k[k];
        //
        std::vector<std::vector<int>> tid_eBestSeq;
        expReward = 0.0;
        for (int e: prob->E_a[a]) {
            p = prob->p_ae[a][e];
            tid_eBestSeq.push_back(std::vector<int>());
            //
            isFeasible = true;
            if (KnM_e[e].find(k) == KnM_e[e].end()) {
                isFeasible = false;
            } else if (volumeLimit < currentVolume + kv) {
                isFeasible = false;
            } else if (weightLimit < currentWeight + kw) {
                isFeasible = false;
            } else {
                int partialSeqSize, minSeqSize;
                partialSeqSize = (int) Sn_e[e].size();
                int *minSeq = new int[partialSeqSize + 2];
                int *seq = new int[partialSeqSize + 2];
                if ( Hn_e[e].find( prob->h_k[k] ) != Hn_e[e].end() ) {
                    minSeqSize = partialSeqSize + 1;
                } else {
                    minSeqSize = partialSeqSize + 2;
                }
                int n0 = prob->h_k[k];
                int n1 = prob->n_k[k];
                set_min_tt_seq(&Sn_e[e][0], partialSeqSize,
                                minSeq, seq,
                                minSeqSize,
                                index4Search_e[e],
                                n0, n1,
                                prob->al_i, prob->be_i,
                                prob->t_ij);
                
                if (minSeq[0] == -1) {
                    isFeasible = false;
                } else {
                    double tt = get_travelTime(minSeq, minSeqSize, prob->al_i, prob->be_i, prob->t_ij);
                    if (tt <= u_e[e]) {
                        if (valid_TWs(minSeq, minSeqSize,
                                     prob->al_i, prob->be_i,
                                      prob->t_ij)) {
                            for (int i = 0; i < minSeqSize; i++) {
                                tid_eBestSeq[e].push_back(minSeq[i]);
                            }
                        } else {
                            isFeasible = false;
                        }
                    } else {
                        isFeasible = false;
                    }
                }
                delete [] minSeq;
                delete [] seq;
            }
            if (isFeasible) {
                expReward += kr * p;
            }
        }
        if (best_expReward < expReward) {
            best_expReward = expReward;
            best_tid = k;
            best_tid_eBestSeq.clear();
            for (std::vector<int> v: tid_eBestSeq) {
                best_tid_eBestSeq.push_back(std::vector<int> (v.begin(), v.end()));
            }
        }
    }
}
