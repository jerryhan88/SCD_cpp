//
//  Problem.hpp
//  SCD_cpp
//
//  Created by Chung-kyun HAN on 16/2/20.
//  Copyright © 2020 Chung-kyun HAN. All rights reserved.
//

#ifndef Problem_hpp
#define Problem_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <map>

#include "Etc.hpp"
#include "nlohmann/json.hpp"

class Problem {
public:
    std::string problemName;
    //
    std::vector<int> K;
    int *h_k, *n_k;
    double *v_k, *w_k, *r_k;
    //
    std::vector<int> A;
    double *v_a, *w_a;
    std::vector<std::vector<int>> E_a;
    //
    double **p_ae, **u_ae, **o_ae, **d_ae;
    std::vector<std::vector<std::set<int>>> R_ae, K_ae, P_ae, D_ae, N_ae;
    //
    std::vector<int> cN;
    double **t_ij;
    double *al_i, *be_i;
    int ****c_aeij;
    double M;
    
    static Problem* read_json(std::string ifpath) {
        std::ifstream is(ifpath);
        nlohmann::json prob_json;
        is >> prob_json;
        is.close();
        //
        Problem *prob = new Problem();
        prob->problemName = prob_json["problemName"];
        //
        std::string delimiter("&");
        int indexCounter = 0;
        for (auto it: prob_json["t_ij"].items()) {
            std::string loc0loc1(it.key());
            std::vector<std::string> tokens = parseWithDelimiter(loc0loc1, delimiter);
            if (tokens.size() != 2) {
                throw "check t_ij key!";
            }
            for (std::string loc: tokens) {
                if (prob->locID.find(loc) == prob->locID.end()) {
                    prob->locID[loc] = indexCounter;
                    prob->idLoc[indexCounter] = loc;
                    indexCounter++;
                }
            }
        }
        //
        size_t numNodes = prob->locID.size();
        prob->al_i = new double[numNodes];
        prob->be_i = new double[numNodes];
        prob->ga_i = new double[numNodes];
        for (auto it: prob_json["al_i"].items()) {
            std::string loc(it.key());
            prob->al_i[prob->locID[loc]] = it.value();
        }
        for (auto it: prob_json["be_i"].items()) {
            std::string loc(it.key());
            prob->be_i[prob->locID[loc]] = it.value();
        }
        for (auto it: prob_json["ga_i"].items()) {
            std::string loc(it.key());
            prob->ga_i[prob->locID[loc]] = it.value();
        }
        prob->t_ij = new double*[numNodes];
        for (int i = 0; i < numNodes; i++) {
            prob->t_ij[i] = new double[numNodes];
        }
        for (auto it: prob_json["t_ij"].items()) {
            std::string loc0loc1(it.key());
            std::vector<std::string> tokens = parseWithDelimiter(loc0loc1, delimiter);
            int n0 = prob->locID[tokens[0]];
            int n1 = prob->locID[tokens[1]];
            prob->t_ij[n0][n1] = it.value();
        }
        //
        for (int h: prob_json["H"]) {
            prob->H.push_back(h);
        }
        for (int k: prob_json["K"]) {
            prob->K.push_back(k);
        }
        for (int a: prob_json["A"]) {
            prob->A.push_back(a);
        }
        for (std::string n: prob_json["N"]) {
            prob->N.push_back(prob->locID[n]);
        }
        //
        size_t numTasks = prob->K.size();
        prob->r_k = new double[numTasks];
        prob->v_k = new double[numTasks];
        prob->w_k = new double[numTasks];
        prob->h_k = new int[numTasks];
        prob->n_k = new int[numTasks];
        for (int k: prob->K) {
            prob->r_k[k] = prob_json["r_k"][k];
            prob->v_k[k] = prob_json["v_k"][k];
            prob->w_k[k] = prob_json["w_k"][k];
            prob->h_k[k] = prob->locID[prob_json["h_k"][k]];
            prob->n_k[k] = prob->locID[prob_json["n_k"][k]];
        }
        //
        size_t numAgents = prob->A.size();
        prob->v_a = new double[numAgents];
        prob->w_a = new double[numAgents];
        prob->p_ae = new double*[numAgents];
        prob->u_ae = new double*[numAgents];
        prob->l_ae = new double*[numAgents];
        prob->c_aeij = new int***[numAgents];
        for (int a: prob->A) {
            prob->v_a[a] = prob_json["v_a"][a];
            prob->w_a[a] = prob_json["w_a"][a];
            std::vector<int> aE;
            std::vector<std::set<int>> a_S, a_N, a_F, a_iF;
            for (int e: prob_json["E_a"][a]) {
                aE.push_back(e);
            }
            prob->E_a.push_back(aE);
            prob->S_ae.push_back(a_S);
            prob->N_ae.push_back(a_N);
            prob->F_ae.push_back(a_F);
            prob->iF_ae.push_back(a_iF);
            //
            size_t numRR = aE.size();
            prob->p_ae[a] = new double[numRR];
            prob->u_ae[a] = new double[numRR];
            prob->l_ae[a] = new double[numRR];
            prob->c_aeij[a] = new int**[numRR];
            for (int e: prob->E_a[a]) {
                prob->c_aeij[a][e] = new int*[numNodes];
                for (int i = 0; i < numNodes; i++) {
                    prob->c_aeij[a][e][i] = new int[numNodes];
                }
                std::set<int> ae_S, ae_N, ae_F, ae_iF;
                prob->S_ae[a].push_back(ae_S);
                prob->N_ae[a].push_back(ae_N);
                prob->F_ae[a].push_back(ae_F);
                prob->iF_ae[a].push_back(ae_iF);
            }
        }
        //
        for (auto it: prob_json["S_ae"].items()) {
            std::string _ae(it.key());
            std::vector<std::string> tokens = parseWithDelimiter(_ae, delimiter);
            if (tokens.size() != 2) {
                throw "check t_ij key!";
            }
            int a = std::stoi(tokens[0]);
            int e = std::stoi(tokens[1]);
            for (std::string loc: prob_json["S_ae"][_ae]) {
                prob->S_ae[a][e].insert(prob->locID[loc]);
            }
            for (std::string loc: prob_json["N_ae"][_ae]) {
                prob->N_ae[a][e].insert(prob->locID[loc]);
            }
            for (int k: prob_json["F_ae"][_ae]) {
                prob->F_ae[a][e].insert(k);
            }
            prob->p_ae[a][e] = prob_json["p_ae"][_ae];
            prob->l_ae[a][e] = prob_json["l_ae"][_ae];
            prob->u_ae[a][e] = prob_json["u_ae"][_ae];
        }
        //
        for (auto it: prob_json["c_aeij"].items()) {
            std::string _aeij(it.key());
            std::vector<std::string> tokens = parseWithDelimiter(_aeij, delimiter);
            if (tokens.size() != 4) {
                throw "check t_ij key!";
            }
            int a = std::stoi(tokens[0]);
            int e = std::stoi(tokens[1]);
            int i = prob->locID[tokens[2]];
            int j = prob->locID[tokens[3]];
            prob->c_aeij[a][e][i][j] = prob_json["c_aeij"][_aeij];
        }
        //
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                for (int k: prob->K) {
                    if (prob->F_ae[a][e].find(k) == prob->F_ae[a][e].end()) {
                        prob->iF_ae[a][e].insert(k);
                    }
                }
            }
        }
        //
        double maxTravelTime = 0.0, maxHandlingTime = 0.0;
        for (int i = 0; i < numNodes; i++) {
            for (int j = 0; j < numNodes; j++) {
                if (maxTravelTime < prob->t_ij[i][j]) {
                    maxTravelTime = prob->t_ij[i][j];
                }
            }
            if (maxHandlingTime < prob->ga_i[i]) {
                maxHandlingTime = prob->ga_i[i];
            }
        }
        prob->M = prob->N.size() * (maxTravelTime + maxHandlingTime);
        //
        return prob;
    }
};
#endif /* Problem_hpp */
