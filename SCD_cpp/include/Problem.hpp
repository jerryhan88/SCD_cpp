//
//  Problem.hpp
//  SCD_cpp
//
//  Created by Chung-kyun HAN on 16/2/20.
//  Copyright © 2020 Chung-kyun HAN. All rights reserved.
//

#ifndef Problem_hpp
#define Problem_hpp

#include <vector>
#include <set>
#include <map>

#include "ck_util/util.hpp"
#include "nlohmann/json.hpp"
#include "ck_route/Other.hpp"     // from BnC_CPLEX

class Problem;

double** new_dbl_ak(Problem *prob);
double*** new_dbl_aek(Problem *prob);
double**** new_dbl_aeij(Problem *prob);
double*** new_dbl_aei(Problem *prob);

void delete_dbl_ak(Problem *prob, double **dbl_ak);
void delete_dbl_aek(Problem *prob, double ***dbl_aek);
void delete_dbl_aeij(Problem *prob, double ****dbl_aeij);
void delete_dbl_aei(Problem *prob, double ***dbl_aei);

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
    double **p_ae, **u_ae;
    int **o_ae, **d_ae;
    std::vector<std::vector<std::vector<int>>> R_ae;
    std::vector<std::vector<std::set<int>>> K_ae, P_ae, D_ae, PD_ae, N_ae;
    std::vector<std::vector<rut::Problem*>> RP_ae;
    std::vector<std::vector<std::map<int, int>>> scdID_rutID_ae, rutID_scdID_ae;
    std::map<int, int> dp_k;
    //
    std::vector<int> cN;
    double **t_ij;
    double *al_i, *be_i;
    int ****c_aeij;
    double M;
    //
    ~Problem() {
        std::vector<int*> ipV = {h_k, n_k};
        for(auto p: ipV) {
            delete [] p;
        }
        //
        std::vector<double*> dpV = {r_k, v_k, w_k, v_a, w_a, al_i, be_i};
        for(auto p: dpV) {
            delete [] p;
        }
        //
        for (int a: A) {
            delete [] p_ae[a];
            delete [] u_ae[a];
            delete [] o_ae[a];
            delete [] d_ae[a];
        }
        delete [] p_ae;
        delete [] u_ae;
        delete [] o_ae;
        delete [] d_ae;
        //
        for (int i = 0; i < cN.size(); i++) {
            delete [] t_ij[i];
        }
        delete [] t_ij;
        //
        for (int a: A) {
            for (int e: E_a[a]) {
                for (int i = 0; i < cN.size(); i++) {
                    delete [] c_aeij[a][e][i];
                }
                delete [] c_aeij[a][e];
                delete RP_ae[a][e];
            }
            delete [] c_aeij[a];
        }
        delete [] c_aeij;
        //
        R_ae.clear(); K_ae.clear(); P_ae.clear(); D_ae.clear(); PD_ae.clear(); N_ae.clear();
        K.clear(); A.clear(); E_a.clear(); cN.clear();
        RP_ae.clear();
        dp_k.clear();
    }
    
    void gen_aeProbs();
    
    static Problem* read_json(std::string ifpath) {
        std::ifstream is(ifpath);
        nlohmann::json prob_json;
        is >> prob_json;
        is.close();
        //
        Problem *prob = new Problem();
        prob->problemName = prob_json["problemName"];
        //
        for (int k: prob_json["K"]) {
            prob->K.push_back(k);
        }
        size_t numTasks = prob->K.size();
        prob->r_k = new double[numTasks];
        prob->v_k = new double[numTasks];
        prob->w_k = new double[numTasks];
        prob->h_k = new int[numTasks];
        prob->n_k = new int[numTasks];
        for (int k: prob->K) {
            int pp = prob_json["h_k"][k];
            int dp = prob_json["n_k"][k];
            prob->r_k[k] = prob_json["r_k"][k];
            prob->v_k[k] = prob_json["v_k"][k];
            prob->w_k[k] = prob_json["w_k"][k];
            prob->h_k[k] = pp;
            prob->n_k[k] = dp;
            //
            prob->dp_k[pp] = -1;
            prob->dp_k[dp] = k;
        }
        //
        for (int a: prob_json["A"]) {
            prob->A.push_back(a);
        }
        size_t numAgents = prob->A.size();
        prob->v_a = new double[numAgents];
        prob->w_a = new double[numAgents];
        for (int a: prob->A) {
            prob->v_a[a] = prob_json["v_a"][a];
            prob->w_a[a] = prob_json["w_a"][a];
            std::vector<int> aE;
            for (int e: prob_json["E_a"][a]) {
                aE.push_back(e);
            }
            prob->E_a.push_back(aE);
        }
        prob->p_ae = new double*[numAgents];
        prob->u_ae = new double*[numAgents];
        prob->o_ae = new int*[numAgents];
        prob->d_ae = new int*[numAgents];
        for (int a: prob->A) {
            size_t numRR = prob->E_a[a].size();
            prob->p_ae[a] = new double[numRR];
            prob->u_ae[a] = new double[numRR];
            prob->o_ae[a] = new int[numRR];
            prob->d_ae[a] = new int[numRR];
            for (int e: prob->E_a[a]) {
                prob->p_ae[a][e] = prob_json["p_ae"][a][e];
                prob->u_ae[a][e] = prob_json["u_ae"][a][e];
                prob->o_ae[a][e] = prob_json["o_ae"][a][e];
                prob->d_ae[a][e] = prob_json["d_ae"][a][e];
            }
        }
        for (int a: prob->A) {
            std::vector<std::vector<int>> aR;
            std::vector<std::set<int>> aK, aP, aD, aPD, aN;
            for (int e: prob->E_a[a]) {
                std::vector<int> aeR;
                std::set<int> aeK, aeP, aeD, aePD, aeN;
                for (int i: prob_json["R_ae"][a][e]) {
                    aeR.push_back(i);
                    prob->dp_k[i] = -1;
                }
                for (int i: prob_json["K_ae"][a][e]) {
                    aeK.insert(i);
                }
                for (int i: prob_json["P_ae"][a][e]) {
                    aeP.insert(i);
                }
                for (int i: prob_json["D_ae"][a][e]) {
                    aeD.insert(i);
                }
                for (int i: prob_json["PD_ae"][a][e]) {
                    aePD.insert(i);
                }
                for (int i: prob_json["N_ae"][a][e]) {
                    aeN.insert(i);
                }
                aR.push_back(aeR);
                aK.push_back(aeK);
                aP.push_back(aeP);
                aD.push_back(aeD);
                aPD.push_back(aePD);
                aN.push_back(aeN);
            }
            prob->R_ae.push_back(aR);
            prob->K_ae.push_back(aK);
            prob->P_ae.push_back(aP);
            prob->D_ae.push_back(aD);
            prob->PD_ae.push_back(aPD);
            prob->N_ae.push_back(aN);
        }
        //
        for (int i: prob_json["cN"]) {
            prob->cN.push_back(i);
        }
        size_t numNodes = prob->cN.size();
        prob->al_i = new double[numNodes];
        prob->be_i = new double[numNodes];
        prob->t_ij = new double*[numNodes];
        for (int i = 0; i < numNodes; i++) {
            prob->t_ij[i] = new double[numNodes];
        }
        prob->c_aeij = new int***[numAgents];
        for (int a: prob->A) {
            size_t numRR = prob->E_a[a].size();
            prob->c_aeij[a] = new int**[numRR];
            for (int e: prob->E_a[a]) {
                prob->c_aeij[a][e] = new int*[numNodes];
                for (int i = 0; i < numNodes; i++) {
                    prob->c_aeij[a][e][i] = new int[numNodes];
                    for (int j = 0; j < numNodes; j++) {
                        prob->c_aeij[a][e][i][j] = 0;
                    }
                }
            }
        }
        for (int i: prob->cN) {
            prob->al_i[i] = prob_json["al_i"][i];
            prob->be_i[i] = prob_json["be_i"][i];
            for (int j: prob->cN) {
                prob->t_ij[i][j] = prob_json["t_ij"][i][j];
            }
        }
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                for (int i0 = 0; i0 < prob->R_ae[a][e].size(); i0++) {
                    int locID0 = prob->R_ae[a][e][i0];
                    for (int i1 = i0; i1 < prob->R_ae[a][e].size(); i1++) {
                        int locID1 = prob->R_ae[a][e][i1];
                        prob->c_aeij[a][e][locID0][locID1] = 1;
                        prob->c_aeij[a][e][locID1][locID0] = 0;
                    }
                }
            }
        }
        //
        prob->M = prob_json["M"];
        //
        return prob;
    }
};
#endif /* Problem_hpp */
