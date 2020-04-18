//
//  PrimalExtractor.hpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 13/4/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef PrimalExtractor_hpp
#define PrimalExtractor_hpp

#include <ilcplex/ilocplex.h>

#include "Problem.hpp"
#include "SolApprBase.hpp"

class PrimalExtractor {
public:
    Problem *prob;
    //
    IloEnv env;
    IloModel *pexModel;
    IloCplex *pexCplex;
    //
    double **y_ak, ***z_aek, ****x_aeij, ***u_aei;
    double F_V;
    //
    PrimalExtractor(Problem *prob) {
        this->prob = prob;
        //
        size_t numAgents, numTasks, numNodes;
        numAgents = prob->A.size();
        numTasks = prob->K.size();
        numNodes = prob->cN.size();
        //
        y_ak = new_dbl_ak(prob);
        z_aek = new_dbl_aek(prob);
        x_aeij = new_dbl_aeij(prob);
        u_aei = new_dbl_aei(prob);
    }
    ~PrimalExtractor() {
        delete_dbl_ak(prob, y_ak);
        delete_dbl_aek(prob, z_aek);
        delete_dbl_aeij(prob, x_aeij);
        delete_dbl_aei(prob, u_aei);
    }
    virtual void build() {
        throw "Should override build()";
    }
    virtual void update(double ****lrh_x_aeij, double ***lrh_u_aei) {
        throw "Should override update()";
    }
    virtual void solve(double ****lrh_x_aeij, double ***lrh_u_aei) {
        throw "Should override solve()";
    }
};

class RouteFixPE: public PrimalExtractor {
public:
    IloNumVar **pex_y_ak, ***pex_z_aek;
    IloRangeArray *pex_COM_cnsts;
    long ***pex_COM_cnsts_index;
    //
    RouteFixPE(Problem *prob): PrimalExtractor(prob) {
        pex_y_ak = gen_y_ak(prob, env, 'I');
        pex_z_aek = gen_z_aek(prob, env, 'I');
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
    ~RouteFixPE() {
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                delete [] pex_z_aek[a][e];
                delete [] pex_COM_cnsts_index[a][e];
            }
            delete [] pex_y_ak[a]; delete [] pex_z_aek[a];
            delete [] pex_COM_cnsts_index[a];
        }
        delete [] pex_y_ak; delete [] pex_z_aek;
        delete [] pex_COM_cnsts_index;
        //
        delete pexCplex; delete pexModel;
        delete pex_COM_cnsts;
        env.end();
    }
    //
    void build();
    void update(double ****lrh_x_aeij, double ***lrh_u_aei);
    void solve(double ****lrh_x_aeij, double ***lrh_u_aei);
};

class ColGenPE: public PrimalExtractor {
public:
    IloNumVar **pex_y_ak;
    IloNumVarArray *pex_th_w;
    IloRangeArray *pex_ONE_cnsts;
    long **pex_ONE_cnsts_index;
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
    ColGenPE(Problem *prob): PrimalExtractor(prob) {
        pex_y_ak = gen_y_ak(prob, env, 'I');
        pex_th_w = new IloNumVarArray(env);
        pex_ONE_cnsts = new IloRangeArray(env);
        pex_ONE_cnsts_index = new long*[prob->A.size()];
        for (int a: prob->A) {
            pex_ONE_cnsts_index[a] = new long[prob->E_a[a].size()];
            std::vector<int> a_og;
            std::vector<std::vector<int>> ae_og;
            for (int e = 0; e < prob->E_a[a].size(); e++) {
                std::vector<int> e_og;
                ae_og.push_back(e_og);
            }
            og_a.push_back(a_og);
            og_ae.push_back(ae_og);
        }
        build();
    }
    void build();
    void update(double ****lrh_x_aeij, double ***lrh_u_aei);
    void solve(double ****lrh_x_aeij, double ***lrh_u_aei);
};
#endif /* PrimalExtractor_hpp */
