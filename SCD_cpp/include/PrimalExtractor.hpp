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
#include "Solution.hpp"
#include "SolApprBase.hpp"

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
    PrimalExtractor(Problem *prob, double ****lrh_x_aeij, double ***lrh_u_aei) {
        this->prob = prob;
        this->lrh_x_aeij = lrh_x_aeij;
        this->lrh_u_aei = lrh_u_aei;
        //
    }
    ~PrimalExtractor() { }
    virtual void build() {
        throw "Should override build()";
    }
    virtual void update() {
        throw "Should override update()";
    }
    virtual void getSol(Solution *sol) {
        throw "Should override getSol()";
    }
};

class RouteFixPE: public PrimalExtractor {
public:
    IloNumVar **pex_y_ak, ***pex_z_aek;
    IloRangeArray *pex_COM_cnsts;
    long ***pex_COM_cnsts_index;
    //
    RouteFixPE(Problem *prob, double ****lrh_x_aeij, double ***lrh_u_aei): PrimalExtractor(prob, lrh_x_aeij, lrh_u_aei) {
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
    ~RouteFixPE() {
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
    //
    void build();
    void update();
    void getSol(Solution *sol);
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
    ColGenPE(Problem *prob, double ****lrh_x_aeij, double ***lrh_u_aei): PrimalExtractor(prob, lrh_x_aeij, lrh_u_aei) {
        pex_y_ak = new_inv_ak(prob, env, 'I');
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
    ~ColGenPE() {
        delete_inv_ak(prob, pex_y_ak);
        delete pex_th_w;
        delete pex_ONE_cnsts;
        for (int a: prob->A) {
            delete [] pex_ONE_cnsts_index[a];
        }
        delete [] pex_ONE_cnsts_index;
        //
        og.clear(); og_a.clear(); og_ae.clear(); p_w.clear(); e_wk.clear();
        og_rut.clear(); og_arT.clear(); og_tsk.clear();
        //
        delete pexCplex; delete pexModel;
        env.end();
    }
    
    void build();
    void update();
    void getSol(Solution *sol);
};
#endif /* PrimalExtractor_hpp */
