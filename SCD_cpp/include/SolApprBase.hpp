//
//  SolBase.hpp
//  SCD_cpp
//
//  Created by Chung-Kyun HAN on 29/3/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef SolApprBase_hpp
#define SolApprBase_hpp

#include <cfloat>
#include <math.h>
#include <set>
#include <map>

#include <ilcplex/ilocplex.h>

#include "Problem.hpp"
#include "Solution.hpp"
#include "ThreadPool.hpp"

#include "ck_util/util.hpp"         // from util

#define DEFAULT_BUFFER_SIZE 2048
//#define LP_TIME_LIMIT 1800

#define LP_TIME_LIMIT 300


class SolApprBase {
public:
    Problem *prob;
    TimeTracker *tt;
    unsigned long time_limit_sec;
    std::string logPath;
    ThreadPool &routerPool = ThreadPool::getInstance(1);
    //
    SolApprBase(Problem *prob, TimeTracker *tt, unsigned long time_limit_sec, int numThreads, std::string logPath) {
        this->prob = prob;
        this->logPath = logPath;
        this->tt = tt;
        this->time_limit_sec = time_limit_sec;
        if (numThreads > 1) {
            routerPool.resize(numThreads);
        }
    }
    ~SolApprBase() {}
    //
    virtual Solution* solve() {
        throw "Should override solve()";
    }
};

class BaseMM : public SolApprBase {
public:
    std::string lpPath;
    std::string pathPrefix;
    std::string lp_algo;
    //
    IloEnv env;
    IloModel *baseModel;
    IloCplex *baseCplex;
    IloNumVar **y_ary, ***z_ary, ****x_ary, ***u_ary;
    std::map<iiTup, IloNumVar> y_map;
    std::map<iiiTup, IloNumVar> z_map;
    std::map<iiiiTup, IloNumVar> x_map;
    std::map<iiiTup, IloNumVar> u_map;
    IloRangeArray *COM_cnsts;
    long ***COM_cnsts_index;
    bool ****bool_x_aeij;
    //
    BaseMM(Problem *prob, TimeTracker *tt,
           unsigned long time_limit_sec, int numThreads,
           std::string logPath, std::string lpPath,
           std::string lp_algo, char vType): SolApprBase(prob, tt, time_limit_sec, numThreads, logPath) {
        if (lp_algo != "DSM" && lp_algo != "IPM") {
            assert(false);
        }
        this->lpPath = lpPath;
        this->pathPrefix = lpPath.substr(0, lpPath.find(".lp"));
        this->lp_algo = lp_algo;
        //
        bool_x_aeij = new bool***[prob->A.size()];
        for (int a: prob->A) {
            bool_x_aeij[a] = new bool**[prob->E_a[a].size()];
            for (int e: prob->E_a[a]) {
                bool_x_aeij[a][e] = new bool*[prob->cN.size()];
                for (int i = 0; i < prob->cN.size(); i++) {
                    bool_x_aeij[a][e][i] = new bool[prob->cN.size()];
                }
                for (int i: prob->N_ae[a][e]) {
                    for (int j: prob->N_ae[a][e]) {
                        bool_x_aeij[a][e][i][j] = true;
                    }
                }
                
            }
        }
        preprocessing();
        //
        generate_y_inv(vType);
        generate_z_inv(vType);
        generate_x_inv(vType);
        generate_u_inv();
        //
        COM_cnsts = new IloRangeArray(env);
        COM_cnsts_index = new long**[prob->A.size()];
        for (int a : prob->A) {
            COM_cnsts_index[a] = new long*[prob->E_a[a].size()];
            for (int e: prob->E_a[a]) {
                COM_cnsts_index[a][e] = new long[prob->K.size()];
            }
        }
        //
        baseModel = new IloModel(env);
        build_baseModel();
        baseCplex = new IloCplex(*baseModel);
        baseCplex->setOut(env.getNullStream());
        if (lp_algo == "IPM") {
            baseCplex->setParam(IloCplex::Param::RootAlgorithm, IloCplex::Barrier);
        }
    }
    ~BaseMM() {
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                for (int i = 0; i < prob->cN.size(); i++) {
                    delete [] bool_x_aeij[a][e][i];
                }
                delete [] bool_x_aeij[a][e];
            }
            delete [] bool_x_aeij[a];
        }
        delete [] bool_x_aeij;
        //
        delete_y_inv();
        delete_z_inv();
        delete_x_inv();
        delete_u_inv();
        //
        delete baseCplex;
        delete baseModel;
    }
    //
    static void def_ETA_cnsts(Problem *prob, IloEnv &env, IloNumVar **y_ary, IloNumVar ***z_ary, IloModel *model) {
        //
        // Evaluation of the Task Assignment
        //
        char buf[2048];
        IloRangeArray cnsts(env);
        IloExpr linExpr(env);
        //
        for (int k : prob->K) {
            linExpr.clear();
            sprintf(buf, "TAS(%d)", k);  // Task Assignment
            for (int a : prob->A) {
                linExpr += y_ary[a][k];
            }
            cnsts.add(linExpr <= 1);
            cnsts[cnsts.getSize() - 1].setName(buf);
        }
        //
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                for (int k: prob->K) {
                    if (prob->K_ae[a][e].find(k) == prob->K_ae[a][e].end()) {
                        linExpr.clear();
                        sprintf(buf, "ITA(%d)(%d)(%d)", a, e, k);  // Infeasible Tasks Assignment
                        linExpr += y_ary[a][k];
                        cnsts.add(linExpr == 0);
                        cnsts[cnsts.getSize() - 1].setName(buf);
                    }
                }
            }
        }
        //
        for (int a: prob->A) {
            linExpr.clear();
            sprintf(buf, "VL(%d)", a);  // Volume Limit
            for (int k: prob->K) {
                linExpr += prob->v_k[k] * y_ary[a][k];
            }
            cnsts.add(linExpr <= prob->v_a[a]);
            cnsts[cnsts.getSize() - 1].setName(buf);
        }
        //
        for (int a: prob->A) {
            linExpr.clear();
            sprintf(buf, "WL(%d)", a);  // Weight Limit
            for (int k: prob->K) {
                linExpr += prob->w_k[k] * y_ary[a][k];
            }
            cnsts.add(linExpr <= prob->w_a[a]);
            cnsts[cnsts.getSize() - 1].setName(buf);
        }
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                for (int k: prob->K) {
                    linExpr.clear();
                    sprintf(buf, "TAC(%d)(%d)(%d)", a, e, k);  // Task Accomplishment
                    linExpr += z_ary[a][e][k];
                    linExpr -= y_ary[a][k];
                    cnsts.add(linExpr <= 0);
                    cnsts[cnsts.getSize() - 1].setName(buf);
                }
            }
        }
        model->add(cnsts);
    }
protected:
    void generate_y_inv(char vType);
    void generate_z_inv(char vType);
    void generate_x_inv(char vType);
    void generate_u_inv();
    //
    void delete_y_inv();
    void delete_z_inv();
    void delete_x_inv();
    void delete_u_inv();
    //
    IloNumVar get_y_inv(int a, int k);
    IloNumVar get_z_inv(int a, int e, int k);
    IloNumVar get_x_inv(int a, int e, int i, int j);
    IloNumVar get_u_inv(int a, int e, int i);
    //
    void build_baseModel();
    void preprocessing();
    void def_objF();
    void def_ETA_cnsts();
    void def_RUT_cnsts();
    void def_FC_cnsts_aeGiven(int a, int e);
    void def_AT_cnsts_aeGiven(int a, int e);
    void def_COM_cnsts();
};


#endif /* SolApprBase_h */
