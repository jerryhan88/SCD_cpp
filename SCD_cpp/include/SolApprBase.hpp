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

#include <ilcplex/ilocplex.h>

#include "Problem.hpp"
#include "Solution.hpp"
#include "ThreadPool.hpp"

#include "ck_util/util.hpp"         // from util

#define DEFAULT_BUFFER_SIZE 2048

class SolApprBase {
public:
    Problem *prob;
    TimeTracker *tt;
    unsigned long time_limit_sec;
    std::string logPath;
    ThreadPool& pool = ThreadPool::getInstance(1);
    //
    SolApprBase(Problem *prob, TimeTracker *tt, unsigned long time_limit_sec, int numThreads, std::string logPath) {
        this->prob = prob;
        this->logPath = logPath;
        this->tt = tt;
        this->time_limit_sec = time_limit_sec;
        if (numThreads > 1) {
            pool.resize(numThreads);
        }
    }
    ~SolApprBase() {}
    //
    virtual Solution* solve() {
        throw "Should override solve()";
    }
};

IloNumVar** gen_y_ak(Problem *prob, IloEnv &env, char vType);
IloNumVar*** gen_z_aek(Problem *prob, IloEnv &env, char vType);
IloNumVar**** gen_x_aeij(Problem *prob, IloEnv &env, char vType);
IloNumVar*** gen_u_aei(Problem *prob, IloEnv &env);

class BaseMM : public SolApprBase {
public:
    std::string lpPath;
    std::string pathPrefix;
    //
    IloEnv env;
    IloModel *baseModel;
    IloCplex *baseCplex;
    IloNumVar **y_ak, ***z_aek, ****x_aeij, ***u_aei;
    IloRangeArray *COM_cnsts;
    long ***COM_cnsts_index;
    //
    BaseMM(Problem *prob, TimeTracker *tt, unsigned long time_limit_sec, int numThreads, std::string logPath, std::string lpPath, char vType): SolApprBase(prob, tt, time_limit_sec, numThreads, logPath) {
        this->lpPath = lpPath;
        this->pathPrefix = lpPath.substr(0, lpPath.find(".lp"));
        //
        y_ak = gen_y_ak(prob, env, vType);
        z_aek = gen_z_aek(prob, env, vType);
        x_aeij = gen_x_aeij(prob, env, vType);
        u_aei = gen_u_aei(prob, env);
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
    }
    ~BaseMM() {
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                for (int i = 0; i < prob->cN.size(); i++) {
                    delete [] x_aeij[a][e][i];
                }
                delete [] z_aek[a][e]; delete [] x_aeij[a][e]; delete [] u_aei[a][e];
            }
            delete [] y_ak[a]; delete [] z_aek[a]; delete [] x_aeij[a]; delete [] u_aei[a];
        }
        delete [] y_ak; delete [] z_aek; delete [] x_aeij; delete [] u_aei;
        //
        delete baseCplex;
        delete baseModel;
    }
    //
    static void def_ETA_cnsts(Problem *prob, IloEnv &env, IloNumVar **y_ak, IloNumVar ***z_aek, IloModel *model) {
        //
        // Evaluation of the Task Assignment
        //
        char buf[DEFAULT_BUFFER_SIZE];
        IloRangeArray cnsts(env);
        IloExpr linExpr(env);
        //
        for (int k : prob->K) {
            linExpr.clear();
            sprintf(buf, "TAS(%d)", k);  // Task Assignment
            for (int a : prob->A) {
                linExpr += y_ak[a][k];
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
                        linExpr += y_ak[a][k];
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
                linExpr += prob->v_k[k] * y_ak[a][k];
            }
            cnsts.add(linExpr <= prob->v_a[a]);
            cnsts[cnsts.getSize() - 1].setName(buf);
        }
        //
        for (int a: prob->A) {
            linExpr.clear();
            sprintf(buf, "WL(%d)", a);  // Weight Limit
            for (int k: prob->K) {
                linExpr += prob->w_k[k] * y_ak[a][k];
            }
            cnsts.add(linExpr <= prob->w_a[a]);
            cnsts[cnsts.getSize() - 1].setName(buf);
        }
        for (int a: prob->A) {
            for (int e: prob->E_a[a]) {
                for (int k: prob->K) {
                    linExpr.clear();
                    sprintf(buf, "TAC(%d)(%d)(%d)", a, e, k);  // Task Accomplishment
                    linExpr += z_aek[a][e][k];
                    linExpr -= y_ak[a][k];
                    cnsts.add(linExpr <= 0);
                    cnsts[cnsts.getSize() - 1].setName(buf);
                }
            }
        }
        model->add(cnsts);
    }
    
protected:
    void build_baseModel();
    void def_objF();
    void def_RUT_cnsts();
    void def_FC_cnsts_aeGiven(int a, int e);
    void def_AT_cnsts_aeGiven(int a, int e);
    void def_COM_cnsts();
};


#endif /* SolApprBase_h */
