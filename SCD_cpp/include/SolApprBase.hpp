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
#include <tuple>

#include <ilcplex/ilocplex.h>

#include "Problem.hpp"
#include "Solution.hpp"
#include "ThreadPool.hpp"

#include "ck_util/util.hpp"         // from util

#define DEFAULT_BUFFER_SIZE 2048
#define LP_TIME_LIMIT 1800

typedef std::tuple<int, int, int, int> iiiiTup;
typedef std::tuple<int, int, int, double> iiidTup;
typedef std::tuple<int, double> idTup;



class SolApprBase {
public:
    Problem *prob;
    TimeTracker *tt;
    unsigned long time_limit_sec;
    std::string logPath;
//    ThreadPool &builderPool = ThreadPool::getInstance(8);
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

IloNumVar** new_inv_ak(Problem *prob, IloEnv &env, char vType);
IloNumVar*** new_inv_aek(Problem *prob, IloEnv &env, char vType);
IloNumVar**** new_inv_aeij(Problem *prob, IloEnv &env, char vType);
IloNumVar*** new_inv_aei(Problem *prob, IloEnv &env);


void delete_inv_ak(Problem *prob, IloNumVar **inv_ak);
void delete_inv_aek(Problem *prob, IloNumVar ***inv_aek);
void delete_inv_aeij(Problem *prob, IloNumVar ****inv_aeij);
void delete_inv_aei(Problem *prob, IloNumVar ***inv_aei);


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
    BaseMM(Problem *prob, TimeTracker *tt,
           unsigned long time_limit_sec, int numThreads,
           std::string logPath, std::string lpPath,
           std::string lp_algo, char vType): SolApprBase(prob, tt, time_limit_sec, numThreads, logPath) {
        if (lp_algo != "DSM" && lp_algo != "IPM") {
            assert(false);
        }
        this->lpPath = lpPath;
        this->pathPrefix = lpPath.substr(0, lpPath.find(".lp"));
        //
        y_ak = new_inv_ak(prob, env, vType);
        z_aek = new_inv_aek(prob, env, vType);
        x_aeij = new_inv_aeij(prob, env, vType);
        u_aei = new_inv_aei(prob, env);
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
        delete_inv_ak(prob, y_ak);
        delete_inv_aek(prob, z_aek);
        delete_inv_aeij(prob, x_aeij);
        delete_inv_aei(prob, u_aei);
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
    void def_preprocessing();
    void def_objF();
    void def_RUT_cnsts();
    void def_FC_cnsts_aeGiven(int a, int e);
    void def_AT_cnsts_aeGiven(int a, int e);
    void def_COM_cnsts();
};


#endif /* SolApprBase_h */
